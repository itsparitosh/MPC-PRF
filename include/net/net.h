// =====================================================================
// net.h -- the TCP link between S0 and S1
//
// Same layout as Remise (network.hpp): S0 listens on a port (9200 by
// default), S1 connects to it and retries until S0 is up, so the two
// servers can be started in two terminals in any order.
//
// Unlike Remise we use plain blocking sockets instead of asio coroutines:
// our protocol code is synchronous and spreads its heavy crypto over many
// threads, which does not mix well with coroutines. All network calls are
// made from the main thread between two computation steps.
//
// Every message is one frame: [type: 1 byte][length: 8 bytes][payload].
// type DATA carries protocol data. type ABORT tells the peer that we stop,
// together with the reason (e.g. which record we blame it for), so the peer
// does not mistake our abort for a crash.
// =====================================================================
#pragma once
#include "base/serial.h"
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <unistd.h>
#include <cerrno>
#include <chrono>
#include <iostream>
#include <thread>

namespace parry {

// The peer closed the link. `reason` is its abort message, or empty if it
// just disappeared (crash or refusal to continue).
struct PeerGone : std::runtime_error {
    std::string reason;
    bool deadline = false;    // the peer sent nothing before the round deadline
    explicit PeerGone(const std::string& r, bool late = false)
        : std::runtime_error(late ? "no message from the peer before the round deadline"
                                  : r.empty() ? "peer closed the connection" : "peer aborted: " + r),
          reason(r), deadline(late) {}
};

struct CommStats {
    uint64_t bytes_sent = 0, bytes_recv = 0, messages = 0;
    uint64_t total() const { return bytes_sent + bytes_recv; }
};

class Channel {
public:
    static constexpr uint8_t kData = 1, kAbort = 0xFF;
    static constexpr uint64_t kMaxFrame = 1ull << 34;   // 16 GiB, only to reject garbage lengths

    Channel() = default;
    explicit Channel(int fd) : fd_(fd) { set_nodelay(); }
    Channel(const Channel&) = delete;
    Channel& operator=(const Channel&) = delete;
    Channel(Channel&& o) noexcept : stats(o.stats), fd_(o.fd_) { o.fd_ = -1; }
    Channel& operator=(Channel&& o) noexcept {
        if (this != &o) { close(); fd_ = o.fd_; stats = o.stats; o.fd_ = -1; }
        return *this;
    }
    ~Channel() { close(); }

    // S0: wait for S1 on `port` (like Remise make_server).
    static Channel listen_on(uint16_t port) {
        int s = ::socket(AF_INET, SOCK_STREAM, 0);
        if (s < 0) throw std::runtime_error("socket() failed");
        int one = 1;
        ::setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
        sockaddr_in a{};
        a.sin_family = AF_INET;
        a.sin_addr.s_addr = htonl(INADDR_ANY);
        a.sin_port = htons(port);
        if (::bind(s, reinterpret_cast<sockaddr*>(&a), sizeof(a)) != 0 || ::listen(s, 1) != 0) {
            ::close(s);
            throw std::runtime_error("cannot listen on port " + std::to_string(port));
        }
        int c = ::accept(s, nullptr, nullptr);
        ::close(s);
        if (c < 0) throw std::runtime_error("accept() failed");
        return Channel(c);
    }

    // S1: connect to S0, retrying every 500 ms (like Remise connect_with_retry).
    static Channel connect_to(const std::string& host, uint16_t port, int max_tries = 600) {
        for (int attempt = 1;; attempt++) {
            addrinfo hints{}, *res = nullptr;
            hints.ai_family = AF_INET;
            hints.ai_socktype = SOCK_STREAM;
            if (::getaddrinfo(host.c_str(), std::to_string(port).c_str(), &hints, &res) == 0) {
                int s = ::socket(res->ai_family, res->ai_socktype, res->ai_protocol);
                if (s >= 0 && ::connect(s, res->ai_addr, res->ai_addrlen) == 0) {
                    ::freeaddrinfo(res);
                    return Channel(s);
                }
                if (s >= 0) ::close(s);
                ::freeaddrinfo(res);
            }
            if (attempt >= max_tries) throw std::runtime_error("cannot connect to " + host + ":" + std::to_string(port));
            if (attempt == 1) std::cerr << "[S1] waiting for S0 at " << host << ":" << port << " ..." << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
    }

    void send(const Bytes& payload) { send_frame(kData, payload); }

    Bytes recv() {
        uint8_t hdr[9];
        read_all(hdr, 9);
        uint64_t len = 0;
        for (int i = 1; i < 9; i++) len = (len << 8) | hdr[i];
        if (len > kMaxFrame) throw PeerGone("frame too large");
        Bytes b(len);
        if (len) read_all(b.data(), len);
        stats.bytes_recv += 9 + len;
        if (hdr[0] == kAbort) throw PeerGone(std::string(b.begin(), b.end()));
        if (hdr[0] != kData) throw PeerGone("unknown frame type");
        return b;
    }

    // Both sides send and receive at the same time. The send runs in a helper
    // thread, so two large messages cannot block each other in the TCP buffers.
    Bytes exchange(const Bytes& mine) {
        std::exception_ptr err;
        std::thread t([&] { try { send(mine); } catch (...) { err = std::current_exception(); } });
        Bytes theirs;
        try {
            theirs = recv();
        } catch (...) {
            t.join();
            throw;
        }
        t.join();
        if (err) std::rethrow_exception(err);
        stats.messages++;
        return theirs;
    }

    // Best effort: tell the peer why we stop. Never throws.
    void send_abort(const std::string& why) noexcept {
        try { send_frame(kAbort, Bytes(why.begin(), why.end())); } catch (...) {}
    }

    // Round deadline (paper: "a message that is not on the log by its deadline counts as
    // missing"). If the peer sends nothing for `seconds`, recv() stops with PeerGone, and
    // the honest server treats it as a missing message. 0 = wait forever.
    void set_deadline(unsigned seconds) {
        timeval tv{};
        tv.tv_sec = static_cast<time_t>(seconds);
        ::setsockopt(fd_, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    }

    void close() {
        if (fd_ >= 0) { ::shutdown(fd_, SHUT_RDWR); ::close(fd_); fd_ = -1; }
    }
    bool open() const { return fd_ >= 0; }

    CommStats stats;

private:
    int fd_ = -1;

    void set_nodelay() {
        int one = 1;
        ::setsockopt(fd_, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));   // fails harmlessly on socketpairs
    }

    void send_frame(uint8_t type, const Bytes& payload) {
        uint8_t hdr[9];
        hdr[0] = type;
        uint64_t len = payload.size();
        for (int i = 8; i >= 1; i--) { hdr[i] = static_cast<uint8_t>(len); len >>= 8; }
        write_all(hdr, 9);
        if (!payload.empty()) write_all(payload.data(), payload.size());
        stats.bytes_sent += 9 + payload.size();
    }

    void write_all(const uint8_t* p, size_t n) {
        while (n) {
            ssize_t w = ::send(fd_, p, n, MSG_NOSIGNAL);
            if (w < 0 && errno == EINTR) continue;
            if (w <= 0) throw PeerGone("");
            p += w;
            n -= static_cast<size_t>(w);
        }
    }
    void read_all(uint8_t* p, size_t n) {
        while (n) {
            ssize_t r = ::recv(fd_, p, n, 0);
            if (r < 0 && errno == EINTR) continue;
            if (r < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) throw PeerGone("", true);
            if (r <= 0) throw PeerGone("");
            p += r;
            n -= static_cast<size_t>(r);
        }
    }
};

}  // namespace parry
