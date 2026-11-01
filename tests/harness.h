// =====================================================================
// harness.h -- runs a two-server test inside one test binary
// S0 runs in the parent process and S1 in a forked child; they talk over a
// socketpair, with exactly the same Channel/Party code as the real servers.
// Each side returns a short outcome string; a Fault becomes "BLAME S<b>: ...".
// =====================================================================
#pragma once
#include "init/keys.h"
#include <sys/wait.h>
#include <functional>

namespace parry {

struct Outcome {
    std::string s[2];
    bool ok() const { return s[0] == "ok" && s[1] == "ok"; }
};

using PartyFn = std::function<std::string(Party&)>;

inline std::string run_side(int b, const Params& P, int fd, const std::string& dir, const PartyFn& f,
                            const std::set<std::string>& cheats, bool init) {
    Channel ch(fd);
    Party c(b, P);
    c.ch = &ch;
    c.adv.points = cheats;
    try {
        if (init) {
            ServerState s = run_init(c);
            save_state(dir, s);
            return "ok";
        }
        install(c, load_state(dir, b));
        c.start_session("test");
        return f(c);
    } catch (const Fault& e) {
        ch.send_abort(e.what());
        return std::string("BLAME S") + std::to_string(e.server) + ": " + e.what();
    } catch (const PeerGone& e) {
        return std::string("GONE: ") + e.what();
    } catch (const std::exception& e) {
        ch.send_abort(e.what());
        return std::string("ERR: ") + e.what();
    }
}

// cheat_server: which server deviates (-1: nobody); cheats: its deviation points.
inline Outcome two_party(const Params& P, const std::string& dir, const PartyFn& f, int cheat_server = -1,
                         const std::set<std::string>& cheats = {}, bool init = false) {
    int sv[2], pipefd[2];
    if (socketpair(AF_UNIX, SOCK_STREAM, 0, sv) != 0 || pipe(pipefd) != 0) throw std::runtime_error("socketpair/pipe");
    std::cout.flush();
    pid_t pid = fork();
    if (pid < 0) throw std::runtime_error("fork failed");
    if (pid == 0) {   // child: S1
        ::close(sv[0]);
        ::close(pipefd[0]);
        std::string r = run_side(1, P, sv[1], dir, f, cheat_server == 1 ? cheats : std::set<std::string>{}, init);
        ssize_t w = ::write(pipefd[1], r.data(), r.size());
        (void)w;
        ::close(pipefd[1]);
        _exit(0);
    }
    ::close(sv[1]);
    ::close(pipefd[1]);
    Outcome o;
    o.s[0] = run_side(0, P, sv[0], dir, f, cheat_server == 0 ? cheats : std::set<std::string>{}, init);
    std::string r;
    char buf[4096];
    ssize_t k;
    while ((k = ::read(pipefd[0], buf, sizeof buf)) > 0) r.append(buf, static_cast<size_t>(k));
    ::close(pipefd[0]);
    int st = 0;
    waitpid(pid, &st, 0);
    o.s[1] = r.empty() ? "child crashed" : r;
    return o;
}

}  // namespace parry
