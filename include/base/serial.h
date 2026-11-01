// =====================================================================
// serial.h -- turning values into bytes and back
// Used for network messages and for the state files each server keeps
// between runs (keys, preprocessing packages, database shares).
// Integers are big-endian and every variable-length item is length-prefixed,
// so a reader never has to guess where an item ends.
// =====================================================================
#pragma once
#include "base/common.h"
#include <emmintrin.h>
#include <cstring>
#include <cstdio>
#include <functional>
#include <thread>
#include <unistd.h>
#include <fstream>

namespace parry {

using Bytes = std::vector<uint8_t>;

struct Writer {
    Bytes buf;

    void u8(uint8_t v) { buf.push_back(v); }
    void u64(uint64_t v) {
        for (int i = 7; i >= 0; i--) buf.push_back(static_cast<uint8_t>(v >> (8 * i)));
    }
    void raw(const void* p, size_t n) {
        const uint8_t* q = static_cast<const uint8_t*>(p);
        buf.insert(buf.end(), q, q + n);
    }
    void bytes(const Bytes& b) { u64(b.size()); raw(b.data(), b.size()); }
    void str(const std::string& s) { u64(s.size()); raw(s.data(), s.size()); }
    void mpz(const mpz_class& x) {
        if (x < 0) throw std::invalid_argument("serial: negative integer");
        size_t n = (bitlen(x) + 7) / 8;
        if (x == 0) n = 0;
        u64(n);
        size_t old = buf.size();
        buf.resize(old + n);
        size_t got = 0;
        if (n) mpz_export(buf.data() + old, &got, 1, 1, 0, 0, x.get_mpz_t());
        if (got != n) throw std::logic_error("serial: mpz_export size");
    }
    void mpzs(const std::vector<mpz_class>& v) { u64(v.size()); for (const mpz_class& x : v) mpz(x); }
    void block(const __m128i& b) { raw(&b, 16); }
    void blocks(const std::vector<__m128i>& v) { u64(v.size()); for (const __m128i& b : v) block(b); }
    void u8s(const std::vector<uint8_t>& v) { bytes(v); }
};

struct Reader {
    const uint8_t* p;
    size_t n, pos = 0;
    explicit Reader(const Bytes& b) : p(b.data()), n(b.size()) {}
    Reader(const uint8_t* d, size_t len) : p(d), n(len) {}

    void need(size_t k) const {
        if (k > n - pos) throw std::runtime_error("serial: truncated input");
    }
    uint8_t u8() { need(1); return p[pos++]; }
    uint64_t u64() {
        need(8);
        uint64_t v = 0;
        for (int i = 0; i < 8; i++) v = (v << 8) | p[pos++];
        return v;
    }
    void raw(void* out, size_t k) { need(k); std::memcpy(out, p + pos, k); pos += k; }
    Bytes bytes() {
        uint64_t k = u64();
        need(k);
        Bytes b(p + pos, p + pos + k);
        pos += k;
        return b;
    }
    std::string str() {
        uint64_t k = u64();
        need(k);
        std::string s(reinterpret_cast<const char*>(p + pos), k);
        pos += k;
        return s;
    }
    mpz_class mpz() {
        uint64_t k = u64();
        need(k);
        mpz_class x = 0;
        if (k) mpz_import(x.get_mpz_t(), k, 1, 1, 0, 0, p + pos);
        pos += k;
        return x;
    }
    std::vector<mpz_class> mpzs() {
        uint64_t k = u64();
        if (k > n) throw std::runtime_error("serial: bad vector length");   // each item takes >= 8 bytes
        std::vector<mpz_class> v;
        v.reserve(k);
        for (uint64_t i = 0; i < k; i++) v.push_back(mpz());
        return v;
    }
    __m128i block() { __m128i b; raw(&b, 16); return b; }
    std::vector<__m128i> blocks() {
        uint64_t k = u64();
        need(k * 16);
        std::vector<__m128i> v(k);
        for (uint64_t i = 0; i < k; i++) v[i] = block();
        return v;
    }
    std::vector<uint8_t> u8s() { return bytes(); }
    bool done() const { return pos == n; }
};

inline void save_file(const std::string& path, const Bytes& b) {
    // unique temporary name: both servers may write the same public file into one folder
    std::string tmp = path + ".tmp" + std::to_string(::getpid()) + "_" + std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()));
    {
        std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("cannot write " + path);
        f.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
        if (!f) throw std::runtime_error("write failed: " + path);
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) throw std::runtime_error("cannot rename " + tmp);
}

inline Bytes load_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return Bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

inline bool file_exists(const std::string& path) { return static_cast<bool>(std::ifstream(path)); }

// Bytes <-> non-negative integers, so that byte strings (hashes, curve points, OT
// ciphertexts) can travel as log values. Fixed-length strings round-trip exactly.
inline mpz_class bytes_to_mpz(const uint8_t* p, size_t n) {
    mpz_class x = 0;
    mpz_import(x.get_mpz_t(), n, 1, 1, 0, 0, p);
    return x;
}
// Writes x into exactly n bytes (big-endian). False if x does not fit.
inline bool mpz_to_bytes(const mpz_class& x, uint8_t* out, size_t n) {
    if (x < 0 || bitlen(x) > 8 * n) return false;
    std::fill(out, out + n, 0);
    size_t got = 0;
    if (x == 0) return true;
    std::vector<uint8_t> tmp((bitlen(x) + 7) / 8);
    mpz_export(tmp.data(), &got, 1, 1, 0, 0, x.get_mpz_t());
    std::copy(tmp.begin(), tmp.begin() + static_cast<long>(got), out + (n - got));
    return true;
}

}  // namespace parry
