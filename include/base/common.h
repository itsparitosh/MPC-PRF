// =====================================================================
// common.h -- basic helpers used everywhere
//   * secure randomness from the operating system
//   * arithmetic in the prime field F_p (all shared values live here)
//   * Fault: the public evidence object used for blame
// =====================================================================
#pragma once
#include <gmpxx.h>
#include <sys/random.h>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <exception>
#include <future>
#include <thread>

namespace parry {

// Fills buf with n bytes from the OS CSPRNG (getrandom).
inline void random_bytes(unsigned char* buf, size_t n) {
    size_t got = 0;
    while (got < n) {
        ssize_t r = getrandom(buf + got, n - got, 0);
        if (r < 0) throw std::runtime_error("getrandom failed");
        got += static_cast<size_t>(r);
    }
}

// Uniform integer in [0, bound). We draw 128 extra random bits before reducing,
// so the bias of the final "mod bound" is below 2^-128.
inline mpz_class rand_below(const mpz_class& bound) {
    if (bound <= 0) throw std::invalid_argument("rand_below: bound must be positive");
    size_t nbytes = (mpz_sizeinbase(bound.get_mpz_t(), 2) + 7) / 8 + 16;
    std::vector<unsigned char> buf(nbytes);
    random_bytes(buf.data(), nbytes);
    mpz_class r;
    mpz_import(r.get_mpz_t(), nbytes, 1, 1, 0, 0, buf.data());
    return r % bound;
}

// Uniform integer with exactly `bits` bits (top bit set).
inline mpz_class rand_bits(unsigned bits) {
    mpz_class top = mpz_class(1) << (bits - 1);
    return top + rand_below(top);
}

// Number of bits of x.
inline unsigned bitlen(const mpz_class& x) { return static_cast<unsigned>(mpz_sizeinbase(x.get_mpz_t(), 2)); }

// Arithmetic modulo the prime p. Every result is reduced into [0, p).
struct Field {
    mpz_class p;
    mpz_class mod(const mpz_class& a) const { mpz_class r = a % p; if (r < 0) r += p; return r; }
    mpz_class add(const mpz_class& a, const mpz_class& b) const { return mod(a + b); }
    mpz_class sub(const mpz_class& a, const mpz_class& b) const { return mod(a - b); }
    mpz_class mul(const mpz_class& a, const mpz_class& b) const { return mod(a * b); }
    mpz_class neg(const mpz_class& a) const { return mod(-a); }
    mpz_class pow(const mpz_class& a, const mpz_class& e) const {
        mpz_class r, base = mod(a);
        mpz_powm(r.get_mpz_t(), base.get_mpz_t(), e.get_mpz_t(), p.get_mpz_t());
        return r;
    }
    mpz_class inv(const mpz_class& a) const {
        mpz_class r, base = mod(a);
        if (mpz_invert(r.get_mpz_t(), base.get_mpz_t(), p.get_mpz_t()) == 0)
            throw std::domain_error("Field::inv: zero has no inverse");
        return r;
    }
    mpz_class rand() const { return rand_below(p); }               // uniform in F_p
    mpz_class rand_nonzero() const { return rand_below(p - 1) + 1; } // uniform in F_p^*
};

// Blame(S_server): public evidence that server `server` deviated.
//   record >= 0 : index of a signed record on the public log that fails a public check.
//                 Anyone can re-verify the signature and the failed check.
//   record = -1 : the server refused to send a message it had to send.
struct Fault : std::runtime_error {
    int server;
    long record;
    Fault(int s, const std::string& reason, long rec)
        : std::runtime_error("S" + std::to_string(s) + ": " + reason), server(s), record(rec) {}
};

// ---- Multi-core helpers ----
// Number of worker threads: PARRY_THREADS if set, otherwise all hardware threads.
inline unsigned& num_threads() {
    static unsigned n = [] {
        if (const char* e = std::getenv("PARRY_THREADS")) {
            int v = std::atoi(e);
            if (v >= 1) return static_cast<unsigned>(v);
        }
        unsigned h = std::thread::hardware_concurrency();
        return h == 0 ? 1u : h;
    }();
    return n;
}

// Thread budget. PARRY_THREADS is a hard limit on the total number of busy threads, also
// when parallel parts are nested (e.g. many PRF calls in parallel, each with parallel proofs).
// Every thread has a budget; when work is split, the budget is split between the parts.
inline unsigned& thread_budget_slot() { thread_local unsigned b = 0; return b; }   // 0 = not set
inline unsigned thread_budget() { unsigned b = thread_budget_slot(); return b ? b : num_threads(); }

// Sets the calling thread's budget for the lifetime of this object.
struct BudgetScope {
    unsigned old;
    explicit BudgetScope(unsigned b) : old(thread_budget_slot()) { thread_budget_slot() = b < 1 ? 1 : b; }
    ~BudgetScope() { thread_budget_slot() = old; }
};

// Runs body(i) for i = 0..n-1 on up to (budget) threads; each worker gets an equal part
// of the budget for any parallel work nested inside. Exceptions (e.g. a Fault) thrown in a
// worker are re-thrown in the caller. With a budget of one it is a plain loop.
template <class Body>
inline void parallel_for(size_t n, Body body) {
    const unsigned T = thread_budget();
    const size_t t = std::min<size_t>(T, n);
    if (t <= 1) { for (size_t i = 0; i < n; i++) body(i); return; }
    std::vector<std::thread> pool;
    std::vector<std::exception_ptr> err(t);
    for (size_t k = 0; k < t; k++) {
        const unsigned share = static_cast<unsigned>(T / t + (k < T % t ? 1 : 0));
        pool.emplace_back([&, k, share] {
            BudgetScope scope(share);
            try { for (size_t i = k; i < n; i += t) body(i); } catch (...) { err[k] = std::current_exception(); }
        });
    }
    for (std::thread& th : pool) th.join();
    for (std::exception_ptr& e : err)
        if (e) std::rethrow_exception(e);
}

// Runs two independent tasks at the same time (if the budget allows) and returns both
// results; the budget is split between them. If both throw, the first exception is reported.
template <class F1, class F2>
inline auto run_both(F1 f1, F2 f2) {
    const unsigned T = thread_budget();
    if (T <= 1) {
        auto a = f1();
        auto b = f2();
        return std::make_pair(std::move(a), std::move(b));
    }
    const unsigned t2 = T / 2, t1 = T - t2;
    auto fut = std::async(std::launch::async, [&f2, t2] { BudgetScope scope(t2); return f2(); });
    std::exception_ptr e1;
    decltype(f1()) a{};
    {
        BudgetScope scope(t1);
        try { a = f1(); } catch (...) { e1 = std::current_exception(); }
    }
    auto b = fut.get();          // re-throws an exception of f2
    if (e1) std::rethrow_exception(e1);
    return std::make_pair(std::move(a), std::move(b));
}

}  // namespace parry
