// =====================================================================
// params.h -- public parameters (paper: Section "Notation", Parameters)
//   * PRF field F_p (Gold PRF): p = 2^lambda * d + 1, lambda = 128, |d| = 256 bits,
//     so |p| = 384 bits. F_k(x) = (k+x)^d is a 2^lambda-th root of unity,
//     i.e. one evaluation gives lambda = 128 output bits.
//   * Commitment group G: the subgroup of order exactly p inside Z_Q^*, where
//     Q = 2 * p * q_0 * ... * q_7 + 1 is a 3072-bit prime and every q_i is a prime
//     larger than 2^128 (needed for batch verification, see below). Generators g and h are derived by hashing,
//     so nobody knows log_g h (needed for Pedersen binding).
// Parameters are generated once, cached in params.txt, and re-validated on load.
// =====================================================================
#pragma once
#include "base/crypto.h"
#include <fstream>
#include <iostream>

namespace parry {

struct Params {
    unsigned lambda = 0;
    mpz_class d;                    // PRF exponent
    Field F;                        // arithmetic mod p
    mpz_class Q, g, h;              // commitment group and generators
    std::vector<mpz_class> qf;      // Q = 2 * p * qf[0] * ... * qf[k-1] + 1, every qf[i] prime > 2^128
    const mpz_class& p() const { return F.p; }
};

// Why the special form of Q: batch verification (sigma.h) checks many equations at once
// with random 128-bit weights. That is only sound if every element involved lies in a
// group whose prime factors are all larger than 2^128 (we require 2^256, see below). With Q - 1 = 2 * p * prod(qf),
// the quadratic residues mod Q form exactly such a group (order p * prod(qf)), and
// membership in it is a cheap Legendre-symbol test instead of an exponentiation.
constexpr unsigned kCofactorPrimes = 8;
// Each cofactor prime must be larger than any difference of two 256-bit Fiat-Shamir
// challenges: then a commitment that passes a proof equation for two challenges has no
// component outside G (special soundness), and the 128-bit batch weights are sound too.
constexpr unsigned kMinCofactorBits = 257;

inline bool is_probable_prime(const mpz_class& n) { return mpz_probab_prime_p(n.get_mpz_t(), 50) > 0; }

// Maps a seed to an element of the order-p subgroup. Used for g and h, so that nobody knows log_g h.
inline mpz_class hash_to_subgroup(const mpz_class& Q, const mpz_class& p, const std::string& seed) {
    mpz_class cof = (Q - 1) / p;
    for (int64_t ctr = 0;; ctr++) {
        Transcript t("hash-to-group");
        t.put_str(seed);
        t.put("ctr", ctr);
        mpz_class a = t.challenge_mod("x", Q - 3) + 2, y;
        mpz_powm(y.get_mpz_t(), a.get_mpz_t(), cof.get_mpz_t(), Q.get_mpz_t());
        if (y != 1) return y;
    }
}

inline mpz_class random_prime(unsigned bits) {
    mpz_class r, s = rand_bits(bits);
    mpz_nextprime(r.get_mpz_t(), s.get_mpz_t());
    return r;
}

inline Params generate_params(unsigned lambda, unsigned q_bits) {
    Params P;
    P.lambda = lambda;
    const mpz_class two_l = mpz_class(1) << lambda;
    do {                                    // p = 2^lambda * d + 1 with d odd, |d| = 2*lambda bits
        P.d = rand_bits(2 * lambda) | 1;
        P.F.p = two_l * P.d + 1;
    } while (!is_probable_prime(P.F.p));
    // Q = 2 * p * q_0 * ... * q_7 + 1 with |Q| = q_bits. The first seven cofactor primes are
    // fixed; the last one is redrawn until Q is prime and has the right length.
    const unsigned rest = q_bits - 1 - bitlen(P.F.p);   // bits left for the cofactor primes
    const unsigned each = rest / kCofactorPrimes;
    for (;;) {
        P.qf.clear();
        mpz_class prod = 2 * P.F.p;
        for (unsigned i = 0; i + 1 < kCofactorPrimes; i++) { P.qf.push_back(random_prime(each)); prod *= P.qf.back(); }
        const unsigned last = q_bits - bitlen(prod);
        bool found = false;
        for (int tries = 0; tries < 20000 && !found; tries++) {
            mpz_class ql = random_prime(last);
            mpz_class Q = prod * ql + 1;
            if (bitlen(Q) == q_bits && is_probable_prime(Q)) { P.qf.push_back(ql); P.Q = Q; found = true; }
        }
        if (found) break;
    }
    P.g = hash_to_subgroup(P.Q, P.F.p, "generator-g");
    P.h = hash_to_subgroup(P.Q, P.F.p, "generator-h");
    return P;
}

// Checks every property the protocols rely on.
inline void validate_params(const Params& P) {
    const mpz_class two_l = mpz_class(1) << P.lambda;
    mpz_class t;
    auto in_sub = [&](const mpz_class& x) {
        mpz_powm(t.get_mpz_t(), x.get_mpz_t(), P.p().get_mpz_t(), P.Q.get_mpz_t());
        return x > 1 && x < P.Q && t == 1;
    };
    mpz_class prod = 2 * P.p();
    for (const mpz_class& q : P.qf) {
        if (bitlen(q) < kMinCofactorBits || !is_probable_prime(q) || q == P.p())
            throw std::runtime_error("invalid parameters: bad cofactor prime");
        prod *= q;
    }
    if (P.qf.empty() || prod + 1 != P.Q || P.p() != two_l * P.d + 1 || !is_probable_prime(P.p()) ||
        !is_probable_prime(P.Q) || !in_sub(P.g) || !in_sub(P.h) || P.g == P.h ||
        P.h != hash_to_subgroup(P.Q, P.p(), "generator-h") || P.g != hash_to_subgroup(P.Q, P.p(), "generator-g"))
        throw std::runtime_error("invalid parameters");
}

inline void save_params(const Params& P, const std::string& file) {
    std::ofstream f(file);
    f << P.lambda << "\n" << P.d.get_str(16) << "\n" << P.F.p.get_str(16) << "\n"
      << P.Q.get_str(16) << "\n" << P.g.get_str(16) << "\n" << P.h.get_str(16) << "\n" << P.qf.size() << "\n";
    for (const mpz_class& q : P.qf) f << q.get_str(16) << "\n";
}

inline bool load_params(Params& P, const std::string& file) {
    std::ifstream f(file);
    std::string d, p, Q, g, h;
    size_t k = 0;
    if (!(f >> P.lambda >> d >> p >> Q >> g >> h >> k) || k == 0 || k > 64) return false;
    P.d.set_str(d, 16); P.F.p.set_str(p, 16); P.Q.set_str(Q, 16); P.g.set_str(g, 16); P.h.set_str(h, 16);
    P.qf.clear();
    for (size_t i = 0; i < k; i++) {
        std::string q;
        if (!(f >> q)) return false;
        P.qf.emplace_back(q, 16);
    }
    return true;
}

// Loads cached parameters, or generates them once (this takes some time).
inline Params get_params(const std::string& file, unsigned lambda = 128, unsigned q_bits = 3072) {
    Params P;
    bool ok = load_params(P, file) && P.lambda == lambda && bitlen(P.Q) == q_bits;
    if (ok) {
        try { validate_params(P); } catch (const std::exception&) { ok = false; }
    }
    if (!ok) {
        std::cout << "Generating parameters (one time, may take a minute)..." << std::endl;
        P = generate_params(lambda, q_bits);
        validate_params(P);
        save_params(P, file);
    }
    return P;
}

}  // namespace parry
