// =====================================================================
// pedersen.h -- Pedersen commitments Com(m; rho) = g^m h^rho in G
//   * hiding: the commitment reveals nothing about m (perfectly)
//   * binding: nobody can open it to two different values (discrete log)
//   * homomorphic: Com(a; r) * Com(b; s) = Com(a + b; r + s), so commitments
//     to sums and scalar multiples can be derived by anyone.
// g^e and h^e use precomputed tables (8-bit windows) for speed.
// =====================================================================
#pragma once
#include "base/params.h"
#include <algorithm>

namespace parry {

// Fixed-base exponentiation: table[i][j] = base^(j * 256^i).
class FixedBase {
public:
    FixedBase(const mpz_class& base, const mpz_class& mod, unsigned max_exp_bits) : mod_(mod) {
        windows_ = (max_exp_bits + 7) / 8;
        table_.assign(windows_, std::vector<mpz_class>(256));
        mpz_class b = base;
        for (unsigned i = 0; i < windows_; i++) {
            table_[i][0] = 1;
            table_[i][1] = b;
            for (unsigned j = 2; j < 256; j++) table_[i][j] = table_[i][j - 1] * b % mod_;
            b = table_[i][255] * b % mod_;   // base^(256^(i+1))
        }
    }
    // Requires 0 <= e < 2^(8*windows).
    mpz_class pow(const mpz_class& e) const {
        if (e < 0 || bitlen(e) > 8 * windows_) throw std::invalid_argument("FixedBase: exponent out of range");
        std::vector<unsigned char> bytes(windows_, 0);
        size_t n = 0;
        if (e != 0) mpz_export(bytes.data(), &n, -1, 1, 0, 0, e.get_mpz_t());   // little-endian bytes
        mpz_class r = 1;
        for (size_t i = 0; i < n; i++)
            if (bytes[i]) r = r * table_[i][bytes[i]] % mod_;
        return r;
    }

private:
    mpz_class mod_;
    unsigned windows_;
    std::vector<std::vector<mpz_class>> table_;
};

class Pedersen {
public:
    explicit Pedersen(const Params& P)
        : P_(P), g_(P.g, P.Q, bitlen(P.p())), h_(P.h, P.Q, bitlen(P.p())) {}

    // Exponents are taken mod p (the group order).
    mpz_class g_pow(const mpz_class& e) const { return g_.pow(P_.F.mod(e)); }
    mpz_class h_pow(const mpz_class& e) const { return h_.pow(P_.F.mod(e)); }
    mpz_class pow(const mpz_class& a, const mpz_class& e) const {
        mpz_class r, ee = P_.F.mod(e);
        mpz_powm(r.get_mpz_t(), a.get_mpz_t(), ee.get_mpz_t(), P_.Q.get_mpz_t());
        return r;
    }
    mpz_class mul(const mpz_class& a, const mpz_class& b) const { return a * b % P_.Q; }
    mpz_class inv(const mpz_class& a) const {
        mpz_class r;
        if (mpz_invert(r.get_mpz_t(), a.get_mpz_t(), P_.Q.get_mpz_t()) == 0) throw std::domain_error("not invertible");
        return r;
    }

    mpz_class commit(const mpz_class& m, const mpz_class& r) const { return mul(g_pow(m), h_pow(r)); }
    bool opens(const mpz_class& C, const mpz_class& m, const mpz_class& r) const { return commit(m, r) == C; }

    // True iff x is an element of the order-p subgroup. Every received group element is checked.
    bool in_group(const mpz_class& x) const {
        if (x <= 0 || x >= P_.Q) return false;
        mpz_class t;
        mpz_powm(t.get_mpz_t(), x.get_mpz_t(), P_.p().get_mpz_t(), P_.Q.get_mpz_t());
        return t == 1;
    }

    // True iff x is a quadratic residue mod Q (Legendre symbol, as cheap as a gcd).
    // The quadratic residues form the group of order p * q_0 * ... * q_7, whose prime
    // factors are all > 2^128; that is all the batch verification needs.
    bool is_qr(const mpz_class& x) const {
        return x > 0 && x < P_.Q && mpz_jacobi(x.get_mpz_t(), P_.Q.get_mpz_t()) == 1;
    }

    // Are all elements of xs in G? Cheaper than one exponentiation per element:
    // every element must be a quadratic residue (Legendre test), and then
    // (prod_i x_i^{w_i})^p must be 1 for random 128-bit weights w_i. An element
    // outside G has a component of order q_j > 2^128, and that component survives
    // the random product except with probability 2^-128.
    bool all_in_group(const std::vector<mpz_class>& xs) const {
        if (xs.size() <= 2) {
            for (const mpz_class& x : xs)
                if (!in_group(x)) return false;
            return true;
        }
        for (const mpz_class& x : xs)
            if (!is_qr(x)) return false;
        const mpz_class two128 = mpz_class(1) << 128;
        std::vector<mpz_class> w(xs.size());
        for (mpz_class& wi : w) wi = rand_below(two128) + 1;
        mpz_class prod = multi_pow(xs, w), t;
        mpz_powm(t.get_mpz_t(), prod.get_mpz_t(), P_.p().get_mpz_t(), P_.Q.get_mpz_t());
        return t == 1;
    }

    // a^x * b^y mod Q in one pass (Shamir's trick, 2-bit joint windows): the squarings are
    // shared, so this costs about one exponentiation plus a quarter instead of two.
    mpz_class pow2(const mpz_class& a, const mpz_class& x, const mpz_class& b, const mpz_class& y) const {
        const mpz_class& Q = P_.Q;
        mpz_class tab[16];                       // tab[i + 4j] = a^i * b^j, i, j in 0..3
        tab[0] = 1; tab[1] = a % Q; tab[4] = b % Q;
        tab[2] = tab[1] * tab[1] % Q; tab[3] = tab[2] * tab[1] % Q;
        tab[8] = tab[4] * tab[4] % Q; tab[12] = tab[8] * tab[4] % Q;
        for (int j = 1; j < 4; j++)
            for (int i = 1; i < 4; i++) tab[i + 4 * j] = tab[i] * tab[4 * j] % Q;
        unsigned bits = std::max(bitlen(x), bitlen(y));
        bits += bits % 2;
        mpz_class r = 1;
        for (long k = static_cast<long>(bits) - 2; k >= 0; k -= 2) {
            r = r * r % Q;
            r = r * r % Q;
            unsigned i = static_cast<unsigned>(mpz_tstbit(x.get_mpz_t(), k) | (mpz_tstbit(x.get_mpz_t(), k + 1) << 1));
            unsigned j = static_cast<unsigned>(mpz_tstbit(y.get_mpz_t(), k) | (mpz_tstbit(y.get_mpz_t(), k + 1) << 1));
            if (i | j) r = r * tab[i + 4 * j] % Q;
        }
        return r;
    }

    // Multi-exponentiation prod_i bases[i]^exps[i] mod Q. With several cores, the bases are
    // split into one chunk per thread, each chunk is done with Pippenger, and the partial
    // results are multiplied.
    mpz_class multi_pow(const std::vector<mpz_class>& bases, const std::vector<mpz_class>& exps) const {
        if (bases.size() != exps.size()) throw std::invalid_argument("multi_pow: size mismatch");
        const size_t t = std::min<size_t>(thread_budget(), bases.size() / 64);
        if (t <= 1) return multi_pow_seq(bases, exps);
        std::vector<mpz_class> part(t);
        parallel_for(t, [&](size_t k) {
            std::vector<mpz_class> b, e;
            for (size_t i = k; i < bases.size(); i += t) { b.push_back(bases[i]); e.push_back(exps[i]); }
            part[k] = multi_pow_seq(b, e);
        });
        mpz_class r = 1;
        for (const mpz_class& x : part) r = r * x % P_.Q;
        return r;
    }

    // Multi-exponentiation on one core (Pippenger's bucket method).
    // Exponents must be >= 0. Much cheaper than separate exponentiations because all
    // bases share the same squarings.
    mpz_class multi_pow_seq(const std::vector<mpz_class>& bases, const std::vector<mpz_class>& exps) const {
        if (bases.size() != exps.size()) throw std::invalid_argument("multi_pow: size mismatch");
        const mpz_class& Q = P_.Q;
        size_t n = bases.size();
        if (n == 0) return 1;
        unsigned maxbits = 1;
        for (const mpz_class& e : exps) {
            if (e < 0) throw std::invalid_argument("multi_pow: negative exponent");
            maxbits = std::max(maxbits, bitlen(e));
        }
        // window size c: about log2(n) - 3, between 2 and 10 (c = 7 for ~1000 bases, tuned)
        unsigned c = 2;
        while (c < 10 && (1u << (c + 4)) < n) c++;
        const unsigned windows = (maxbits + c - 1) / c;
        mpz_class result = 1;
        std::vector<mpz_class> bucket(static_cast<size_t>(1) << c);
        for (long w = static_cast<long>(windows) - 1; w >= 0; w--) {
            for (unsigned i = 0; i < c; i++) result = result * result % Q;     // result^(2^c)
            for (auto& bk : bucket) bk = 1;
            bool any = false;
            for (size_t i = 0; i < n; i++) {
                unsigned digit = 0;
                for (unsigned j = 0; j < c; j++)
                    if (mpz_tstbit(exps[i].get_mpz_t(), static_cast<mp_bitcnt_t>(w * c + j))) digit |= 1u << j;
                if (digit) { bucket[digit] = bucket[digit] * bases[i] % Q; any = true; }
            }
            if (!any) continue;
            // sum_{k} k * bucket[k] (in the exponent) via running products
            mpz_class running = 1, acc = 1;
            for (size_t k = bucket.size() - 1; k >= 1; k--) {
                running = running * bucket[k] % Q;
                acc = acc * running % Q;
            }
            result = result * acc % Q;
        }
        return result;
    }

private:
    const Params& P_;
    FixedBase g_, h_;
};

}  // namespace parry
