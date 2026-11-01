// =====================================================================
// paillier.h -- Paillier encryption (used by MtA in preprocessing)
//
//   Enc_N(m; eta) = (1 + N)^m * eta^N mod N^2      (additively homomorphic)
//
// N is a product of two safe primes. At init each server proves that its N is a
// Paillier-Blum modulus with no small factors (zk/cggmp.h: pi^mod, pi^fac). Then
// gcd(N, phi(N)) = 1, so every ciphertext has exactly one plaintext in Z_N, and
// pi^dec can pin down which integer a server decrypted.
//
// The key owner knows P and Q and uses the CRT for everything it computes
// modulo N or N^2: decryption, encryption under its own key, and the
// exponentiations it does when it verifies proofs made under its own key.
// The results are exactly the same numbers, only computed faster.
// =====================================================================
#pragma once
#include "base/crypto.h"

namespace parry {

struct PaillierPub {
    mpz_class N, N2;

    // Enc(m; r) = (1 + mN) r^N mod N^2, for 0 <= m < N and r in Z_N^*.
    mpz_class enc(const mpz_class& m, const mpz_class& r) const {
        if (m < 0 || m >= N) throw std::invalid_argument("Paillier: plaintext out of range");
        mpz_class rn;
        mpz_powm(rn.get_mpz_t(), r.get_mpz_t(), N.get_mpz_t(), N2.get_mpz_t());
        return (1 + m * N) % N2 * rn % N2;
    }
    mpz_class rand_r() const {
        mpz_class r;
        do { r = rand_below(N); } while (r == 0 || gcd(r, N) != 1);
        return r;
    }
    bool is_unit_mod_N2(const mpz_class& c) const { return c > 0 && c < N2 && gcd(c, N) == 1; }
    // Homomorphic operations.
    mpz_class add(const mpz_class& c1, const mpz_class& c2) const { return c1 * c2 % N2; }
    mpz_class scale(const mpz_class& c, const mpz_class& k) const {
        mpz_class r;
        mpz_powm(r.get_mpz_t(), c.get_mpz_t(), k.get_mpz_t(), N2.get_mpz_t());
        return r;
    }
};

class PaillierKey {
public:
    PaillierPub pub;

    // A modulus N = P Q with P, Q safe primes (P = 2P' + 1 with P' prime). Such an N is a
    // Paillier-Blum modulus (P, Q = 3 mod 4, gcd(N, phi(N)) = 1), which the proofs in zk/cggmp.h
    // need, and it can also serve as the server's ring-Pedersen modulus.
    static PaillierKey generate(unsigned modulus_bits) {
        const unsigned half = modulus_bits / 2;
        mpz_class P, Q;
        do {
            P = random_safe_prime(half);
            Q = random_safe_prime(half);
        } while (P == Q || bitlen(P * Q) != modulus_bits);
        return from_primes(P, Q);
    }

    // Safe prime with exactly `bits` bits and its top two bits set (so that the product of
    // two of them has exactly 2*bits bits). Candidates are sieved by small primes first, for
    // both P' and 2P' + 1, which makes the search take seconds instead of minutes.
    static mpz_class random_safe_prime(unsigned bits) {
        static const std::vector<unsigned> small = [] {
            const unsigned B = 1u << 16;
            std::vector<bool> comp(B, false);
            std::vector<unsigned> v;
            for (unsigned i = 3; i < B; i += 2) {
                if (comp[i]) continue;
                v.push_back(i);
                for (unsigned long j = static_cast<unsigned long>(i) * i; j < B; j += 2 * i) comp[j] = true;
            }
            return v;
        }();
        for (;;) {
            mpz_class q = rand_bits(bits - 1);                  // P' has bits-1 bits, top bit set
            mpz_setbit(q.get_mpz_t(), bits - 3);                // second bit, so P has its top two bits set
            mpz_setbit(q.get_mpz_t(), 0);                       // odd
            std::vector<unsigned> res(small.size());
            for (size_t i = 0; i < small.size(); i++) res[i] = static_cast<unsigned>(mpz_fdiv_ui(q.get_mpz_t(), small[i]));
            for (unsigned long delta = 0; delta < (1ul << 20); delta += 2) {
                bool ok = true;
                for (size_t i = 0; i < small.size() && ok; i++) {
                    const unsigned m = small[i];
                    const unsigned long r = (res[i] + delta) % m;
                    ok = r != 0 && (2 * r + 1) % m != 0;            // m divides neither P' nor 2P' + 1
                }
                if (!ok) continue;
                mpz_class qq = q + delta, PP = 2 * qq + 1;
                if (bitlen(qq) != bits - 1) break;
                if (mpz_probab_prime_p(qq.get_mpz_t(), 1) && mpz_probab_prime_p(PP.get_mpz_t(), 1) &&
                    mpz_probab_prime_p(qq.get_mpz_t(), 40) && mpz_probab_prime_p(PP.get_mpz_t(), 40))
                    return PP;
            }
        }
    }

    // Decryption read as a signed integer in (-N/2, N/2). The range proofs only bound values
    // in absolute value, so a (dishonest) plaintext may be negative.
    mpz_class dec_signed(const mpz_class& c) const {
        mpz_class t = dec(c);
        if (2 * t > pub.N) t -= pub.N;
        return t;
    }

    // Rebuilds the key from its two primes (used when a server reloads its state file).
    static PaillierKey from_primes(const mpz_class& P, const mpz_class& Q) {
        PaillierKey k;
        k.P_ = P; k.Q_ = Q;
        k.pub.N = P * Q;
        k.phi_ = (P - 1) * (Q - 1);
        if (P == Q || gcd(k.pub.N, k.phi_) != 1) throw std::invalid_argument("Paillier: bad primes");
        k.pub.N2 = k.pub.N * k.pub.N;
        mpz_invert(k.mu_.get_mpz_t(), k.phi_.get_mpz_t(), k.pub.N.get_mpz_t());     // phi^{-1} mod N
        mpz_invert(k.ninv_.get_mpz_t(), k.pub.N.get_mpz_t(), k.phi_.get_mpz_t());   // N^{-1} mod phi
        k.init_crt();
        return k;
    }
    const mpz_class& prime_P() const { return P_; }
    const mpz_class& prime_Q() const { return Q_; }

    // ---- CRT speed-up: the key owner knows P and Q, so it works mod P^2 and Q^2 ----
    // (numbers half as long, so each exponentiation is about 4x cheaper)

    // Dec(c), computed separately mod P and mod Q:
    //   m_P = L_P(c^(P-1) mod P^2) * hP mod P,  L_P(u) = (u - 1) / P,  and likewise mod Q.
    mpz_class dec(const mpz_class& c) const {
        mpz_class up, uq, cp = c % P2_, cq = c % Q2_;
        mpz_powm(up.get_mpz_t(), cp.get_mpz_t(), Pm1_.get_mpz_t(), P2_.get_mpz_t());
        mpz_powm(uq.get_mpz_t(), cq.get_mpz_t(), Qm1_.get_mpz_t(), Q2_.get_mpz_t());
        mpz_class mp = (up - 1) / P_ * hP_ % P_;
        mpz_class mq = (uq - 1) / Q_ * hQ_ % Q_;
        return crt(mp, mq, P_, Q_, Pinv_Q_);
    }

    // Dec without CRT (reference, used in the tests to check the CRT version).
    mpz_class dec_plain(const mpz_class& c) const {
        mpz_class u;
        mpz_powm(u.get_mpz_t(), c.get_mpz_t(), phi_.get_mpz_t(), pub.N2.get_mpz_t());
        return (u - 1) / pub.N * mu_ % pub.N;
    }

    // Enc(m; r) under the owner's own key: r^N mod N^2 via CRT mod P^2 and Q^2.
    // Gives exactly the same ciphertext as pub.enc(m, r).
    mpz_class enc(const mpz_class& m, const mpz_class& r) const {
        if (m < 0 || m >= pub.N) throw std::invalid_argument("Paillier: plaintext out of range");
        mpz_class rp, rq, bp = r % P2_, bq = r % Q2_;
        mpz_powm(rp.get_mpz_t(), bp.get_mpz_t(), eP_.get_mpz_t(), P2_.get_mpz_t());   // r^N mod P^2
        mpz_powm(rq.get_mpz_t(), bq.get_mpz_t(), eQ_.get_mpz_t(), Q2_.get_mpz_t());   // r^N mod Q^2
        mpz_class rn = crt(rp, rq, P2_, Q2_, P2inv_Q2_);
        return (1 + m * pub.N) % pub.N2 * rn % pub.N2;
    }

    // The unique r with c = Enc(m; r). Throws if c does not encrypt m.
    mpz_class randomness(const mpz_class& c, const mpz_class& m) const {
        mpz_class cr = c * ((1 - m * pub.N) % pub.N2 + pub.N2) % pub.N2;   // c * (1+N)^{-m} = r^N mod N^2
        mpz_class base = cr % pub.N, r;
        mpz_powm(r.get_mpz_t(), base.get_mpz_t(), ninv_.get_mpz_t(), pub.N.get_mpz_t());
        if (pub.enc(m, r) != c) throw std::invalid_argument("Paillier: ciphertext does not encrypt m");
        return r;
    }

    // base^e mod N and base^e mod N^2 via the CRT. base must be a unit; e may be negative.
    // The exponent is reduced mod the order of the group mod P (or P^2) and mod Q (or Q^2).
    mpz_class pow_mod_n(const mpz_class& base, const mpz_class& e) const {
        return crt(pow_reduced(base, e, P_, Pm1_), pow_reduced(base, e, Q_, Qm1_), P_, Q_, Pinv_Q_);
    }
    mpz_class pow_mod_n2(const mpz_class& base, const mpz_class& e) const {
        return crt(pow_reduced(base, e, P2_, P_ * Pm1_), pow_reduced(base, e, Q2_, Q_ * Qm1_), P2_, Q2_, P2inv_Q2_);
    }

    // N-th root mod N (exists for every unit because gcd(N, phi(N)) = 1).
    mpz_class nth_root(const mpz_class& x) const {
        mpz_class r;
        mpz_powm(r.get_mpz_t(), x.get_mpz_t(), ninv_.get_mpz_t(), pub.N.get_mpz_t());
        return r;
    }

private:
    mpz_class P_, Q_, phi_, mu_, ninv_;
    mpz_class P2_, Q2_, Pm1_, Qm1_, hP_, hQ_, Pinv_Q_, P2inv_Q2_, eP_, eQ_;

    // base^e mod m, where `order` is a multiple of the order of base mod m.
    static mpz_class pow_reduced(const mpz_class& base, const mpz_class& e, const mpz_class& m, const mpz_class& order) {
        mpz_class b = base % m, ee = e % order, r;
        if (ee < 0) ee += order;
        mpz_powm(r.get_mpz_t(), b.get_mpz_t(), ee.get_mpz_t(), m.get_mpz_t());
        return r;
    }

    // x = a mod m1, x = b mod m2  ->  x mod m1*m2   (inv = m1^{-1} mod m2)
    static mpz_class crt(const mpz_class& a, const mpz_class& b, const mpz_class& m1, const mpz_class& m2,
                         const mpz_class& inv) {
        mpz_class t = (b - a) % m2;
        if (t < 0) t += m2;
        return a + m1 * (t * inv % m2);
    }

    void init_crt() {
        P2_ = P_ * P_; Q2_ = Q_ * Q_; Pm1_ = P_ - 1; Qm1_ = Q_ - 1;
        const mpz_class g = 1 + pub.N;
        mpz_class u;
        mpz_class gp = g % P2_, gq = g % Q2_;
        mpz_powm(u.get_mpz_t(), gp.get_mpz_t(), Pm1_.get_mpz_t(), P2_.get_mpz_t());
        hP_ = (u - 1) / P_;
        mpz_invert(hP_.get_mpz_t(), hP_.get_mpz_t(), P_.get_mpz_t());
        mpz_powm(u.get_mpz_t(), gq.get_mpz_t(), Qm1_.get_mpz_t(), Q2_.get_mpz_t());
        hQ_ = (u - 1) / Q_;
        mpz_invert(hQ_.get_mpz_t(), hQ_.get_mpz_t(), Q_.get_mpz_t());
        mpz_invert(Pinv_Q_.get_mpz_t(), P_.get_mpz_t(), Q_.get_mpz_t());
        mpz_invert(P2inv_Q2_.get_mpz_t(), P2_.get_mpz_t(), Q2_.get_mpz_t());
        eP_ = pub.N % (P_ * Pm1_);   // exponent N reduced mod phi(P^2) = P(P-1)
        eQ_ = pub.N % (Q_ * Qm1_);
    }
};

}  // namespace parry
