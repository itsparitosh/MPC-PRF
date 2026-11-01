// =====================================================================
// sigma.h -- zero-knowledge proofs over Pedersen commitments
// (Sigma protocols made non-interactive with Fiat-Shamir; 256-bit challenges)
//
//   OpeningProof : "I know (x, r) with C = g^x h^r"                   (Okamoto)
//                  used for the PRF key shares and the random values zeta.
//   ProductProof : "C commits to a*b", where A commits to a and B to b.
//                  Proves A = g^a h^ra and C = B^a h^s with the same a.
//                  Used for pi^prod and for "n * t = 1" / "u * u^-1 = 1"
//                  (take C = g = Com(1; 0)).
//   PowerProof   : pi^pow, "Cy commits to x^e" for a public e >= 2. Built as a
//                  square-and-multiply chain: commit to each intermediate power
//                  and prove each step with a ProductProof.
//
// Each proof binds the prover's identity and a context string, so a proof
// cannot be copied by the other server or reused in another protocol step.
// The caller must check that statement elements are in G.
//
// Written for this project; these are the textbook Sigma protocols (Schnorr /
// Okamoto for openings, Cramer-Damgard CRYPTO '98 for products). The batch
// check of a power proof is the small-exponent test of Bellare-Garay-Rabin
// (EUROCRYPT '98); it is sound here because of the special form of Q (params.h).
// =====================================================================
#pragma once
#include "net/party.h"
#include <map>

namespace parry {

constexpr unsigned kChallengeBits = 256;

// ---------------- Opening proof ----------------
struct OpeningProof { mpz_class T, z1, z2; };

inline OpeningProof prove_opening(const Party& c, const std::string& ctx, int prover, const mpz_class& C,
                                  const mpz_class& x, const mpz_class& r) {
    mpz_class u = c.F().rand(), v = c.F().rand();
    OpeningProof pf;
    pf.T = c.ped.commit(u, v);
    Transcript t("opening");
    t.put_str(ctx); t.put("prover", prover); t.put("C", C); t.put("T", pf.T);
    mpz_class e = t.challenge_bits("e", kChallengeBits);
    pf.z1 = c.F().add(u, e * x);
    pf.z2 = c.F().add(v, e * r);
    return pf;
}

// Responses must be canonical (in [0, p)), so every valid proof has one encoding.
// The first messages T need no separate membership test: the verification equation
// T = g^z1 h^z2 * C^-e forces T into G whenever the statement C is in G.
// (Statement elements are checked by the caller; intermediate commitments of a
// power proof are checked in verify_power.)
inline bool canonical(const Party& c, const mpz_class& z) { return z >= 0 && z < c.P.p(); }

inline bool verify_opening(const Party& c, const std::string& ctx, int prover, const mpz_class& C, const OpeningProof& pf) {
    if (pf.T <= 0 || pf.T >= c.P.Q || !canonical(c, pf.z1) || !canonical(c, pf.z2)) return false;
    Transcript t("opening");
    t.put_str(ctx); t.put("prover", prover); t.put("C", C); t.put("T", pf.T);
    mpz_class e = t.challenge_bits("e", kChallengeBits);
    return c.ped.commit(pf.z1, pf.z2) == c.ped.mul(pf.T, c.ped.pow(C, e));
}

// ---------------- Product proof ----------------
struct ProductProof { mpz_class T1, T2, za, zr, zs; };

// Witness: a, ra with A = g^a h^ra, and s with C = B^a h^s.
inline ProductProof prove_product(const Party& c, const std::string& ctx, int prover, const mpz_class& A,
                                  const mpz_class& B, const mpz_class& C, const mpz_class& a,
                                  const mpz_class& ra, const mpz_class& s) {
    mpz_class u = c.F().rand(), v = c.F().rand(), w = c.F().rand();
    ProductProof pf;
    pf.T1 = c.ped.commit(u, v);
    pf.T2 = c.ped.mul(c.ped.pow(B, u), c.ped.h_pow(w));
    Transcript t("product");
    t.put_str(ctx); t.put("prover", prover);
    t.put("A", A); t.put("B", B); t.put("C", C); t.put("T1", pf.T1); t.put("T2", pf.T2);
    mpz_class e = t.challenge_bits("e", kChallengeBits);
    pf.za = c.F().add(u, e * a);
    pf.zr = c.F().add(v, e * ra);
    pf.zs = c.F().add(w, e * s);
    return pf;
}

inline mpz_class product_challenge(const std::string& ctx, int prover, const mpz_class& A, const mpz_class& B,
                                   const mpz_class& C, const mpz_class& T1, const mpz_class& T2) {
    Transcript t("product");
    t.put_str(ctx); t.put("prover", prover);
    t.put("A", A); t.put("B", B); t.put("C", C); t.put("T1", T1); t.put("T2", T2);
    return t.challenge_bits("e", kChallengeBits);
}

// Checks  g^za h^zr = T1 * A^e  and  B^za h^zs = T2 * C^e.
// The second check is computed as B^za * C^(p-e) * h^zs in one pass (Shamir's trick);
// this uses that B and C are in G (checked by the caller), so exponents work mod p.
inline bool verify_product(const Party& c, const std::string& ctx, int prover, const mpz_class& A,
                           const mpz_class& B, const mpz_class& C, const ProductProof& pf) {
    if (pf.T1 <= 0 || pf.T1 >= c.P.Q || pf.T2 <= 0 || pf.T2 >= c.P.Q || !canonical(c, pf.za) ||
        !canonical(c, pf.zr) || !canonical(c, pf.zs))
        return false;
    mpz_class e = product_challenge(ctx, prover, A, B, C, pf.T1, pf.T2);
    bool ok1 = c.ped.commit(pf.za, pf.zr) == c.ped.mul(pf.T1, c.ped.pow(A, e));
    bool ok2 = c.ped.mul(c.ped.pow2(B, pf.za, C, c.F().neg(e)), c.ped.h_pow(pf.zs)) == pf.T2;
    return ok1 && ok2;
}

// ---------------- Power proof ----------------
struct PowerProof {
    std::vector<mpz_class> inter;      // commitments to intermediate powers
    std::vector<ProductProof> steps;   // one product proof per square / multiply
};

// Square-and-multiply schedule for e: true = square, false = multiply by x.
inline std::vector<bool> power_schedule(const mpz_class& e) {
    if (e < 2) throw std::invalid_argument("power proof needs e >= 2");
    std::vector<bool> ops;
    for (long i = static_cast<long>(bitlen(e)) - 2; i >= 0; i--) {
        ops.push_back(true);
        if (mpz_tstbit(e.get_mpz_t(), static_cast<mp_bitcnt_t>(i))) ops.push_back(false);
    }
    return ops;
}

// Proves Cy = Com(x^e; ry) given Cx = Com(x; rx).
inline PowerProof prove_power(const Party& c, const std::string& ctx, int prover, const mpz_class& e,
                              const mpz_class& Cx, const mpz_class& x, const mpz_class& rx,
                              const mpz_class& Cy, const mpz_class& ry) {
    std::vector<bool> ops = power_schedule(e);
    const size_t n = ops.size();
    PowerProof pf;
    // Pass 1 (cheap, sequential): the chain of intermediate values and their commitments.
    struct Step { mpz_class A, B, C, a, ra, s; };
    std::vector<Step> st(n);
    mpz_class cur = c.F().mod(x), rcur = rx, Ccur = Cx;
    for (size_t j = 0; j < n; j++) {
        const bool last = (j + 1 == n);
        const mpz_class& bval = ops[j] ? cur : x;     // second factor
        const mpz_class& rb = ops[j] ? rcur : rx;
        const mpz_class& B = ops[j] ? Ccur : Cx;
        mpz_class nxt = c.F().mul(cur, bval);
        mpz_class rn = last ? ry : c.F().rand();
        mpz_class Cn = last ? Cy : c.ped.commit(nxt, rn);
        if (!last) pf.inter.push_back(Cn);
        st[j] = {Ccur, B, Cn, cur, rcur, c.F().sub(rn, c.F().mul(cur, rb))};   // Cn = B^cur * h^s
        cur = nxt; rcur = rn; Ccur = Cn;
    }
    // Pass 2 (expensive, independent): one product proof per step, spread over all cores.
    pf.steps.resize(n);
    parallel_for(n, [&](size_t j) {
        const Step& t = st[j];
        pf.steps[j] = prove_product(c, ctx + ".step" + std::to_string(j), prover, t.A, t.B, t.C, t.a, t.ra, t.s);
    });
    return pf;
}

inline bool verify_power(const Party& c, const std::string& ctx, int prover, const mpz_class& e,
                         const mpz_class& Cx, const mpz_class& Cy, const PowerProof& pf) {
    std::vector<bool> ops = power_schedule(e);
    if (pf.steps.size() != ops.size() || pf.inter.size() + 1 != ops.size()) return false;
    for (const mpz_class& C : pf.inter)
        if (!c.ped.in_group(C)) return false;
    mpz_class Ccur = Cx;
    for (size_t j = 0; j < ops.size(); j++) {
        const mpz_class& B = ops[j] ? Ccur : Cx;
        const mpz_class& Cn = (j + 1 == ops.size()) ? Cy : pf.inter[j];
        if (!verify_product(c, ctx + ".step" + std::to_string(j), prover, Ccur, B, Cn, pf.steps[j])) return false;
        Ccur = Cn;
    }
    return true;
}

// ---- Batch verification of a power proof plus one product proof ----
// Same result as verify_power + verify_product, but all ~2 * 380 equations are checked
// at once: each equation is raised to a fresh random 128-bit weight and everything is
// multiplied together, so that a single multi-exponentiation suffices (small-exponent
// batch test, Bellare-Garay-Rabin 1998).
// Soundness: every element is first checked to be a quadratic residue mod Q (cheap
// Legendre test). The quadratic residues form a group whose prime factors are all
// > 2^128 (see params.h), so a false equation survives the batch with probability at
// most 2^-128. An honest proof always passes, so a failure proves the proof is invalid.
// The intermediate commitments need no separate order-p test: each one appears as A in
// a product equation, which forces it into G (T1 is fixed before the challenge e).
inline bool verify_power_and_product_batch(const Party& c, const std::string& ctx_pow, const std::string& ctx_prod,
                                           int prover, const mpz_class& e_pub, const mpz_class& Cx,
                                           const mpz_class& Cy, const PowerProof& pf, const mpz_class& A2,
                                           const mpz_class& B2, const mpz_class& C2, const ProductProof& pr) {
    std::vector<bool> ops = power_schedule(e_pub);
    if (pf.steps.size() != ops.size() || pf.inter.size() + 1 != ops.size()) return false;
    for (const mpz_class& C : pf.inter)
        if (!c.ped.is_qr(C)) return false;
    auto well_formed = [&](const ProductProof& s) {
        return c.ped.is_qr(s.T1) && c.ped.is_qr(s.T2) && canonical(c, s.za) && canonical(c, s.zr) && canonical(c, s.zs);
    };
    if (!well_formed(pr)) return false;
    for (const ProductProof& s : pf.steps)
        if (!well_formed(s)) return false;

    // Collect every equation as (statement A, B, C, proof, context).
    struct Eq { const mpz_class *A, *B, *C; const ProductProof* pf; std::string ctx; };
    std::vector<Eq> eqs;
    const mpz_class* Ccur = &Cx;
    for (size_t j = 0; j < ops.size(); j++) {
        const mpz_class* B = ops[j] ? Ccur : &Cx;
        const mpz_class* Cn = (j + 1 == ops.size()) ? &Cy : &pf.inter[j];
        eqs.push_back({Ccur, B, Cn, &pf.steps[j], ctx_pow + ".step" + std::to_string(j)});
        Ccur = Cn;
    }
    eqs.push_back({&A2, &B2, &C2, &pr, ctx_prod});

    // Accumulate exponents per base. For each equation, with weights w1, w2:
    //   g^(w1 za) h^(w1 zr + w2 zs) B^(w2 za)  =  T1^w1 A^(w1 e) T2^w2 C^(w2 e)
    // g and h are in G, so their exponents are reduced mod p; other exponents stay integers.
    std::map<mpz_class, mpz_class> net;   // base -> (left exponent) - (right exponent)
    mpz_class eg = 0, eh = 0;
    const mpz_class two128 = mpz_class(1) << 128;
    std::vector<mpz_class> ch(eqs.size());   // Fiat-Shamir challenges (hashing, done in parallel)
    parallel_for(eqs.size(), [&](size_t i) {
        const Eq& q = eqs[i];
        ch[i] = product_challenge(q.ctx, prover, *q.A, *q.B, *q.C, q.pf->T1, q.pf->T2);
    });
    for (size_t i = 0; i < eqs.size(); i++) {
        const Eq& q = eqs[i];
        const ProductProof& s = *q.pf;
        const mpz_class& e = ch[i];
        mpz_class w1 = rand_below(two128) + 1, w2 = rand_below(two128) + 1;
        eg += w1 * s.za;
        eh += w1 * s.zr + w2 * s.zs;
        net[*q.B] += w2 * s.za;
        net[s.T1] -= w1;
        net[*q.A] -= w1 * e;
        net[s.T2] -= w2;
        net[*q.C] -= w2 * e;
    }
    std::vector<mpz_class> lb, le, rb, re;   // left and right side of the final check
    for (auto& [base, ex] : net) {
        if (ex > 0) { lb.push_back(base); le.push_back(ex); }
        else if (ex < 0) { rb.push_back(base); re.push_back(-ex); }
    }
    mpz_class left = c.ped.mul(c.ped.commit(c.F().mod(eg), c.F().mod(eh)), c.ped.multi_pow(lb, le));
    return left == c.ped.multi_pow(rb, re);
}

// ---- Encoding proofs as values on the log, and reading them back ----
// Verifiers always read the proof from the signed record, so the record itself is the evidence.
inline std::vector<mpz_class> flatten(const PowerProof& pf) {
    std::vector<mpz_class> v(pf.inter);
    for (const ProductProof& s : pf.steps) v.insert(v.end(), {s.T1, s.T2, s.za, s.zr, s.zs});
    return v;
}
inline ProductProof product_at(const std::vector<mpz_class>& v, size_t i) {
    return {v[i], v[i + 1], v[i + 2], v[i + 3], v[i + 4]};
}
// Reads a power proof with `nsteps` steps followed by one product proof.
// Returns false if the record has the wrong length.
inline bool parse_power_and_product(const std::vector<mpz_class>& v, size_t nsteps, PowerProof& pw, ProductProof& pr) {
    const size_t ninter = nsteps - 1;
    if (v.size() != ninter + 5 * nsteps + 5) return false;
    pw.inter.assign(v.begin(), v.begin() + static_cast<long>(ninter));
    pw.steps.clear();
    for (size_t j = 0; j < nsteps; j++) pw.steps.push_back(product_at(v, ninter + 5 * j));
    pr = product_at(v, ninter + 5 * nsteps);
    return true;
}

}  // namespace parry
