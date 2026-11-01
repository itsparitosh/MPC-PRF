// =====================================================================
// mpc_prf.h -- key registration and the accountable MPC-PRF (online)
// (paper: "The PRF and the Key", "Evaluating the PRF", Protocol MPC-PRF)
//
//   F_k(x) = (k + x)^d mod p,   k = K_0 + K_1,   S_b holds K_b
//
// One Beaver multiplication per evaluation: the servers compute
// omega = r (k + x) with a fresh mask r, and then remove r^d:
//   * shared mode (write phase, write mask (r, v = r^-d)): omega is opened,
//     everyone computes theta = omega^d, and S_b sets y_b = theta * v_b, so
//     y = (k + x)^d and nobody learns it. C_{y,b} = C_{v,b}^theta.
//   * delivery mode (setup phase, setup mask with r^d = 1): omega is NOT opened.
//     Each S_b sends (i*, omega_b, rho_{omega,b}) to the writer with hashed
//     ElGamal; the writer checks both shares against the commitments and gets
//     omega^d = (k + x)^d. A bad share leads to a public complaint.
//
// Three rules (paper): a mask is used once; a setup mask is never used in shared
// mode (its omega^d is the PRF output itself); omega = 0 aborts.
// The packages in prep/packages.h are single use, which enforces the first two.
//
// The setup and write phases put the PRF openings into rounds they share with
// other openings, so they call prf_begin / prf_shared_end / delivery_* directly.
// prf_shared_batch runs complete evaluations on its own (tests, benchmarks).
// =====================================================================
#pragma once
#include "crypto/elgamal.h"
#include "prep/mta.h"

namespace parry {

// ---- key registration (once, at init) ----
// S_b picks K_b, posts Com(K_b) with a proof of knowledge; k = K_0 + K_1 is uniform
// as long as one server is honest, and neither server knows it.
inline CShared register_key(Party& c) { return random_with_pok(c, 1, "key")[0]; }

// ---- both modes: s = k + x (local), then the first half of Mult(s, r) ----
inline MulJob prf_begin(const Party& c, const CShared& key, const CShared& x, const CShared& r, const Triple& t) {
    return mul_begin(c, add(c, key, x), r, t);   // C_{s,b} = C_{k,b} C_{x,b}
}

// ---- shared mode: omega is opened (and nonzero); y = omega^d * v ----
inline CShared prf_shared_end(const Party& c, const CShared& v, const mpz_class& omega) {
    const mpz_class theta = c.F().pow(omega, c.P.d);   // public
    return scale(c, v, theta);                          // C_{y,b} = C_{v,b}^theta
}

// Shared mode for a batch of inputs: <y_i> with y_i = F_k(x_i). Two rounds.
inline std::vector<CShared> prf_shared_batch(Party& c, const CShared& key, const std::vector<CShared>& xs,
                                             const std::vector<WriteMask>& masks, const std::vector<Triple>& T) {
    const size_t n = xs.size();
    if (masks.size() < n || T.size() < n) throw std::invalid_argument("prf_shared_batch: one mask and one triple per input");
    std::vector<MulJob> jobs(n);
    for (size_t i = 0; i < n; i++) jobs[i] = prf_begin(c, key, xs[i], masks[i].r, T[i]);
    const std::string id = c.fresh_id("prf");
    Opener op;
    for (const MulJob& j : jobs) op.add(j);
    std::vector<mpz_class> ev = op.run(c, id + ".mul");
    std::vector<CShared> omega(n);
    for (size_t i = 0; i < n; i++) omega[i] = mul_end(c, jobs[i], ev[2 * i], ev[2 * i + 1]);
    std::vector<mpz_class> w = open_many(c, omega, id + ".omega");
    std::vector<CShared> y(n);
    for (size_t i = 0; i < n; i++) {
        if (w[i] == 0) throw std::runtime_error("MPC-PRF: omega = 0, i.e. k + x = 0 (negligible probability); abort");
        y[i] = prf_shared_end(c, masks[i].v, w[i]);
    }
    return y;
}

// ---- delivery mode ----
// The message S_b sends to the writer: (i*, omega_b, rho_{omega,b}), where rho_{omega,b}
// opens the derived commitment C_{omega,b} (S_b knows it from its triple and mask shares).
inline const std::string kDeliveryLabel = "parry/setup-delivery";

inline HCipher delivery_encrypt(const Party& c, const EcPoint& pk, const mpz_class& istar, const CShared& omega) {
    mpz_class w = omega.x;
    if (c.cheat("setup.deliver")) w = c.F().add(w, 1);   // test hook: wrong share to the writer
    return hel_encrypt(c.F(), kDeliveryLabel, pk, {istar, w, omega.rho});
}

// The writer's (and, in a complaint, the public's) check of one decrypted message.
inline bool delivery_ok(const Pedersen& ped, const mpz_class& istar, const mpz_class& Comega_b, const std::vector<mpz_class>& m) {
    return m.size() == 3 && m[0] == istar && ped.opens(Comega_b, m[1], m[2]);
}

// What the writer does with the two ciphertexts of its ACCESSDB entry.
struct DeliveryResult {
    bool ok = false;
    int bad_server = -1;          // the server whose share failed the check (then a complaint follows)
    mpz_class wp;                 // WP(i*) = omega^d, if ok
};

inline DeliveryResult delivery_finish(const Params& P, const Pedersen& ped, const EcKeyPair& kp, const mpz_class& istar,
                                      const HCipher ct[2], const mpz_class Comega[2]) {
    DeliveryResult res;
    mpz_class omega = 0;
    for (int b = 0; b < 2; b++) {
        EcPoint K;
        if (!hel_key(kp, ct[b], K)) { res.bad_server = b; return res; }
        std::vector<mpz_class> m = hel_decrypt_with(P.F, kDeliveryLabel, kp.pk, ct[b], K);
        if (!delivery_ok(ped, istar, Comega[b], m)) { res.bad_server = b; return res; }
        omega = P.F.add(omega, m[1]);
    }
    if (omega == 0) return res;                  // k + i* = 0: negligible; abort
    res.ok = true;
    res.wp = P.F.pow(omega, P.d);                // omega^d = r^d (k + i*)^d = (k + i*)^d, since r^d = 1
    return res;
}

// Public judgement of a writer's complaint about server b's ciphertext.
inline bool delivery_judge(const Params& P, const Pedersen& ped, const std::string& bind, const EcPoint& pk,
                           const mpz_class& istar, const mpz_class& Comega_b, const HCipher& ct, const HComplaint& cp) {
    return judge_complaint(P.F, kDeliveryLabel, bind, pk, ct, cp,
                           [&](const std::vector<mpz_class>& m) { return delivery_ok(ped, istar, Comega_b, m); });
}

}  // namespace parry
