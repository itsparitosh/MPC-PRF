// =====================================================================
// elgamal.h -- sending a value to one party, with a public complaint
// (paper: Building Blocks, "Sending a value to one party")
//
// Hashed ElGamal (Abdalla-Bellare-Rogaway, DHIES) in E = Ristretto255:
//   receiver key   sk random, pk = g_hat^sk
//   encryption     R = g_hat^u, K = pk^u (one-time key),
//                  c_j = m_j + H(label, pk, R, K, j) mod p      for each element m_j of F_p
//   decryption     K = R^sk, then subtract the same pads.
// The ciphertext is posted on the log, so it is fixed and public.
//
// Complaint: if the decrypted value fails a public check, the receiver posts K
// and a Chaum-Pedersen proof that log_g_hat(pk) = log_R(K). Anyone can then
// decrypt that one ciphertext and redo the check. A failed check is evidence
// against the sender; a correct one rejects the complaint. The receiver cannot
// frame the sender (the ciphertext and K fix the plaintext), and the complaint
// reveals only this one message, never sk.
//
// Our own code on top of libsodium's group operations. Chaum-Pedersen
// (CRYPTO '92) is the standard Sigma protocol for equality of discrete logs,
// made non-interactive with Fiat-Shamir.
// =====================================================================
#pragma once
#include "base/crypto.h"
#include "crypto/ristretto.h"

namespace parry {

struct EcKeyPair {
    EcScalar sk;
    EcPoint pk;
    static EcKeyPair generate() {
        EcKeyPair k;
        k.sk = ec_random_scalar();
        k.pk = ec_base(k.sk);
        return k;
    }
};

struct HCipher {
    EcPoint R;
    std::vector<mpz_class> c;   // one element of F_p per message element
};

// Pad j for a message under `label`. The pad depends on pk and R as well as K, so a
// ciphertext cannot be moved to another receiver or combined with another R.
inline mpz_class hel_pad(const Field& F, const std::string& label, const EcPoint& pk, const EcPoint& R,
                         const EcPoint& K, size_t j) {
    Transcript t("hashed-elgamal");
    t.put_str(label);
    t.put_bytes(pk.data(), pk.size());
    t.put_bytes(R.data(), R.size());
    t.put_bytes(K.data(), K.size());
    t.put("j", static_cast<int64_t>(j));
    return t.challenge_mod("pad", F.p);
}

inline HCipher hel_encrypt(const Field& F, const std::string& label, const EcPoint& pk, const std::vector<mpz_class>& m) {
    if (!ec_valid(pk)) throw std::invalid_argument("hel_encrypt: invalid receiver key");
    const EcScalar u = ec_random_scalar();
    HCipher ct;
    ct.R = ec_base(u);
    EcPoint K;
    if (!ec_mul(K, pk, u)) throw std::logic_error("hel_encrypt: identity key");   // impossible for a valid pk
    for (size_t j = 0; j < m.size(); j++) ct.c.push_back(F.add(m[j], hel_pad(F, label, pk, ct.R, K, j)));
    return ct;
}

// Decryption with the one-time key K (the receiver gets K = R^sk; in a complaint K is public).
inline std::vector<mpz_class> hel_decrypt_with(const Field& F, const std::string& label, const EcPoint& pk,
                                               const HCipher& ct, const EcPoint& K) {
    std::vector<mpz_class> m;
    for (size_t j = 0; j < ct.c.size(); j++) m.push_back(F.sub(ct.c[j], hel_pad(F, label, pk, ct.R, K, j)));
    return m;
}

// K = R^sk. False if R is not a valid point (then the ciphertext is malformed).
inline bool hel_key(const EcKeyPair& kp, const HCipher& ct, EcPoint& K) {
    return ec_valid(ct.R) && ec_mul(K, ct.R, kp.sk);
}

// ---- Chaum-Pedersen: log_g_hat(pk) = log_R(K) ----
struct DleqProof {
    EcPoint T1, T2;
    EcScalar z;
};

// `bind` names the ciphertext the complaint is about, so the proof cannot be reused for another one.
inline EcScalar dleq_challenge(const std::string& bind, const EcPoint& pk, const EcPoint& R, const EcPoint& K,
                               const EcPoint& T1, const EcPoint& T2) {
    Transcript t("dleq-ristretto255");
    t.put_str(bind);
    for (const EcPoint* P : {&pk, &R, &K, &T1, &T2}) t.put_bytes(P->data(), P->size());
    mpz_class e = t.challenge_bits("e", 512);
    std::array<uint8_t, 64> h{};
    mpz_to_bytes(e, h.data(), h.size());
    return ec_scalar_from_hash(h.data());
}

inline DleqProof prove_dleq(const std::string& bind, const EcKeyPair& kp, const EcPoint& R, const EcPoint& K) {
    const EcScalar u = ec_random_scalar();
    DleqProof pf;
    pf.T1 = ec_base(u);
    if (!ec_mul(pf.T2, R, u)) throw std::invalid_argument("prove_dleq: invalid R");
    const EcScalar e = dleq_challenge(bind, kp.pk, R, K, pf.T1, pf.T2);
    pf.z = ec_scalar_add(u, ec_scalar_mul(e, kp.sk));
    return pf;
}

// Checks g_hat^z = T1 + e pk and R^z = T2 + e K (additive notation).
inline bool verify_dleq(const std::string& bind, const EcPoint& pk, const EcPoint& R, const EcPoint& K, const DleqProof& pf) {
    for (const EcPoint* P : {&pk, &R, &K, &pf.T1, &pf.T2})
        if (!ec_valid(*P)) return false;
    const EcScalar e = dleq_challenge(bind, pk, R, K, pf.T1, pf.T2);
    EcPoint l1, l2, ePk, eK, r1, r2;
    try { l1 = ec_base(pf.z); } catch (const std::invalid_argument&) { return false; }
    if (!ec_mul(l2, R, pf.z) || !ec_mul(ePk, pk, e) || !ec_mul(eK, K, e)) return false;
    if (!ec_add(r1, pf.T1, ePk) || !ec_add(r2, pf.T2, eK)) return false;
    return l1 == r1 && l2 == r2;
}

// A complaint about one ciphertext: the one-time key and the proof that it is right.
struct HComplaint {
    EcPoint K;
    DleqProof proof;
};

inline HComplaint make_complaint(const std::string& bind, const EcKeyPair& kp, const HCipher& ct) {
    HComplaint cp;
    if (!hel_key(kp, ct, cp.K)) throw std::invalid_argument("make_complaint: malformed ciphertext");
    cp.proof = prove_dleq(bind, kp, ct.R, cp.K);
    return cp;
}

// Public judgement of a complaint from public data only: the receiver's pk, the posted
// ciphertext (whose signature the caller has checked), the complaint, and the check the
// plaintext must pass. Returns true if the SENDER is at fault.
//   * a ciphertext with an invalid R is the sender's fault (no complaint needed);
//   * a complaint with a wrong K or proof is rejected (the receiver lied);
//   * otherwise the sender is at fault exactly when the decrypted value fails the check.
template <class Check>
inline bool judge_complaint(const Field& F, const std::string& label, const std::string& bind, const EcPoint& pk,
                            const HCipher& ct, const HComplaint& cp, Check plaintext_ok) {
    if (!ec_valid(pk)) return false;          // senders never encrypt to an invalid key
    if (!ec_valid(ct.R)) return true;
    if (!verify_dleq(bind, pk, ct.R, cp.K, cp.proof)) return false;
    return !plaintext_ok(hel_decrypt_with(F, label, pk, ct, cp.K));
}

}  // namespace parry
