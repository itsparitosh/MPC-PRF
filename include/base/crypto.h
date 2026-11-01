// =====================================================================
// crypto.h -- hashing and signatures (both from OpenSSL)
//   * Transcript: SHA-256 based Fiat-Shamir hashing. Used to turn interactive
//     proofs into non-interactive ones, and to derive public random challenges.
//   * SigningKey / ed25519_verify: every server message on the public log is
//     signed, so a bad message is evidence that cannot be denied.
// =====================================================================
#pragma once
#include "base/common.h"
#include <openssl/evp.h>
#include <memory>

namespace parry {

// Fiat-Shamir transcript. Every appended item is length-prefixed, so encodings are unambiguous.
class Transcript {
public:
    explicit Transcript(const std::string& domain) { put_str("parry-v1"); put_str(domain); }

    void put_str(const std::string& s) {
        put_u64(s.size());
        buf_.insert(buf_.end(), s.begin(), s.end());
    }
    void put_bytes(const uint8_t* p, size_t n) {
        put_u64(n);
        buf_.insert(buf_.end(), p, p + n);
    }
    void put_u64(uint64_t v) {
        for (int i = 7; i >= 0; i--) buf_.push_back(static_cast<unsigned char>(v >> (8 * i)));
    }
    void put(const std::string& label, const mpz_class& x) {
        if (x < 0) throw std::invalid_argument("Transcript: negative value");
        put_str(label);
        size_t n = 0;
        std::vector<unsigned char> b((bitlen(x) + 7) / 8 + 1);
        mpz_export(b.data(), &n, 1, 1, 0, 0, x.get_mpz_t());
        put_u64(n);
        buf_.insert(buf_.end(), b.begin(), b.begin() + n);
    }
    void put(const std::string& label, int64_t v) { put_str(label); put_u64(static_cast<uint64_t>(v)); }

    // Challenge uniform in [0, 2^bits). The challenge is absorbed back, so later challenges differ.
    mpz_class challenge_bits(const std::string& label, unsigned bits) {
        mpz_class c = expand(label, bits) % (mpz_class(1) << bits);
        put(label + ".out", c);
        return c;
    }
    // Challenge (statistically) uniform in [0, m).
    mpz_class challenge_mod(const std::string& label, const mpz_class& m) {
        mpz_class c = expand(label, bitlen(m) + 128) % m;
        put(label + ".out", c);
        return c;
    }

private:
    std::vector<unsigned char> buf_;

    mpz_class expand(const std::string& label, unsigned bits) {
        put_str(label);
        std::vector<unsigned char> out;
        for (uint64_t ctr = 0; out.size() * 8 < bits; ctr++) {
            std::vector<unsigned char> in = buf_;
            for (int i = 7; i >= 0; i--) in.push_back(static_cast<unsigned char>(ctr >> (8 * i)));
            unsigned char h[32];
            unsigned int hl = 0;
            if (EVP_Digest(in.data(), in.size(), h, &hl, EVP_sha256(), nullptr) != 1 || hl != 32)
                throw std::runtime_error("SHA-256 failed");
            out.insert(out.end(), h, h + 32);
        }
        mpz_class r;
        mpz_import(r.get_mpz_t(), out.size(), 1, 1, 0, 0, out.data());
        return r;
    }
};

// SHA-256 of a byte string.
inline std::vector<unsigned char> sha256(const unsigned char* data, size_t n) {
    std::vector<unsigned char> h(32);
    unsigned int hl = 0;
    if (EVP_Digest(data, n, h.data(), &hl, EVP_sha256(), nullptr) != 1 || hl != 32)
        throw std::runtime_error("SHA-256 failed");
    return h;
}
inline std::vector<unsigned char> sha256(const std::vector<unsigned char>& b) { return sha256(b.data(), b.size()); }

// Ed25519 signing key (one per server). Public key is 32 bytes.
class SigningKey {
public:
    // Loads a key from its 32-byte private seed (kept in the server's state file).
    static SigningKey from_seed(const std::vector<unsigned char>& seed) {
        SigningKey k(0);
        EVP_PKEY* p = EVP_PKEY_new_raw_private_key(EVP_PKEY_ED25519, nullptr, seed.data(), seed.size());
        if (!p) throw std::runtime_error("Ed25519: bad private key");
        k.key_.reset(p);
        k.export_pub();
        return k;
    }
    std::vector<unsigned char> seed() const {
        std::vector<unsigned char> s(32);
        size_t len = 32;
        if (EVP_PKEY_get_raw_private_key(key_.get(), s.data(), &len) != 1 || len != 32)
            throw std::runtime_error("Ed25519 private key export failed");
        return s;
    }
    SigningKey(SigningKey&&) = default;
    SigningKey& operator=(SigningKey&&) = default;

    SigningKey() {
        std::unique_ptr<EVP_PKEY_CTX, decltype(&EVP_PKEY_CTX_free)> ctx(
            EVP_PKEY_CTX_new_id(EVP_PKEY_ED25519, nullptr), EVP_PKEY_CTX_free);
        EVP_PKEY* k = nullptr;
        if (!ctx || EVP_PKEY_keygen_init(ctx.get()) != 1 || EVP_PKEY_keygen(ctx.get(), &k) != 1)
            throw std::runtime_error("Ed25519 keygen failed");
        key_.reset(k);
        export_pub();
    }
    const std::vector<unsigned char>& pub() const { return pub_; }

    std::vector<unsigned char> sign(const std::vector<unsigned char>& msg) const {
        std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> m(EVP_MD_CTX_new(), EVP_MD_CTX_free);
        size_t len = 0;
        if (!m || EVP_DigestSignInit(m.get(), nullptr, nullptr, nullptr, key_.get()) != 1 ||
            EVP_DigestSign(m.get(), nullptr, &len, msg.data(), msg.size()) != 1)
            throw std::runtime_error("Ed25519 sign failed");
        std::vector<unsigned char> sig(len);
        if (EVP_DigestSign(m.get(), sig.data(), &len, msg.data(), msg.size()) != 1)
            throw std::runtime_error("Ed25519 sign failed");
        sig.resize(len);
        return sig;
    }

private:
    std::unique_ptr<EVP_PKEY, decltype(&EVP_PKEY_free)> key_{nullptr, EVP_PKEY_free};
    std::vector<unsigned char> pub_;

    explicit SigningKey(int) {}
    void export_pub() {
        size_t len = 32;
        pub_.resize(32);
        if (EVP_PKEY_get_raw_public_key(key_.get(), pub_.data(), &len) != 1 || len != 32)
            throw std::runtime_error("Ed25519 public key export failed");
    }
};

inline bool ed25519_verify(const std::vector<unsigned char>& pub, const std::vector<unsigned char>& msg,
                           const std::vector<unsigned char>& sig) {
    std::unique_ptr<EVP_PKEY, decltype(&EVP_PKEY_free)> k(
        EVP_PKEY_new_raw_public_key(EVP_PKEY_ED25519, nullptr, pub.data(), pub.size()), EVP_PKEY_free);
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> m(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    if (!k || !m || EVP_DigestVerifyInit(m.get(), nullptr, nullptr, nullptr, k.get()) != 1) return false;
    return EVP_DigestVerify(m.get(), sig.data(), sig.size(), msg.data(), msg.size()) == 1;
}

}  // namespace parry
