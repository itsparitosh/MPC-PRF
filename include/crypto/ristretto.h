// =====================================================================
// ristretto.h -- the prime-order group E (Ristretto255) for encryption keys
// (paper: Building Blocks, "Parameters": E with generator g_hat)
//
// The writers' keys and the servers' encryption keys live in E, not in the
// commitment group G: same 128-bit security with 32-byte keys instead of 384
// bytes. All group operations are libsodium's (crypto_core_ristretto255_* and
// crypto_scalarmult_ristretto255_*, RFC 9496); this file is only a thin wrapper.
// We write the group additively, as libsodium does: the paper's g_hat^x is
// base(x), and pk^u is mul(pk, u).
//
// Points are checked on the way in: the encoding must be canonical (libsodium
// rejects anything else) and the point must not be the identity.
// =====================================================================
#pragma once
#include "base/serial.h"
#include <sodium.h>
#include <array>

namespace parry {

using EcPoint = std::array<uint8_t, crypto_core_ristretto255_BYTES>;          // 32 bytes
using EcScalar = std::array<uint8_t, crypto_core_ristretto255_SCALARBYTES>;   // 32 bytes, mod the group order

inline void sodium_ready() {
    static const bool ok = sodium_init() >= 0;
    if (!ok) throw std::runtime_error("libsodium failed to initialise");
}

inline bool ec_is_identity(const EcPoint& P) {
    for (uint8_t b : P)
        if (b) return false;
    return true;
}

// A point we accept from someone else: canonical encoding and not the identity.
inline bool ec_valid(const EcPoint& P) {
    sodium_ready();
    return crypto_core_ristretto255_is_valid_point(P.data()) == 1 && !ec_is_identity(P);
}

// Uniform nonzero scalar.
inline EcScalar ec_random_scalar() {
    sodium_ready();
    EcScalar s;
    crypto_core_ristretto255_scalar_random(s.data());   // libsodium never returns 0 here
    return s;
}

// Scalar from 64 bytes of hash output (reduced mod the group order; bias ~2^-260).
inline EcScalar ec_scalar_from_hash(const uint8_t* h64) {
    EcScalar s;
    crypto_core_ristretto255_scalar_reduce(s.data(), h64);
    return s;
}

inline EcScalar ec_scalar_add(const EcScalar& a, const EcScalar& b) {
    EcScalar r;
    crypto_core_ristretto255_scalar_add(r.data(), a.data(), b.data());
    return r;
}
inline EcScalar ec_scalar_mul(const EcScalar& a, const EcScalar& b) {
    EcScalar r;
    crypto_core_ristretto255_scalar_mul(r.data(), a.data(), b.data());
    return r;
}

// g_hat^s. Throws if s = 0 (libsodium refuses to return the identity).
inline EcPoint ec_base(const EcScalar& s) {
    sodium_ready();
    EcPoint P;
    if (crypto_scalarmult_ristretto255_base(P.data(), s.data()) != 0) throw std::invalid_argument("ec_base: zero scalar");
    return P;
}

// P^s. Returns false if the result is the identity (only for s = 0, since E has prime order).
inline bool ec_mul(EcPoint& out, const EcPoint& P, const EcScalar& s) {
    sodium_ready();
    return crypto_scalarmult_ristretto255(out.data(), s.data(), P.data()) == 0;
}

inline bool ec_add(EcPoint& out, const EcPoint& P, const EcPoint& Q) {
    return crypto_core_ristretto255_add(out.data(), P.data(), Q.data()) == 0;
}

// Points and scalars travel in log records as non-negative integers (32 bytes, big-endian).
inline mpz_class ec_to_mpz(const EcPoint& P) { return bytes_to_mpz(P.data(), P.size()); }
inline bool ec_from_mpz(const mpz_class& x, EcPoint& P) { return mpz_to_bytes(x, P.data(), P.size()); }
inline mpz_class ec_scalar_to_mpz(const EcScalar& s) { return bytes_to_mpz(s.data(), s.size()); }
// The byte string is kept as is (libsodium's little-endian encoding), so the round trip is
// exact. A scalar from someone else must be canonical, i.e. below the group order:
// reducing it (zero-extended to 64 bytes) must give it back unchanged.
inline bool ec_scalar_from_mpz(const mpz_class& x, EcScalar& s) {
    if (!mpz_to_bytes(x, s.data(), s.size())) return false;
    std::array<uint8_t, 64> wide{};
    std::copy(s.begin(), s.end(), wide.begin());
    EcScalar red;
    crypto_core_ristretto255_scalar_reduce(red.data(), wide.data());
    return red == s;
}

}  // namespace parry
