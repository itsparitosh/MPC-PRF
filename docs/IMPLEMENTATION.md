# The MPC-PRF implementation: what is where, and where it comes from

This package is the accountable MPC-PRF of Parry on its own: the building
blocks, preprocessing, and the PRF evaluation in both modes (paper: "Parry in
Detail", the subsections Building Blocks to Evaluating the PRF). The full
Parry code (setup and write phases, DPFs) uses exactly these headers.

## 1. Folder layout

Everything is header-only C++17 under `include/`. One folder per layer; a
folder only uses the folders above it in this list.

| folder | phase | files | paper |
|---|---|---|---|
| `base/` | all | `common.h` (randomness, F_p, `Fault`, threads), `serial.h` (bytes), `crypto.h` (SHA-256 transcripts, Ed25519), `params.h` (p, d, Q, g, h) | Parameters |
| `crypto/` | all | `pedersen.h` (commitments in G), `paillier.h`, `ristretto.h` (the group E), `elgamal.h` (hashed ElGamal + Chaum–Pedersen complaint) | Commitments; Sending a value to one party |
| `net/` | all | `net.h` (TCP link), `party.h` (one server, the signed public log) | Notation and the public log |
| `zk/` | prep + init | `sigma.h` (opening, product, power proofs over Pedersen), `cggmp.h` (π^prm, π^mod, π^fac, π^enc, π^aff, π^dec) | Commitments and proofs; Range proofs; Proofs used |
| `mpc/` | online | `sharing.h` (committed sharing ⟨x⟩, checked opening Rec, Beaver Mult) | Committed sharing; Checked opening; Secure multiplication |
| `init/` | init (once) | `keys.h` (keys, π^prm/π^mod/π^fac, PRF key registration, state file) | Paillier keys; Key registration |
| `prep/` | **preprocessing** | `mta.h` (MtA, triples, masks, ζ in one 3-round batch) | Preprocessing |
| `prf/` | **online** | `mpc_prf.h` (key registration, shared mode, delivery mode) | The PRF and the Key; Evaluating the PRF |

## 2. Preprocessing in one picture

`gen_prep_batch(c, n_triples, n_write_masks, n_setup_masks, n_zeta)` in
`prep/mta.h`. Every item, whatever its kind, goes through the same three
rounds, so a whole batch costs three rounds:

| round | each server posts, per item | the other server checks |
|---|---|---|
| A | commitments to its own values; local proof (π^prod for a triple, π^pow + π^prod for a mask, PoK for ζ); for every MtA run where it is S_A: E and π^enc | commitments in G, local proof, π^enc |
| B | for every run where it is S_B: E′, Com(q_B), π^aff | Com(q_B) in G, π^aff |
| C | for every run where it is S_A: Com(q_A), π^dec | Com(q_A) in G, π^dec |

| item | values of S_b | MtA runs (S_A first) | result |
|---|---|---|---|
| triple | α_b, β_b, ψ_b = α_b β_b | (α_0, β_1), (α_1, β_0) | γ = ψ_0 + ψ_1 + both outputs |
| write mask | m_b, t_b = m_b^d, n_b = t_b^-1 | (m_0, m_1), (n_1, n_0) | r = m_0 m_1, v = n_0 n_1 = r^-d |
| setup mask | m_b, u_b = m_b^(2^λ), u_b^-1 | (u_0, u_1) | r = u_0 u_1, r^d = 1 |
| ζ | ζ_b | none | ζ = ζ_0 + ζ_1 |

In Parry, a setup package (one setup mask + one triple) needs 3 MtA runs and a
write package (one write mask + ζ + four triples) needs 10; every mask and
triple is used once.

There is no sacrifice and no replay any more: every MtA message carries a
proof that the other server checks in the round it arrives, so a cheating
server is blamed from its own signed message and nobody reveals anything.

## 3. Parameters used

| name | value | where |
|---|---|---|
| ℓ (field) | 384 bits, p = 2^128·d + 1 | `params.h` |
| λ | 128 | `params.h` |
| G | order-p subgroup of Z_Q^*, Q 3072-bit, Q − 1 = 2·p·q_0⋯q_7 | `params.h` |
| E | Ristretto255 (libsodium) | `crypto/ristretto.h` |
| Paillier N | 3072 bits, two safe primes | `crypto/paillier.h`, `net/party.h` |
| ℓ_π (slack) | 384 (≥ 256 + κ) | `ZkSizes` in `zk/cggmp.h` |
| challenge | 256 bits | `ZkSizes` |
| κ | 80 | `ZkSizes` |
| ℓ_ν | 2ℓ + ℓ_π + κ = 1232 | `ZkSizes` |
| bound in π^dec | honest χ < 2^1618; proven \|χ\| ≤ 2^2003 < N/2 | `ZkSizes::chi`, `mta_modulus_ok` |

## 4. Where every proof and primitive comes from

"Written here" means the code was written for this project from the
description in the cited paper; no code was copied from another library or
project. "Library" means we call a well-known library and do not
reimplement the primitive.

| part | file | origin |
|---|---|---|
| π^prm (ring-Pedersen parameters) | `zk/cggmp.h` | written here, from CGGMP20 (Canetti, Gennaro, Goldfeder, Makriyannis, Peled, CCS 2020, ePrint 2021/060), Π^prm |
| π^mod (Paillier–Blum modulus) | `zk/cggmp.h` | written here, from CGGMP20, Π^mod |
| π^fac (no small factors) | `zk/cggmp.h` | written here, from CGGMP20, Π^fac |
| π^enc (encryption in range) | `zk/cggmp.h` | written here, from CGGMP20 Π^enc; **changed**: the plaintext is tied to a Pedersen commitment g^a h^ρ instead of a curve point a·G (one extra response z4 for ρ) |
| π^aff (affine operation in range) | `zk/cggmp.h` | written here, from CGGMP20 Π^aff-g; **changed** in the same way for x and for B's output share −y (responses v1, v2) |
| π^dec | `zk/cggmp.h` | **ours**: π^enc applied by S_A to E′ with plaintext χ, a larger bound (`ZkSizes::chi`) and its own proof label |
| fast verification of the above | `zk/cggmp.h` (`RpVerifier`), `crypto/paillier.h` (`pow_mod_n`, `pow_mod_n2`) | written here: the verifier uses its own factorization and λ = log_t s to compute the same values with the CRT; the accepted proofs are exactly the same |
| opening proof (PoK) | `zk/sigma.h` | written here; the standard Schnorr/Okamoto proof for a Pedersen commitment (Okamoto, CRYPTO '92) |
| product proof π^prod | `zk/sigma.h` | written here; standard Sigma protocol for a multiplicative relation between Pedersen commitments (Cramer–Damgård, CRYPTO '98) |
| power proof π^pow | `zk/sigma.h` | written here; a chain of product proofs along square-and-multiply |
| batch verification of π^pow | `zk/sigma.h` | written here; small-exponent batch test (Bellare–Garay–Rabin, EUROCRYPT '98), sound because Q − 1 has only large prime factors (`params.h`) |
| Chaum–Pedersen (complaint) | `crypto/elgamal.h` | written here on top of libsodium group operations (Chaum–Pedersen, CRYPTO '92) |
| hashed ElGamal | `crypto/elgamal.h` | written here (Abdalla–Bellare–Rogaway, DHIES) |
| Fiat–Shamir transcripts | `base/crypto.h` | written here; SHA-256 from OpenSSL |
| Ristretto255 group | `crypto/ristretto.h` | library: libsodium (RFC 9496) |
| Ed25519 signatures, SHA-256 | `base/crypto.h` | library: OpenSSL |
| big integers | everywhere | library: GMP |
| Pippenger multi-exponentiation, fixed-base tables, Shamir's trick | `crypto/pedersen.h` | written here, textbook algorithms |

## 5. Optimizations and why they are safe

| optimization | where | why it does not change security |
|---|---|---|
| all items of a batch in the same 3 rounds | `prep/mta.h` | the messages and proofs are exactly those of separate runs; each proof has its own Fiat–Shamir context (batch id, item, run, prover), so nothing can be moved between items |
| MtA outputs are directly the mask shares (no extra Mult) | `prep/mta.h` | the MtA inputs are the committed values the power proofs talk about, and π^enc/π^aff prove that the ciphertexts use exactly those values |
| CRT for everything a server computes under its own key | `crypto/paillier.h`, `zk/cggmp.h` | same numbers as the plain computation (unit test "CRT verifier ..."); only the owner of the factors can do it, and it only helps that owner |
| trapdoor λ for ring-Pedersen checks | `RpVerifier` | s^x t^m = t^(λx + m) holds by construction; the verifier checks its own parameters, so it cannot be fooled by its own trapdoor |
| batch group-membership test and batch verification of π^pow | `crypto/pedersen.h`, `zk/sigma.h` | random 128-bit weights; every element is first checked to be a quadratic residue mod Q, whose group has only prime factors > 2^128, so a false statement passes with probability ≤ 2^-128 |
| multi-threading | `base/common.h` | only local computation is parallel; all network rounds and ids are created on the main thread in the same order on both servers |

What is not done (and not claimed): constant-time arithmetic. GMP's
`mpz_powm` is not constant-time, so a server's own secrets could leak through
timing to a local attacker on the same machine. This is common for research
prototypes; the protocol-level security does not depend on it.

## 6. Test hooks (`--cheat <point>`)

| point | what the server does | caught by |
|---|---|---|
| `triple.psi` | ψ_b ≠ α_b β_b | π^prod, round A |
| `mta.E` | E encrypts a + 1 | π^enc, round A |
| `mta.big` | E encrypts a + p·2^(ℓ+ℓ_π) (Alpha-Rays) | π^enc range check, round A |
| `mta.nu` | ν far too large | π^aff range check, round B |
| `mta.Ep` | E′ uses b + 1 | π^aff, round B |
| `mta.dec` | Com(q_A) holds χ + 1 | π^dec, round C |
| `mask.pow` | a wrong mask factor | π^pow / π^prod, round A |
| `zeta.pok`, `key.pok` | invalid proof of knowledge | round A / init |
| `omega`, `delivery.mul` | wrong share in an opening | checked opening |
| `setup.deliver` | wrong share sent to the writer | writer's complaint |
