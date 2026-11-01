// =====================================================================
// mta.h -- preprocessing: Beaver triples, masks and zeta, all in three rounds
// (paper: "Preprocessing": MtA, Beaver triples, Masks, Rounds and packages)
//
// Everything the online phases consume is built from values that each server
// samples and commits to on its own. Wherever a value of S_0 must be multiplied
// with a value of S_1, the product is turned into shares with MtA:
//
//   S_A knows a_A, S_B knows a_B, both committed on the log.
//   round A  S_A: E  = Enc_A(a_A; eta)                         + pi^enc
//   round B  S_B: E' = E^a_B * Enc_A(nu; xi), Com(q_B = -nu)   + pi^aff
//   round C  S_A: chi = Dec_A(E') (signed), Com(q_A = chi mod p) + pi^dec
//   then q_A + q_B = a_A a_B mod p, and both commitments are on the log.
//
// The proofs (zk/cggmp.h) are checked by the other server in the round they
// arrive, so a server that cheats is blamed from its own signed message, in the
// round in which it cheats. Without them MtA leaks: S_A could encrypt a_A + j p
// for a huge j and read a_B out of chi (Alpha-Rays, ePrint 2021/1621), and S_B
// could make chi wrap around N_A and learn a bit of a_A from the result.
//
// What one batch makes ("items"); each server commits to its own values in round A:
//   Triple     values (alpha_b, beta_b, psi_b = alpha_b beta_b), pi^prod for psi_b.
//              gamma = psi_0 + psi_1 + MtA(alpha_0, beta_1) + MtA(alpha_1, beta_0).
//   WriteMask  values (m_b, t_b = m_b^d, n_b = t_b^-1), pi^pow for t_b, pi^prod for n_b t_b = 1.
//              r = MtA(m_0, m_1), v = MtA(n_1, n_0) = (m_0 m_1)^-d = r^-d.
//   SetupMask  values (m_b, u_b = m_b^(2^lambda), u_b^-1), pi^pow for u_b, pi^prod for u_b u_b^-1 = 1.
//              r = MtA(u_0, u_1), so r^d = (m_0 m_1)^(p-1) = 1.
//   Zeta       value zeta_b with a proof of knowledge of its opening.
// The output shares of MtA are directly the shares of the result (no extra Mult),
// and since pi^enc / pi^aff refer to the same commitments as pi^pow / pi^prod,
// MtA multiplies exactly the values that were proven.
// All items of a batch share the same three rounds, whatever their number.
// =====================================================================
#pragma once
#include "mpc/sharing.h"
#include "zk/sigma.h"
#include <atomic>
#include <chrono>

namespace parry {

struct WriteMask { CShared r, v; };   // v = r^-d

// Is N_A large enough for MtA? Two numbers must stay below N_A / 2:
//  * the bound that pi^dec proves on chi (so the proven chi is the true decryption), and
//  * the true integer a_A a_B + nu for any values that pass pi^enc and pi^aff
//    (one extra bit on each proven bound, from the extraction).
inline bool mta_modulus_ok(const ZkSizes& Z, const mpz_class& N) {
    const mpz_class dec_bound = pow2k(Z.chi + Z.slack + 1);
    const mpz_class true_bound = pow2k(2 * (Z.ell + Z.slack + 1)) + pow2k(Z.nu + Z.slack + 1);
    return 2 * dec_bound < N && 2 * true_bound < N;
}

enum class ItemKind : uint8_t { Triple, WriteMask, SetupMask, Zeta };

// One MtA run inside an item: a_A is value ia of S_A, a_B is value ib of S_B.
struct RunSpec { int A; size_t ia, ib; };

inline size_t item_values(ItemKind k) { return k == ItemKind::Zeta ? 1 : 3; }

inline std::vector<RunSpec> item_runs(ItemKind k) {
    switch (k) {
        case ItemKind::Triple: return {{0, 0, 1}, {1, 0, 1}};      // alpha_0 beta_1, alpha_1 beta_0
        case ItemKind::WriteMask: return {{0, 0, 0}, {1, 2, 2}};   // m_0 m_1, n_1 n_0
        case ItemKind::SetupMask: return {{0, 1, 1}};              // u_0 u_1
        case ItemKind::Zeta: return {};
    }
    return {};
}

inline const char* item_name(ItemKind k) {
    switch (k) {
        case ItemKind::Triple: return "triple";
        case ItemKind::WriteMask: return "wmask";
        case ItemKind::SetupMask: return "smask";
        case ItemKind::Zeta: return "zeta";
    }
    return "?";
}

// Public exponent of the power proof: d for a write mask, 2^lambda for a setup mask.
inline mpz_class pow_exponent(const Party& c, ItemKind k) {
    return k == ItemKind::WriteMask ? c.P.d : mpz_class(mpz_class(1) << c.P.lambda);
}

// Length of the local proof of an item, as it appears in the round-A record.
inline size_t local_proof_len(const Party& c, ItemKind k) {
    switch (k) {
        case ItemKind::Triple: return 5;
        case ItemKind::Zeta: return 3;
        default: {
            const size_t n = power_schedule(pow_exponent(c, k)).size();
            return (n - 1) + 5 * n + 5;     // intermediate commitments, product steps, final product proof
        }
    }
}

// Number of MtA runs of an item in which server b is S_A.
template <class It>
inline size_t runs_of(const It& it, int b) {
    size_t k = 0;
    for (const auto& run : it.runs) k += run.spec.A == b;
    return k;
}

// Fiat-Shamir context of item i, for proofs made by server `b`.
inline std::string item_ctx(const std::string& id, size_t i, int b) { return id + ".i" + std::to_string(i) + ".S" + std::to_string(b); }

// ---- optional timing (benchmarks only) ----
// Preprocessing runs in steps; each step does ONE kind of operation for all items of
// the batch, spread over all cores. The time of a step is wall-clock time on this
// server, so the steps (plus the network) add up exactly to the time of the batch.
enum class Op : int { MtA, TripleLocal, MaskLocal, Commit, ProofEnc, ProofAff, ProofDec, ProofProd, ProofPow, ProofPok,
                      Network, kCount };

struct PrepStats {
    double secs[static_cast<int>(Op::kCount)] = {};
    uint64_t sent = 0;            // bytes this server sent
    double get(Op o) const { return secs[static_cast<int>(o)]; }
    double total() const {
        double t = 0;
        for (double s : secs) t += s;
        return t;
    }
};

// Output of one batch.
struct PrepBatch {
    std::vector<Triple> triples;
    std::vector<WriteMask> wmasks;
    std::vector<CShared> smasks;     // r with r^d = 1
    std::vector<CShared> zetas;
};

// One server's view of an item while it is being built.
struct Item {
    ItemKind kind;
    std::vector<CShared> val;                      // my values; val[k].C[peer] comes from the peer's record
    std::vector<mpz_class> pf_pow, pf_prod, pf_pok; // my local proofs, as they go into the round-A record
    struct Run {
        RunSpec spec;
        mpz_class E, Ep;                           // the two ciphertexts (public)
        mpz_class eta, nu, xi, chi;                // my secrets in this run (as S_A: eta, chi; as S_B: nu, xi)
        mpz_class a_used, b_used;                  // the integers I really used (differ only in test hooks)
        std::vector<mpz_class> pf;                 // my proof of this round (enc, aff or dec)
        CShared out;                               // my share of a_A a_B, with both commitments
    };
    std::vector<Run> runs;
};

// Makes nt triples, nw write masks, ns setup masks and nz values zeta in three rounds.
inline PrepBatch gen_prep_batch(Party& c, size_t nt, size_t nw, size_t ns, size_t nz, PrepStats* stats = nullptr) {
    const int me = c.me, pe = c.peer();
    const ZkSizes& Z = c.Z;
    const Field& F = c.F();
    const PaillierPub& pkme = c.pal[me];
    const PaillierPub& pkpe = c.pal[pe];
    if (!mta_modulus_ok(Z, pkme.N) || !mta_modulus_ok(Z, pkpe.N)) throw std::runtime_error("Paillier modulus too small for MtA");
    const RpVerifier rv = c.rp_verifier();
    const uint64_t sent_before = c.ch ? c.ch->stats.bytes_sent : 0;

    std::vector<Item> I;
    for (size_t i = 0; i < nt; i++) I.push_back({ItemKind::Triple, {}, {}, {}, {}, {}});
    for (size_t i = 0; i < nw; i++) I.push_back({ItemKind::WriteMask, {}, {}, {}, {}, {}});
    for (size_t i = 0; i < ns; i++) I.push_back({ItemKind::SetupMask, {}, {}, {}, {}, {}});
    for (size_t i = 0; i < nz; i++) I.push_back({ItemKind::Zeta, {}, {}, {}, {}, {}});
    const size_t n = I.size();
    PrepBatch out;
    if (n == 0) return out;
    for (Item& it : I)
        for (const RunSpec& s : item_runs(it.kind)) it.runs.push_back({s, 0, 0, 0, 0, 0, 0, 0, 0, {}, CShared{}});
    const std::string id = c.fresh_id("prep");
    auto run_ctx = [&](size_t i, size_t r) { return id + ".i" + std::to_string(i) + ".run" + std::to_string(r); };

    // Index lists, so that each step loops over exactly the items or runs it needs.
    std::vector<size_t> triples, masks, zetas;
    for (size_t i = 0; i < n; i++)
        (I[i].kind == ItemKind::Triple ? triples : I[i].kind == ItemKind::Zeta ? zetas : masks).push_back(i);
    struct RunRef { size_t i, r; };
    std::vector<RunRef> runs_A_me, runs_A_peer;   // runs where I am S_A / where the peer is S_A
    for (size_t i = 0; i < n; i++)
        for (size_t r = 0; r < I[i].runs.size(); r++) (I[i].runs[r].spec.A == me ? runs_A_me : runs_A_peer).push_back({i, r});
    // Every item has at most one run per S_A, so a step over runs touches each item once
    // (the checks below write err[item] from parallel threads).
    for (const Item& it : I)
        if (runs_of(it, 0) > 1 || runs_of(it, 1) > 1) throw std::logic_error("gen_prep_batch: two MtA runs with the same S_A in one item");

    // step(op, list, body): runs body(x) for every x in list on all cores, and adds the
    // wall-clock time to op.
    auto step = [&](Op op, const auto& list, auto body) {
        const auto t0 = std::chrono::steady_clock::now();
        parallel_for(list.size(), [&](size_t k) { body(list[k]); });
        if (stats) stats->secs[static_cast<int>(op)] += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    // Same for work done once on the main thread (group checks, the network).
    auto timed = [&](Op op, auto body) {
        const auto t0 = std::chrono::steady_clock::now();
        auto r = body();
        if (stats) stats->secs[static_cast<int>(op)] += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return r;
    };
    // Throws the first error found by a check step (in item order, so both runs agree).
    auto throw_first = [&](const std::vector<std::string>& err, const Round& rd) {
        for (size_t i = 0; i < n; i++)
            if (!err[i].empty()) throw Fault(pe, err[i], rd.ids[pe][i]);
    };

    // ================= Round A: my values, local proofs, and E + pi^enc as S_A =================
    // Local values. Triple: alpha, beta, psi = alpha beta. Mask: m, m^e, (m^e)^-1 with e = d for a
    // write mask (m, t, n) and e = 2^lambda for a setup mask (m, u, u^-1). Zeta: one random value.
    step(Op::TripleLocal, triples, [&](size_t i) {
        std::vector<mpz_class> x = {F.rand(), F.rand(), 0};
        x[2] = F.mul(x[0], x[1]);
        if (c.cheat("triple.psi")) x[2] = F.add(x[2], 1);                 // test hook
        I[i].val.resize(3);
        for (int k = 0; k < 3; k++) I[i].val[k].x = x[k];
    });
    step(Op::MaskLocal, masks, [&](size_t i) {
        const bool write = I[i].kind == ItemKind::WriteMask;
        std::vector<mpz_class> x(3);
        x[0] = F.rand_nonzero();
        x[1] = F.pow(x[0], pow_exponent(c, I[i].kind));
        if (!write && c.cheat("mask.pow")) x[1] = F.add(x[1], 1);         // test hook (setup mask)
        x[2] = F.inv(x[1]);
        if (write && c.cheat("mask.pow")) x[2] = F.add(x[2], 1);          // test hook (write mask)
        I[i].val.resize(3);
        for (int k = 0; k < 3; k++) I[i].val[k].x = x[k];
    });
    std::vector<size_t> all(n);
    for (size_t i = 0; i < n; i++) all[i] = i;
    step(Op::Commit, all, [&](size_t i) {
        Item& it = I[i];
        if (it.kind == ItemKind::Zeta) { it.val.resize(1); it.val[0].x = F.rand(); }
        for (CShared& v : it.val) {
            v.rho = F.rand();
            v.C[me] = c.ped.commit(v.x, v.rho);
        }
    });
    // pi^prod. Triple: psi = alpha * beta, i.e. C_psi = C_beta^alpha h^s with s = rho_psi - alpha rho_beta.
    // Mask: value 2 times value 1 is 1, i.e. C_1^x2 h^s = g = Com(1; 0) with s = -x2 rho_1; this
    // also shows that the mask factor is not 0.
    std::vector<size_t> prod_items = triples;
    prod_items.insert(prod_items.end(), masks.begin(), masks.end());
    step(Op::ProofProd, prod_items, [&](size_t i) {
        Item& it = I[i];
        const CShared* v = it.val.data();
        const std::string ctx = item_ctx(id, i, me) + ".prod";
        ProductProof p = it.kind == ItemKind::Triple
            ? prove_product(c, ctx, me, v[0].C[me], v[1].C[me], v[2].C[me], v[0].x, v[0].rho, F.sub(v[2].rho, F.mul(v[0].x, v[1].rho)))
            : prove_product(c, ctx, me, v[2].C[me], v[1].C[me], c.P.g, v[2].x, v[2].rho, F.neg(F.mul(v[2].x, v[1].rho)));
        it.pf_prod = {p.T1, p.T2, p.za, p.zr, p.zs};
    });
    // pi^pow: value 1 is value 0 to the power e.
    step(Op::ProofPow, masks, [&](size_t i) {
        Item& it = I[i];
        const CShared* v = it.val.data();
        it.pf_pow = flatten(prove_power(c, item_ctx(id, i, me) + ".pow", me, pow_exponent(c, it.kind), v[0].C[me], v[0].x,
                                        v[0].rho, v[1].C[me], v[1].rho));
    });
    // proof of knowledge of zeta_b
    step(Op::ProofPok, zetas, [&](size_t i) {
        CShared& v = I[i].val[0];
        OpeningProof p = prove_opening(c, item_ctx(id, i, me) + ".pok", me, v.C[me], v.x, v.rho);
        if (c.cheat("zeta.pok")) p.z1 = F.add(p.z1, 1);                    // test hook
        I[i].pf_pok = {p.T, p.z1, p.z2};
    });
    // MtA as S_A: E = Enc_me(a_A; eta).
    // Test hooks: mta.E encrypts a + 1 (does not match the commitment); mta.big encrypts
    // a + p 2^(ell+slack) (matches mod p, but is far too large).
    step(Op::MtA, runs_A_me, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        run.a_used = I[rr.i].val[run.spec.ia].x;
        if (c.cheat("mta.E")) run.a_used = F.add(run.a_used, 1);
        if (c.cheat("mta.big")) run.a_used += c.P.p() * pow2k(Z.ell + Z.slack);
        run.eta = pkme.rand_r();
        run.E = c.pal_key.enc(run.a_used, run.eta);
    });
    step(Op::ProofEnc, runs_A_me, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        const CShared& a = I[rr.i].val[run.spec.ia];
        run.pf = flatten(prove_enc(Z, c.ped, F, c.pal_key, c.rp[pe], run.E, a.C[me], run.a_used, run.eta, a.rho, Z.ell,
                                   run_ctx(rr.i, rr.r) + ".enc", me));
    });
    // Record per item: my commitments, local proofs (pow before prod), then E + pi^enc per run where I am S_A.
    std::vector<std::vector<mpz_class>> recA(n);
    for (size_t i = 0; i < n; i++) {
        const Item& it = I[i];
        std::vector<mpz_class>& rec = recA[i];
        for (const CShared& v : it.val) rec.push_back(v.C[me]);
        for (const std::vector<mpz_class>* p : {&it.pf_pow, &it.pf_prod, &it.pf_pok}) rec.insert(rec.end(), p->begin(), p->end());
        for (const Item::Run& run : it.runs)
            if (run.spec.A == me) {
                rec.push_back(run.E);
                rec.insert(rec.end(), run.pf.begin(), run.pf.end());
            }
    }
    Round ra = timed(Op::Network, [&] { return c.post(id + ".A", recA, n); });

    // ---- check round A of the peer ----
    std::vector<mpz_class> allC;
    for (size_t i = 0; i < n; i++) {
        Item& it = I[i];
        const std::vector<mpz_class>& th = c.peer_values(ra, i);
        const size_t want = it.val.size() + local_proof_len(c, it.kind) + runs_of(it, pe) * (1 + kEncProofLen);
        if (th.size() != want) throw Fault(pe, std::string("prep round A: malformed ") + item_name(it.kind) + " record", ra.ids[pe][i]);
        for (size_t k = 0; k < it.val.size(); k++) {
            it.val[k].C[pe] = th[k];
            allC.push_back(th[k]);
        }
    }
    timed(Op::Commit, [&] { check_peer_commitments(c, allC, ra.first(pe), id + ".A"); return 0; });
    std::vector<std::string> err(n);
    auto peer_local = [&](size_t i) {   // the peer's local proof of item i
        const std::vector<mpz_class>& th = c.peer_values(ra, i);
        const size_t nv = I[i].val.size();
        return std::vector<mpz_class>(th.begin() + static_cast<long>(nv), th.begin() + static_cast<long>(nv + local_proof_len(c, I[i].kind)));
    };
    step(Op::ProofProd, triples, [&](size_t i) {
        const CShared* v = I[i].val.data();
        if (!verify_product(c, item_ctx(id, i, pe) + ".prod", pe, v[0].C[pe], v[1].C[pe], v[2].C[pe], product_at(peer_local(i), 0)))
            err[i] = "invalid proof pi^prod for psi";
    });
    // A mask's power proof and product proof are checked together in one batch (counted as pi^pow).
    step(Op::ProofPow, masks, [&](size_t i) {
        const Item& it = I[i];
        const CShared* v = it.val.data();
        PowerProof pw;
        ProductProof pr;
        const std::string ctx = item_ctx(id, i, pe);
        if (!parse_power_and_product(peer_local(i), power_schedule(pow_exponent(c, it.kind)).size(), pw, pr) ||
            !verify_power_and_product_batch(c, ctx + ".pow", ctx + ".prod", pe, pow_exponent(c, it.kind), v[0].C[pe],
                                            v[1].C[pe], pw, v[2].C[pe], v[1].C[pe], c.P.g, pr))
            err[i] = std::string("invalid pi^pow / pi^prod for a ") + item_name(it.kind);
    });
    step(Op::ProofPok, zetas, [&](size_t i) {
        const std::vector<mpz_class> pf = peer_local(i);
        if (!verify_opening(c, item_ctx(id, i, pe) + ".pok", pe, I[i].val[0].C[pe], {pf[0], pf[1], pf[2]}))
            err[i] = "invalid proof of knowledge for zeta";
    });
    throw_first(err, ra);
    // pi^enc of the peer: its E encrypts the value in its commitment, and that value is small
    auto enc_at = [&](size_t i, size_t r) {   // position of E for run r in the peer's round-A record
        size_t at = I[i].val.size() + local_proof_len(c, I[i].kind);
        for (size_t q = 0; q < r; q++) at += (I[i].runs[q].spec.A == pe) * (1 + kEncProofLen);
        return at;
    };
    step(Op::ProofEnc, runs_A_peer, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        const std::vector<mpz_class>& th = c.peer_values(ra, rr.i);
        const size_t at = enc_at(rr.i, rr.r);
        run.E = th[at];
        EncProof ep;
        if (!parse_enc(th, at + 1, ep) ||
            !verify_enc(Z, c.ped, F, pkpe, rv, run.E, I[rr.i].val[run.spec.ia].C[pe], ep, Z.ell, run_ctx(rr.i, rr.r) + ".enc", pe))
            err[rr.i] = "MtA: invalid proof pi^enc for E";
    });
    throw_first(err, ra);

    // ================= Round B: as S_B, E' + Com(q_B) + pi^aff =================
    step(Op::MtA, runs_A_peer, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        const CShared& b = I[rr.i].val[run.spec.ib];
        run.nu = rand_below(pow2k(Z.nu));
        if (c.cheat("mta.nu")) run.nu = pow2k(Z.nu + Z.slack + 1);                 // test hook: far too large
        run.b_used = c.cheat("mta.Ep") ? F.add(b.x, 1) : b.x;                        // test hook
        run.xi = pkpe.rand_r();
        run.Ep = pkpe.add(pkpe.scale(run.E, run.b_used), pkpe.enc(run.nu, run.xi));
    });
    step(Op::Commit, runs_A_peer, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        run.out.x = F.neg(run.nu);                                                   // my share as S_B: q_B = -nu
        run.out.rho = F.rand();
        run.out.C[me] = c.ped.commit(run.out.x, run.out.rho);
    });
    step(Op::ProofAff, runs_A_peer, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        const CShared& b = I[rr.i].val[run.spec.ib];
        run.pf = flatten(prove_aff(Z, c.ped, F, pkpe, c.rp[pe], run.E, run.Ep, b.C[me], run.out.C[me], run.b_used, run.nu,
                                   run.xi, b.rho, run.out.rho, run_ctx(rr.i, rr.r) + ".aff", me));
    });
    std::vector<std::vector<mpz_class>> recB(n);
    for (size_t i = 0; i < n; i++)
        for (const Item::Run& run : I[i].runs)
            if (run.spec.A == pe) {
                recB[i].push_back(run.Ep);
                recB[i].push_back(run.out.C[me]);
                recB[i].insert(recB[i].end(), run.pf.begin(), run.pf.end());
            }
    Round rb = timed(Op::Network, [&] { return c.post(id + ".B", recB, n); });

    // ---- check round B of the peer ----
    allC.clear();
    for (size_t i = 0; i < n; i++) {
        const std::vector<mpz_class>& th = c.peer_values(rb, i);
        if (th.size() != runs_of(I[i], me) * (2 + kAffProofLen))
            throw Fault(pe, std::string("prep round B: malformed ") + item_name(I[i].kind) + " record", rb.ids[pe][i]);
        size_t at = 0;
        for (Item::Run& run : I[i].runs)
            if (run.spec.A == me) {
                run.Ep = th[at];
                run.out.C[pe] = th[at + 1];                                          // Com(q_B) of the peer
                allC.push_back(th[at + 1]);
                at += 2 + kAffProofLen;
            }
    }
    timed(Op::Commit, [&] { check_peer_commitments(c, allC, rb.first(pe), id + ".B"); return 0; });
    auto aff_at = [&](size_t i, size_t r) {
        size_t at = 0;
        for (size_t q = 0; q < r; q++) at += (I[i].runs[q].spec.A == me) * (2 + kAffProofLen);
        return at;
    };
    step(Op::ProofAff, runs_A_me, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        AffProof ap;
        // E' = E^b (1+N)^nu xi^N for the b the peer committed to, Com(q_B) holds -nu, both small
        if (!parse_aff(c.peer_values(rb, rr.i), aff_at(rr.i, rr.r) + 2, ap) ||
            !verify_aff(Z, c.ped, F, c.pal_key, rv, run.E, run.Ep, I[rr.i].val[run.spec.ib].C[pe], run.out.C[pe], ap,
                        run_ctx(rr.i, rr.r) + ".aff", pe))
            err[rr.i] = "MtA: invalid proof pi^aff for E'";
    });
    throw_first(err, rb);

    // ================= Round C: as S_A, decrypt, Com(q_A) + pi^dec =================
    step(Op::MtA, runs_A_me, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        run.chi = c.pal_key.dec_signed(run.Ep);
    });
    step(Op::Commit, runs_A_me, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        run.out.x = F.mod(c.cheat("mta.dec") ? run.chi + 1 : run.chi);              // my share as S_A: q_A = chi mod p
        run.out.rho = F.rand();
        run.out.C[me] = c.ped.commit(run.out.x, run.out.rho);
    });
    step(Op::ProofDec, runs_A_me, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        if (c.cheat("mta.dec")) {
            // test hook: try to prove the wrong share chi + 1 (the proof cannot verify)
            mpz_class m = run.chi % pkme.N;
            if (m < 0) m += pkme.N;
            run.pf = flatten(prove_enc(Z, c.ped, F, c.pal_key, c.rp[pe], run.Ep, run.out.C[me], run.chi + 1,
                                       c.pal_key.randomness(run.Ep, m), run.out.rho, Z.chi, run_ctx(rr.i, rr.r) + ".dec", me));
        } else {
            run.pf = flatten(prove_dec(Z, c.ped, F, c.pal_key, c.rp[pe], run.Ep, run.out.C[me], run.chi, run.out.rho,
                                       run_ctx(rr.i, rr.r), me));
        }
    });
    std::vector<std::vector<mpz_class>> recC(n);
    for (size_t i = 0; i < n; i++)
        for (const Item::Run& run : I[i].runs)
            if (run.spec.A == me) {
                recC[i].push_back(run.out.C[me]);
                recC[i].insert(recC[i].end(), run.pf.begin(), run.pf.end());
            }
    Round rc = timed(Op::Network, [&] { return c.post(id + ".C", recC, n); });

    // ---- check round C of the peer ----
    allC.clear();
    for (size_t i = 0; i < n; i++) {
        const std::vector<mpz_class>& th = c.peer_values(rc, i);
        if (th.size() != runs_of(I[i], pe) * (1 + kEncProofLen))
            throw Fault(pe, std::string("prep round C: malformed ") + item_name(I[i].kind) + " record", rc.ids[pe][i]);
        size_t at = 0;
        for (Item::Run& run : I[i].runs)
            if (run.spec.A == pe) {
                run.out.C[pe] = th[at];                                              // Com(q_A) of the peer
                allC.push_back(th[at]);
                at += 1 + kEncProofLen;
            }
    }
    timed(Op::Commit, [&] { check_peer_commitments(c, allC, rc.first(pe), id + ".C"); return 0; });
    auto dec_at = [&](size_t i, size_t r) {
        size_t at = 0;
        for (size_t q = 0; q < r; q++) at += (I[i].runs[q].spec.A == pe) * (1 + kEncProofLen);
        return at;
    };
    step(Op::ProofDec, runs_A_peer, [&](RunRef rr) {
        Item::Run& run = I[rr.i].runs[rr.r];
        EncProof dp;
        // Com(q_A) holds chi mod p for the chi that E' really encrypts
        if (!parse_enc(c.peer_values(rc, rr.i), dec_at(rr.i, rr.r) + 1, dp) ||
            !verify_dec(Z, c.ped, F, pkpe, rv, run.Ep, run.out.C[pe], dp, run_ctx(rr.i, rr.r), pe))
            err[rr.i] = "MtA: invalid proof pi^dec for S_A's share";
    });
    throw_first(err, rc);
    if (stats && c.ch) stats->sent += c.ch->stats.bytes_sent - sent_before;

    // ================= assemble =================
    for (Item& it : I) {
        auto run_out = [&](size_t r) -> const CShared& { return it.runs[r].out; };
        switch (it.kind) {
            case ItemKind::Triple: {
                Triple t;
                t.a = it.val[0];
                t.b = it.val[1];
                t.c = add(c, add(c, it.val[2], run_out(0)), run_out(1));   // gamma = psi + q(A=0) + q(A=1)
                out.triples.push_back(t);
                break;
            }
            case ItemKind::WriteMask: out.wmasks.push_back({run_out(0), run_out(1)}); break;
            case ItemKind::SetupMask: out.smasks.push_back(run_out(0)); break;
            case ItemKind::Zeta: out.zetas.push_back(it.val[0]); break;
        }
    }
    return out;
}

// n random values, each with a proof of knowledge of its opening, in one round.
// Used once at init for the PRF key shares K_b. The proof stops a server from
// choosing its commitment as a function of the other server's commitment.
inline std::vector<CShared> random_with_pok(Party& c, size_t n, const std::string& what) {
    std::vector<CShared> out(n);
    if (n == 0) return out;
    const std::string id = c.fresh_id(what);
    std::vector<mpz_class> mine(4 * n);
    parallel_for(n, [&](size_t i) {
        CShared& z = out[i];
        z.x = c.F().rand();
        z.rho = c.F().rand();
        z.C[c.me] = c.ped.commit(z.x, z.rho);
        OpeningProof pf = prove_opening(c, id + "." + std::to_string(i), c.me, z.C[c.me], z.x, z.rho);
        if (c.cheat(what + ".pok")) pf.z1 = c.F().add(pf.z1, 1);   // test hook
        mine[4 * i] = z.C[c.me]; mine[4 * i + 1] = pf.T; mine[4 * i + 2] = pf.z1; mine[4 * i + 3] = pf.z2;
    });
    Round r = c.post1(id + ".pok", mine);
    const long rec = r.first(c.peer());
    const std::vector<mpz_class>& th = c.peer_values(r);
    if (th.size() != 4 * n) throw Fault(c.peer(), what + ": malformed record", rec);
    std::vector<mpz_class> Cs(n);
    for (size_t i = 0; i < n; i++) Cs[i] = th[4 * i];
    check_peer_commitments(c, Cs, rec, id);
    std::vector<char> bad(n, 0);
    parallel_for(n, [&](size_t i) {
        bad[i] = !verify_opening(c, id + "." + std::to_string(i), c.peer(), th[4 * i], {th[4 * i + 1], th[4 * i + 2], th[4 * i + 3]});
    });
    for (size_t i = 0; i < n; i++)
        if (bad[i]) throw Fault(c.peer(), what + ": invalid proof of knowledge", rec);
    for (size_t i = 0; i < n; i++) out[i].C[c.peer()] = Cs[i];
    return out;
}

}  // namespace parry
