// =====================================================================
// keys.h -- the init phase (once) and what a server keeps on disk
// (paper: "Paillier keys", "Key registration")
//
// <dir>/s<b>.state   my signing key, my Paillier key and ring-Pedersen
//                    trapdoor, the peer's public keys, and my share of the
//                    PRF key with both commitments. Private to S_b.
// =====================================================================
#pragma once
#include "prf/mpc_prf.h"
#include <chrono>

namespace parry {

struct ServerState {
    int me = 0;
    Bytes sig_seed;              // my Ed25519 private key
    Bytes sig_pub[2];
    mpz_class pal_P, pal_Q;      // my Paillier primes
    mpz_class pal_N[2];          // both Paillier moduli
    mpz_class rp_s[2], rp_t[2];  // both ring-Pedersen parameters (their modulus is pal_N)
    mpz_class rp_lam;            // my ring-Pedersen trapdoor (log_t s), secret
    CShared key;                 // my share of the PRF key
};

inline std::string state_path(const std::string& dir, int b) { return dir + "/s" + std::to_string(b) + ".state"; }

inline void save_state(const std::string& dir, const ServerState& s) {
    Writer w;
    w.str("parry-mpcprf-state-v1");
    w.u64(static_cast<uint64_t>(s.me));
    w.bytes(s.sig_seed); w.bytes(s.sig_pub[0]); w.bytes(s.sig_pub[1]);
    w.mpz(s.pal_P); w.mpz(s.pal_Q); w.mpz(s.pal_N[0]); w.mpz(s.pal_N[1]);
    for (int b = 0; b < 2; b++) { w.mpz(s.rp_s[b]); w.mpz(s.rp_t[b]); }
    w.mpz(s.rp_lam);
    write_cs(w, s.key);
    save_file(state_path(dir, s.me), w.buf);
}

inline ServerState load_state(const std::string& dir, int b) {
    Bytes buf = load_file(state_path(dir, b));
    Reader r(buf);
    if (r.str() != "parry-mpcprf-state-v1") throw std::runtime_error("unknown state file format (run init again)");
    ServerState s;
    s.me = static_cast<int>(r.u64());
    if (s.me != b) throw std::runtime_error("state file belongs to the other server");
    s.sig_seed = r.bytes(); s.sig_pub[0] = r.bytes(); s.sig_pub[1] = r.bytes();
    s.pal_P = r.mpz(); s.pal_Q = r.mpz(); s.pal_N[0] = r.mpz(); s.pal_N[1] = r.mpz();
    for (int k = 0; k < 2; k++) { s.rp_s[k] = r.mpz(); s.rp_t[k] = r.mpz(); }
    s.rp_lam = r.mpz();
    s.key = read_cs(r);
    return s;
}

// Puts the keys of a loaded state into a Party.
inline void install(Party& c, const ServerState& s) {
    c.sig = SigningKey::from_seed(s.sig_seed);
    c.sig_pub[0] = s.sig_pub[0];
    c.sig_pub[1] = s.sig_pub[1];
    if (c.sig.pub() != s.sig_pub[c.me]) throw std::runtime_error("state: signing key does not match");
    c.pal_key = PaillierKey::from_primes(s.pal_P, s.pal_Q);
    for (int b = 0; b < 2; b++) {
        c.pal[b].N = s.pal_N[b];
        c.pal[b].N2 = s.pal_N[b] * s.pal_N[b];
        c.rp[b] = {s.pal_N[b], s.rp_s[b], s.rp_t[b]};
    }
    if (c.pal[c.me].N != c.pal_key.pub.N) throw std::runtime_error("state: Paillier key does not match");
    c.rp_lam = s.rp_lam;
    if (powm(c.rp[c.me].t, c.rp_lam, c.rp[c.me].N) != c.rp[c.me].s) throw std::runtime_error("state: ring-Pedersen trapdoor does not match");
}

// ---- Init phase (once): keys, proofs about the Paillier keys, PRF key registration ----
// The new signing key is exchanged in the session handshake; from then on every
// message is signed. Then, as in CGGMP20, each server
//   round 1: posts its Paillier modulus N_b and ring-Pedersen parameters (s_b, t_b),
//            with pi^prm (the parameters are well formed) and pi^mod (N_b is a
//            Paillier-Blum modulus, so gcd(N_b, phi(N_b)) = 1);
//   round 2: posts pi^fac (N_b has no small factors). This proof uses the OTHER
//            server's ring-Pedersen parameters, which is why it needs a second round.
// The peer checks every proof; a bad one is signed evidence against the sender.
// Optional timing of init (benchmarks only), in seconds.
struct InitStats {
    double keygen = 0, prove_prm_mod = 0, check_prm_mod = 0, prove_fac = 0, check_fac = 0, key_reg = 0, network = 0;
};

inline ServerState run_init(Party& c, InitStats* stats = nullptr) {
    using clock = std::chrono::steady_clock;
    clock::time_point tick = clock::now();
    auto lap = [&](double InitStats::*field) {
        const clock::time_point now = clock::now();
        if (stats) stats->*field += std::chrono::duration<double>(now - tick).count();
        tick = now;
    };
    ServerState s;
    s.me = c.me;
    c.sig = SigningKey();
    c.sig_pub[c.me] = c.sig.pub();
    c.sig_pub[c.peer()].clear();
    c.pal_key = PaillierKey::generate(kPaillierBits);
    c.pal[c.me] = c.pal_key.pub;
    const RingPedersenSecret rps = ring_pedersen_from(c.pal_key);
    c.rp[c.me] = rps.pub;
    c.rp_lam = rps.lam;
    lap(&InitStats::keygen);
    c.start_session("init");
    const int me = c.me, pe = c.peer();
    const std::string id = c.fresh_id("init");
    auto ctx = [&](const std::string& what, int b) { return id + "." + what + ".S" + std::to_string(b); };

    // round 1
    Transcript pt("params");
    pt.put("p", c.P.p()); pt.put("Q", c.P.Q); pt.put("g", c.P.g); pt.put("h", c.P.h);
    std::vector<mpz_class> msg = {c.pal_key.pub.N, rps.pub.s, rps.pub.t, pt.challenge_bits("digest", 256)};
    std::vector<mpz_class> prm, mod;
    std::tie(prm, mod) = run_both([&] { return prove_prm(rps, ctx("prm", me), me); },
                                  [&] { return prove_mod(c.pal_key, ctx("mod", me), me); });
    msg.insert(msg.end(), prm.begin(), prm.end());
    msg.insert(msg.end(), mod.begin(), mod.end());
    lap(&InitStats::prove_prm_mod);
    Round r1 = c.post1("init.keys", msg);
    lap(&InitStats::network);
    const std::vector<mpz_class>& th = c.peer_values(r1);
    const long rec1 = r1.first(pe);
    const size_t nprm = 2 * kPrmRounds, nmod = 1 + 3 * kModRounds;
    if (th.size() != 4 + nprm + nmod) throw Fault(pe, "init: malformed record", rec1);
    if (th[3] != msg[3]) throw std::runtime_error("the two servers use different public parameters (params.txt)");
    RingPedersen prp{th[0], th[1], th[2]};
    std::vector<mpz_class> pprm(th.begin() + 4, th.begin() + 4 + static_cast<long>(nprm));
    std::vector<mpz_class> pmod(th.begin() + 4 + static_cast<long>(nprm), th.end());
    bool ok_prm = false, ok_mod = false;
    std::tie(ok_prm, ok_mod) = run_both([&] { return verify_prm(prp, pprm, ctx("prm", pe), pe, kPaillierBits); },
                                        [&] { return verify_mod(prp.N, pmod, ctx("mod", pe), pe); });
    if (!ok_prm) throw Fault(pe, "init: invalid proof pi^prm (ring-Pedersen parameters)", rec1);
    if (!ok_mod) throw Fault(pe, "init: invalid proof pi^mod (not a Paillier-Blum modulus)", rec1);
    lap(&InitStats::check_prm_mod);
    c.pal[pe].N = prp.N;
    c.pal[pe].N2 = prp.N * prp.N;
    c.rp[pe] = prp;

    // round 2: no small factors, proven with the peer's ring-Pedersen parameters
    std::vector<mpz_class> fac = prove_fac(c.Z, c.pal_key, c.rp[pe], ctx("fac", me), me);
    if (c.cheat("init.fac")) fac[6] = zz(unzz(fac[6]) + 1);   // test hook
    lap(&InitStats::prove_fac);
    Round r2 = c.post1("init.fac", fac);
    lap(&InitStats::network);
    if (!verify_fac(c.Z, c.pal[pe].N, c.rp[me], c.peer_values(r2), ctx("fac", pe), pe))
        throw Fault(pe, "init: invalid proof pi^fac (modulus may have small factors)", r2.first(pe));

    lap(&InitStats::check_fac);
    s.key = register_key(c);
    lap(&InitStats::key_reg);
    s.sig_seed = c.sig.seed();
    s.sig_pub[0] = c.sig_pub[0];
    s.sig_pub[1] = c.sig_pub[1];
    s.rp_lam = c.rp_lam;
    s.pal_P = c.pal_key.prime_P();
    s.pal_Q = c.pal_key.prime_Q();
    for (int b = 0; b < 2; b++) {
        s.pal_N[b] = c.pal[b].N;
        s.rp_s[b] = c.rp[b].s;
        s.rp_t[b] = c.rp[b].t;
    }
    return s;
}

}  // namespace parry
