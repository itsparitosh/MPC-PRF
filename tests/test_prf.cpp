// =====================================================================
// test_prf.cpp -- unit tests of the MPC-PRF and everything it is built on
//
// Every test runs S0 and S1 as two processes (fork + socketpair) with the same
// Party/Channel code as the real servers. Cheating tests start one side with a
// --cheat point and check that BOTH sides name the right server (or, where the
// check is done only by the honest server, that the honest one does).
// The first run creates test keys in tests/keys (takes a few seconds).
// =====================================================================
#include "harness.h"
#include "prf/mpc_prf.h"
#include <chrono>

using namespace parry;

static int g_pass = 0, g_fail = 0;

static void check(bool ok, const std::string& name, const std::string& detail = "") {
    if (ok) { g_pass++; std::cout << "  PASS  " << name << std::endl; }
    else { g_fail++; std::cout << "  FAIL  " << name << (detail.empty() ? "" : "   [" + detail + "]") << std::endl; }
}
static std::string both(const Outcome& o) { return o.s[0].substr(0, 90) + " | " + o.s[1].substr(0, 90); }
static bool blames(const std::string& s, int b) { return s.rfind("BLAME S" + std::to_string(b), 0) == 0; }

int main() {
    Params P = get_params("params.txt");
    const std::string dir = "tests/keys";
    std::system(("mkdir -p " + dir).c_str());
    if (!file_exists(dir + "/s0.state") || !file_exists(dir + "/s1.state")) {
        Outcome o = two_party(P, dir, nullptr, -1, {}, true);
        if (!o.ok()) { std::cout << "creating test keys failed: " << both(o) << std::endl; return 1; }
    }
    Pedersen ped(P);

    std::cout << "== init: proofs about the Paillier keys ==" << std::endl;
    {
        const std::string d2 = "tests/keys_tmp";
        std::system(("rm -rf " + d2 + " && mkdir -p " + d2).c_str());
        Outcome o = two_party(P, d2, nullptr, 0, {"init.fac"}, true);
        check(blames(o.s[1], 0), "init: a bad pi^fac by S0 is blamed by S1", both(o));
        std::system(("rm -rf " + d2).c_str());
        // the proofs themselves, on the keys made above
        ServerState s0 = load_state(dir, 0);
        PaillierKey k = PaillierKey::from_primes(s0.pal_P, s0.pal_Q);
        RingPedersenSecret rps = ring_pedersen_from(k);
        RingPedersen other{s0.pal_N[1], s0.rp_s[1], s0.rp_t[1]};
        ZkSizes Z(P);
        std::vector<mpz_class> pm = prove_prm(rps, "t", 0), md = prove_mod(k, "t", 0), fc = prove_fac(Z, k, other, "t", 0);
        bool good = verify_prm(rps.pub, pm, "t", 0, kPaillierBits) && verify_mod(k.pub.N, md, "t", 0) &&
                    verify_fac(Z, k.pub.N, other, fc, "t", 0);
        bool ctx_bound = !verify_mod(k.pub.N, md, "u", 0) && !verify_prm(rps.pub, pm, "t", 1, kPaillierBits);
        // a proof made for N must not pass for another modulus
        mpz_class Nbad = k.pub.N + 2;
        bool rejects = !verify_mod(Nbad, md, "t", 0) && !verify_fac(Z, Nbad, other, fc, "t", 0);
        std::vector<mpz_class> fc2 = fc;
        fc2[6] = zz(unzz(fc2[6]) + 1);
        check(good && ctx_bound && rejects && !verify_fac(Z, k.pub.N, other, fc2, "t", 0),
              "pi^prm, pi^mod, pi^fac: accept honest, reject wrong modulus, context or response");
        // The fast verifier (CRT + trapdoor) must give exactly the slow results.
        RpVerifier rv{rps.pub, rps.lam, &k};
        bool fast_ok = true;
        for (int i = 0; i < 20 && fast_ok; i++) {
            mpz_class x = rand_signed(pow2k(2000)), m = rand_signed(pow2k(3500)), base = k.pub.rand_r();
            mpz_class e = rand_signed(pow2k(800)), b2 = k.pub.rand_r() * k.pub.rand_r() % k.pub.N2;
            fast_ok = rv.commit(x, m) == rp_commit(rps.pub, x, m) && rv.pow(base, e) == powm(base, e, k.pub.N) &&
                      k.pow_mod_n2(b2, e) == powm(b2, e, k.pub.N2);
        }
        check(fast_ok, "CRT verifier gives exactly the same values as the plain computation");
        // pi^enc, pi^aff and pi^dec directly. In real runs the verifier is the other server;
        // here one key plays both roles, which is enough to test the algebra and the bounds.
        mpz_class a = P.F.rand(), rho = P.F.rand(), C = ped.commit(a, rho), eta = k.pub.rand_r();
        mpz_class E = k.enc(a, eta);
        EncProof pe = prove_enc(Z, ped, P.F, k, rps.pub, E, C, a, eta, rho, Z.ell, "t", 0);
        bool enc_ok = verify_enc(Z, ped, P.F, k.pub, rv, E, C, pe, Z.ell, "t", 0);
        mpz_class abig = a + P.p() * pow2k(Z.ell + Z.slack);
        mpz_class Ebig = k.enc(abig, eta);
        EncProof pb = prove_enc(Z, ped, P.F, k, rps.pub, Ebig, C, abig, eta, rho, Z.ell, "t", 0);
        bool enc_big = !verify_enc(Z, ped, P.F, k.pub, rv, Ebig, C, pb, Z.ell, "t", 0);
        check(enc_ok && enc_big && !verify_enc(Z, ped, P.F, k.pub, rv, E, ped.commit(a + 1, rho), pe, Z.ell, "t", 0) &&
                  !verify_enc(Z, ped, P.F, k.pub, rv, E, C, pe, Z.ell, "t", 1),
              "pi^enc: accepts honest, rejects too large plaintext, wrong commitment, wrong prover");
        mpz_class b = P.F.rand(), rb = P.F.rand(), Cb = ped.commit(b, rb);
        mpz_class nu = rand_below(pow2k(Z.nu)), rq = P.F.rand(), Cq = ped.commit(P.F.neg(nu), rq), xi = k.pub.rand_r();
        mpz_class Ep = k.pub.add(k.pub.scale(E, b), k.pub.enc(nu, xi));
        AffProof pa = prove_aff(Z, ped, P.F, k.pub, rps.pub, E, Ep, Cb, Cq, b, nu, xi, rb, rq, "t", 1);
        bool aff_ok = verify_aff(Z, ped, P.F, k, rv, E, Ep, Cb, Cq, pa, "t", 1);
        mpz_class chi = k.dec_signed(Ep);
        bool sum_ok = P.F.add(P.F.mod(chi), P.F.neg(nu)) == P.F.mul(a, b);
        bool aff_bad = !verify_aff(Z, ped, P.F, k, rv, E, Ep, Cb, ped.commit(P.F.neg(nu + 1), rq), pa, "t", 1);
        check(aff_ok && sum_ok && aff_bad, "pi^aff: accepts honest, rejects wrong output commitment; q_A + q_B = a b");
        // pi^dec: Com(chi mod p) with the true chi passes; a different integer with the same
        // ciphertext (chi + N) or a wrong commitment does not; a pi^enc is not a pi^dec.
        mpz_class rqa = P.F.rand(), CqA = ped.commit(P.F.mod(chi), rqa);
        EncProof pd = prove_dec(Z, ped, P.F, k, rps.pub, Ep, CqA, chi, rqa, "t", 0);
        bool dec_ok = verify_dec(Z, ped, P.F, k.pub, rv, Ep, CqA, pd, "t", 0);
        mpz_class chi2 = chi + k.pub.N, Cq2 = ped.commit(P.F.mod(chi2), rqa);
        EncProof pd2 = prove_enc(Z, ped, P.F, k, rps.pub, Ep, Cq2, chi2, k.randomness(Ep, chi % k.pub.N), rqa, Z.chi, "t.dec", 0);
        bool dec_wrap = !verify_dec(Z, ped, P.F, k.pub, rv, Ep, Cq2, pd2, "t", 0);
        bool dec_wrongC = !verify_dec(Z, ped, P.F, k.pub, rv, Ep, ped.commit(P.F.mod(chi + 1), rqa), pd, "t", 0);
        EncProof pe2 = prove_enc(Z, ped, P.F, k, rps.pub, Ep, CqA, chi, k.randomness(Ep, chi % k.pub.N), rqa, Z.chi, "t", 0);
        bool not_enc = !verify_dec(Z, ped, P.F, k.pub, rv, Ep, CqA, pe2, "t", 0);
        check(dec_ok && dec_wrap && dec_wrongC && not_enc && mta_modulus_ok(Z, k.pub.N),
              "pi^dec: accepts the true decryption, rejects chi + N, a wrong commitment and a relabelled pi^enc");
        // Every single field of a valid proof, changed by one, must make it fail.
        bool enc_tamper = true, aff_tamper = true;
        std::vector<mpz_class> fe = flatten(pe), fa = flatten(pa);
        for (size_t j = 0; j < fe.size(); j++) {
            std::vector<mpz_class> t = fe;
            t[j] += 1;
            EncProof q;
            enc_tamper = enc_tamper && (!parse_enc(t, 0, q) || !verify_enc(Z, ped, P.F, k.pub, rv, E, C, q, Z.ell, "t", 0));
        }
        for (size_t j = 0; j < fa.size(); j++) {
            std::vector<mpz_class> t = fa;
            t[j] += 1;
            AffProof q;
            aff_tamper = aff_tamper && (!parse_aff(t, 0, q) || !verify_aff(Z, ped, P.F, k, rv, E, Ep, Cb, Cq, q, "t", 1));
        }
        check(enc_tamper && aff_tamper, "pi^enc, pi^aff: changing any one of the 8 / 14 proof values makes the proof fail");
        // A dishonest B with nu just over the bound of pi^aff, or b just over its bound, is rejected.
        mpz_class nu_big = pow2k(Z.nu + Z.slack) + 1, b_big = b + P.p() * pow2k(Z.ell + Z.slack);
        mpz_class Ep2 = k.pub.add(k.pub.scale(E, b), k.pub.enc(nu_big, xi)), Cq2b = ped.commit(P.F.neg(nu_big), rq);
        AffProof pa2 = prove_aff(Z, ped, P.F, k.pub, rps.pub, E, Ep2, Cb, Cq2b, b, nu_big, xi, rb, rq, "t", 1);
        mpz_class Ep3 = k.pub.add(k.pub.scale(E, b_big), k.pub.enc(nu, xi));
        AffProof pa3 = prove_aff(Z, ped, P.F, k.pub, rps.pub, E, Ep3, Cb, Cq, b_big, nu, xi, rb, rq, "t", 1);
        check(!verify_aff(Z, ped, P.F, k, rv, E, Ep2, Cb, Cq2b, pa2, "t", 1) && !verify_aff(Z, ped, P.F, k, rv, E, Ep3, Cb, Cq, pa3, "t", 1),
              "pi^aff: a mask nu or an input b beyond the proven range is rejected");
    }

    std::cout << "== soundness: tampered or false Sigma proofs are rejected ==" << std::endl;
    {
        Party pc(0, P);   // parameters only, no network
        const Field& F = P.F;
        // product proof
        mpz_class a = F.rand(), b = F.rand(), ra = F.rand(), rb = F.rand(), rc = F.rand();
        mpz_class A = ped.commit(a, ra), B = ped.commit(b, rb), Cc = ped.commit(F.mul(a, b), rc);
        ProductProof pp = prove_product(pc, "t", 0, A, B, Cc, a, ra, F.sub(rc, F.mul(a, rb)));
        bool prod_ok = verify_product(pc, "t", 0, A, B, Cc, pp);
        mpz_class Cw = ped.commit(F.add(F.mul(a, b), 1), rc);
        ProductProof pw_ = prove_product(pc, "t", 0, A, B, Cw, a, ra, F.sub(rc, F.mul(a, rb)));
        bool prod_false = !verify_product(pc, "t", 0, A, B, Cw, pw_);
        bool prod_tamper = true;
        for (mpz_class ProductProof::*f : {&ProductProof::T1, &ProductProof::T2, &ProductProof::za, &ProductProof::zr, &ProductProof::zs}) {
            ProductProof q = pp;
            q.*f = (q.*f + 1) % (f == &ProductProof::T1 || f == &ProductProof::T2 ? P.Q : P.p());
            prod_tamper = prod_tamper && !verify_product(pc, "t", 0, A, B, Cc, q);
        }
        check(prod_ok && prod_false && prod_tamper && !verify_product(pc, "u", 0, A, B, Cc, pp) && !verify_product(pc, "t", 1, A, B, Cc, pp),
              "pi^prod: accepts a true product; rejects a false one, any changed value, another context or prover");
        // power proof for e = d, with the product proof n * t = 1 of a write mask, checked in a batch
        mpz_class m = F.rand_nonzero(), rm = F.rand(), t = F.pow(m, P.d), rt = F.rand(), nn = F.inv(t), rn = F.rand();
        mpz_class Cm = ped.commit(m, rm), Ct = ped.commit(t, rt), Cn = ped.commit(nn, rn);
        PowerProof pw = prove_power(pc, "pw", 0, P.d, Cm, m, rm, Ct, rt);
        ProductProof pr = prove_product(pc, "pr", 0, Cn, Ct, P.g, nn, rn, F.neg(F.mul(nn, rt)));
        bool pow_ok = verify_power(pc, "pw", 0, P.d, Cm, Ct, pw) &&
                      verify_power_and_product_batch(pc, "pw", "pr", 0, P.d, Cm, Ct, pw, Cn, Ct, P.g, pr);
        mpz_class Ct2 = ped.commit(F.add(t, 1), rt);
        PowerProof pw2 = prove_power(pc, "pw", 0, P.d, Cm, m, rm, Ct2, rt);
        bool pow_false = !verify_power(pc, "pw", 0, P.d, Cm, Ct2, pw2) &&
                         !verify_power_and_product_batch(pc, "pw", "pr", 0, P.d, Cm, Ct2, pw2, Cn, Ct2, P.g, pr);
        // one changed value anywhere among the ~380 steps must make the batch check fail
        bool batch_tamper = true;
        for (int trial = 0; trial < 12 && batch_tamper; trial++) {
            PowerProof q = pw;
            const size_t j = rand_below(q.steps.size()).get_ui();
            if (trial % 3 == 0) q.steps[j].za = (q.steps[j].za + 1) % P.p();
            else if (trial % 3 == 1) q.steps[j].T2 = ped.mul(q.steps[j].T2, P.g);
            else q.inter[j % q.inter.size()] = ped.mul(q.inter[j % q.inter.size()], P.h);
            batch_tamper = !verify_power_and_product_batch(pc, "pw", "pr", 0, P.d, Cm, Ct, q, Cn, Ct, P.g, pr);
        }
        ProductProof prb = pr;
        prb.zs = (prb.zs + 1) % P.p();
        mpz_class nbad = ped.commit(F.add(nn, 1), rn);   // n * t != 1
        ProductProof prn = prove_product(pc, "pr", 0, nbad, Ct, P.g, F.add(nn, 1), rn, F.neg(F.mul(F.add(nn, 1), rt)));
        check(pow_ok && pow_false && batch_tamper && !verify_power_and_product_batch(pc, "pw", "pr", 0, P.d, Cm, Ct, pw, Cn, Ct, P.g, prb) &&
                  !verify_power_and_product_batch(pc, "pw", "pr", 0, P.d, Cm, Ct, pw, nbad, Ct, P.g, prn),
              "pi^pow (380 steps, batch check): accepts m^d; rejects a wrong power, a changed step, or n * t != 1");
        // an element of order outside G (a square mod Q of order != p) is caught by the batch membership test
        mpz_class z;
        do { z = rand_below(P.Q - 2) + 2; z = z * z % P.Q; } while (ped.in_group(z));
        check(ped.is_qr(z) && !ped.all_in_group({ped.commit(1, 2), P.g, z, P.h}) && ped.all_in_group({ped.commit(1, 2), P.g, Cm, P.h}),
              "commitments outside G are rejected (batch membership test)");
        // a random "opening proof" does not verify; a valid one does not move to another prover
        OpeningProof op = prove_opening(pc, "o", 0, Cm, m, rm);
        OpeningProof fake{ped.commit(F.rand(), F.rand()), F.rand(), F.rand()};
        check(verify_opening(pc, "o", 0, Cm, op) && !verify_opening(pc, "o", 1, Cm, op) && !verify_opening(pc, "o", 0, Cm, fake) &&
                  !verify_opening(pc, "o", 0, Ct, op),
              "proof of knowledge: accepts the real one; rejects a forgery, another prover or another commitment");
    }

    std::cout << "== delivery in E and complaints (writer side) ==" << std::endl;
    {
        EcKeyPair k = EcKeyPair::generate();
        mpz_class w0 = P.F.rand(), r0 = P.F.rand(), w1 = P.F.rand(), r1 = P.F.rand();
        mpz_class Cw[2] = {ped.commit(w0, r0), ped.commit(w1, r1)};
        HCipher ct[2] = {hel_encrypt(P.F, kDeliveryLabel, k.pk, {77, w0, r0}), hel_encrypt(P.F, kDeliveryLabel, k.pk, {77, w1, r1})};
        DeliveryResult d = delivery_finish(P, ped, k, 77, ct, Cw);
        check(d.ok && d.wp == P.F.pow(P.F.add(w0, w1), P.d), "hashed ElGamal round trip; WP = omega^d");
        HComplaint c0 = make_complaint("ct0", k, ct[0]);
        check(!delivery_judge(P, ped, "ct0", k.pk, 77, Cw[0], ct[0], c0), "complaint about a correct share is rejected");
        HCipher bad = hel_encrypt(P.F, kDeliveryLabel, k.pk, {77, P.F.add(w1, 1), r1});
        HCipher both_bad[2] = {ct[0], bad};
        DeliveryResult d2 = delivery_finish(P, ped, k, 77, both_bad, Cw);
        HComplaint cb = make_complaint("ct1", k, bad);
        check(!d2.ok && d2.bad_server == 1 && delivery_judge(P, ped, "ct1", k.pk, 77, Cw[1], bad, cb),
              "a wrong share is found by the writer and the complaint blames the server");
        HCipher wrong_idx = hel_encrypt(P.F, kDeliveryLabel, k.pk, {78, w1, r1});
        HComplaint ci = make_complaint("ct2", k, wrong_idx);
        check(delivery_judge(P, ped, "ct2", k.pk, 77, Cw[1], wrong_idx, ci), "a wrong index in the ciphertext blames the server");
        HComplaint lie = cb;
        EcPoint g1 = ec_base(ec_random_scalar());
        ec_add(lie.K, lie.K, g1);
        check(!delivery_judge(P, ped, "ct1", k.pk, 77, Cw[1], bad, lie), "a writer cannot frame a server with a false K");
        check(!delivery_judge(P, ped, "other", k.pk, 77, Cw[1], bad, cb), "a complaint cannot be moved to another ciphertext");
        HComplaint lie2 = c0;
        lie2.proof.z = ec_random_scalar();
        check(!delivery_judge(P, ped, "ct0", k.pk, 77, Cw[0], ct[0], lie2) && !verify_dleq("ct0", k.pk, ct[0].R, c0.K, lie2.proof),
              "a forged Chaum-Pedersen proof is rejected");
        HCipher noR = ct[0];
        noR.R = EcPoint{};
        check(delivery_judge(P, ped, "ct0", k.pk, 77, Cw[0], noR, c0), "a ciphertext with an invalid R is the sender's fault");
        EcScalar s1 = ec_random_scalar(), s2, s3, big;
        big.fill(0xff);                                    // far above the group order
        check(ec_scalar_from_mpz(ec_scalar_to_mpz(s1), s2) && s2 == s1 && !ec_scalar_from_mpz(ec_scalar_to_mpz(big), s3),
              "scalars: exact round trip through the log, non-canonical scalars rejected");
    }

    std::cout << "== log and checked opening ==" << std::endl;
    auto f_open = [](Party& c) -> std::string {
        std::vector<CShared> s = share_fresh(c, {mpz_class(c.me + 5)}, "t.share");
        mpz_class v = open1(c, s[0], "t.open");
        CShared t = add_const(c, scale(c, s[0], 3), 4);                 // local ops: 3x + 4
        return v == 11 && open1(c, t, "t.open2") == 37 ? "ok" : "wrong value";
    };
    Outcome o = two_party(P, dir, f_open);
    check(o.ok(), "honest opening and local operations", both(o));
    o = two_party(P, dir, f_open, 1, {"t.open"});
    check(blames(o.s[0], 1), "wrong opening by S1 is blamed by S0", both(o));
    o = two_party(P, dir, f_open, 0, {"t.open2"});
    check(blames(o.s[1], 0), "wrong opening of a derived value by S0 is blamed by S1", both(o));

    std::cout << "== preprocessing: triples, masks, zeta (MtA with pi^enc, pi^aff, pi^dec) ==" << std::endl;
    auto f_tri = [](Party& c) -> std::string {
        PrepBatch B = gen_prep_batch(c, 3, 0, 0, 0);
        if (B.triples.size() != 3) return "wrong number of triples";
        for (const Triple& t : B.triples) {
            mpz_class a = open1(c, t.a, "chk.a"), b = open1(c, t.b, "chk.b"), g = open1(c, t.c, "chk.c");
            if (c.F().mul(a, b) != g) return "gamma != alpha * beta";
        }
        std::vector<CShared> x = share_fresh(c, {c.F().rand()}, "x"), y = share_fresh(c, {c.F().rand()}, "y");
        CShared z = mul_many(c, x, y, {B.triples[0]}, "mul")[0];
        mpz_class xv = open1(c, x[0], "ox"), yv = open1(c, y[0], "oy");
        return open1(c, z, "oz") == c.F().mul(xv, yv) ? "ok" : "Mult wrong";
    };
    o = two_party(P, dir, f_tri);
    check(o.ok(), "honest triples are correct, Mult works", both(o));
    // every deviation is caught by a proof in the round it happens; the honest server blames
    for (const char* pt : {"triple.psi", "mta.E", "mta.big", "mta.nu", "mta.Ep", "mta.dec"})
        for (int cs = 0; cs < 2; cs++) {
            o = two_party(P, dir, f_tri, cs, {pt});
            check(blames(o.s[1 - cs], cs), std::string("S") + std::to_string(1 - cs) + " blames S" + std::to_string(cs) +
                                               " for " + pt, both(o));
        }

    std::cout << "== masks, key, MPC-PRF in shared mode ==" << std::endl;
    auto f_prf = [](Party& c) -> std::string {
        PrepBatch B = gen_prep_batch(c, 2, 2, 1, 1);
        for (const WriteMask& w : B.wmasks) {
            mpz_class r = open1(c, w.r, "chk.r"), v = open1(c, w.v, "chk.v");
            if (r == 0 || c.F().mul(v, c.F().pow(r, c.P.d)) != 1) return "write mask: v != r^-d";
        }
        mpz_class r = open1(c, B.smasks[0], "chk.s");
        if (r == 0 || c.F().pow(r, c.P.d) != 1) return "setup mask: r^d != 1";
        if (B.zetas.size() != 1) return "no zeta";
        CShared key = register_key(c);
        std::vector<CShared> xs = {share_public(c, 7), share_fresh(c, {mpz_class(c.me ? 3 : 9)}, "x")[0]};
        std::vector<CShared> y = prf_shared_batch(c, key, xs, B.wmasks, B.triples);
        mpz_class k = open1(c, key, "chk.k");
        if (open1(c, y[0], "chk.y0") != c.F().pow(k + 7, c.P.d)) return "PRF on a public input wrong";
        if (open1(c, y[1], "chk.y1") != c.F().pow(k + 12, c.P.d)) return "PRF on a shared input wrong";
        return "ok";
    };
    o = two_party(P, dir, f_prf);
    check(o.ok(), "write mask v = r^-d, setup mask r^d = 1, y = (k + x)^d", both(o));
    for (int cs = 0; cs < 2; cs++) {
        o = two_party(P, dir, f_prf, cs, {"mask.pow"});
        check(blames(o.s[1 - cs], cs), "a wrong mask factor of S" + std::to_string(cs) + " fails pi^pow / pi^prod and is blamed", both(o));
    }
    o = two_party(P, dir, f_prf, 1, {"zeta.pok"});
    check(blames(o.s[0], 1), "invalid proof of knowledge for zeta is blamed", both(o));
    o = two_party(P, dir, f_prf, 0, {"key.pok"});
    check(blames(o.s[1], 0), "invalid proof of knowledge for the key share is blamed", both(o));
    o = two_party(P, dir, f_prf, 1, {"omega"});
    check(blames(o.s[0], 1), "wrong opening of omega is blamed", both(o));

    std::cout << "== MPC-PRF in delivery mode (setup) ==" << std::endl;
    {
        const EcKeyPair writer = EcKeyPair::generate();   // made before the fork, so both test processes know it
        auto f_del = [&writer](Party& c) -> std::string {
            PrepBatch B = gen_prep_batch(c, 1, 0, 1, 0);
            CShared key = register_key(c);
            const mpz_class x = 41;
            MulJob j = prf_begin(c, key, share_public(c, x), B.smasks[0], B.triples[0]);
            Opener op;
            op.add(j);
            std::vector<mpz_class> ev = op.run(c, "del.mul");
            CShared omega = mul_end(c, j, ev[0], ev[1]);
            HCipher mine = delivery_encrypt(c, writer.pk, x, omega);
            Round r = c.post1("del.ct", [&] { std::vector<mpz_class> v = {ec_to_mpz(mine.R)}; v.insert(v.end(), mine.c.begin(), mine.c.end()); return v; }());
            // the writer's side, from the log only
            HCipher ct[2];
            for (int b = 0; b < 2; b++) {
                const std::vector<mpz_class>& v = c.rec(r.first(b)).values;
                if (v.size() != 4 || !ec_from_mpz(v[0], ct[b].R)) return "malformed delivery record";
                ct[b].c = {v[1], v[2], v[3]};
            }
            DeliveryResult d = delivery_finish(c.P, c.ped, writer, x, ct, omega.C);
            if (!d.ok) {
                const int b = d.bad_server;
                HComplaint cp = make_complaint("del", writer, ct[b]);
                return delivery_judge(c.P, c.ped, "del", writer.pk, x, omega.C[b], ct[b], cp) ? "BLAME S" + std::to_string(b) + ": complaint upheld"
                                                                                              : "complaint rejected";
            }
            mpz_class k = open1(c, key, "chk.k");
            return d.wp == c.F().pow(k + x, c.P.d) ? "ok" : "writer got a wrong WP";
        };
        o = two_party(P, dir, f_del);
        check(o.ok(), "the writer gets omega^d = (k + x)^d; omega is never opened", both(o));
        o = two_party(P, dir, f_del, 0, {"setup.deliver"});
        check(blames(o.s[0], 0) && blames(o.s[1], 0), "a wrong share to the writer: the complaint blames S0", both(o));
    }

    std::cout << std::endl << "unit tests: " << g_pass << " passed, " << g_fail << " failed" << std::endl;
    return g_fail == 0 ? 0 : 1;
}
