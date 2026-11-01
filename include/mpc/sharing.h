// =====================================================================
// sharing.h -- committed sharings <x> (paper: "Committed sharing")
//
//   <x> means:  x = x_0 + x_1 mod p,
//               S_b privately holds its share x_b and randomness rho_b,
//               the public log holds C_b = Com(x_b; rho_b) for b = 0, 1.
//
// In this two-process version a CShared object lives in ONE server: x and rho
// are that server's own share and randomness, and C[0], C[1] are the two
// public commitments (both servers derive the same C[0], C[1]).
//
// Local operations need no communication: each server updates its share and
// randomness, and both derive the new commitments from the old ones, so a
// server can never commit to a wrong intermediate value.
// Checked opening: both servers post (x_b, rho_b); each checks the peer's pair
// against its commitment. A mismatch is signed evidence against the peer.
// Beaver multiplication: two checked openings plus local work.
// =====================================================================
#pragma once
#include "net/party.h"

namespace parry {

struct CShared {
    mpz_class x = 0, rho = 0;   // my share and my randomness
    mpz_class C[2] = {1, 1};    // public commitments of S_0's and S_1's share
};

// Checks that every commitment the peer posted lies in G. A value outside G could
// break the binding argument, so it counts as cheating.
inline void check_peer_commitments(Party& c, const std::vector<mpz_class>& Cs, long rec, const std::string& what) {
    if (!c.ped.all_in_group(Cs)) throw Fault(c.peer(), "commitment not in G: " + what, rec);
}

// k fresh sharings in which each server picks its own shares (mine[i] is my share
// of the i-th value) and commits to them. One round.
inline std::vector<CShared> share_fresh(Party& c, const std::vector<mpz_class>& mine, const std::string& label) {
    std::vector<CShared> out(mine.size());
    std::vector<mpz_class> Cs(mine.size());
    parallel_for(mine.size(), [&](size_t i) {
        out[i].x = c.F().mod(mine[i]);
        out[i].rho = c.F().rand();
        Cs[i] = c.ped.commit(out[i].x, out[i].rho);
    });
    Round r = c.post1(label, Cs);
    const std::vector<mpz_class>& theirs = c.peer_values(r);
    if (theirs.size() != mine.size()) throw Fault(c.peer(), "wrong number of commitments: " + label, r.first(c.peer()));
    check_peer_commitments(c, theirs, r.first(c.peer()), label);
    for (size_t i = 0; i < mine.size(); i++) {
        out[i].C[c.me] = Cs[i];
        out[i].C[c.peer()] = theirs[i];
    }
    return out;
}

// A value known only to S_owner, seen as a sharing: the other server's share is 0
// with randomness 0, whose commitment is Com(0; 0) = 1.
inline CShared share_of_one(const Party& c, int owner, const mpz_class& x, const mpz_class& rho, const mpz_class& Cx) {
    CShared s;
    if (c.me == owner) { s.x = x; s.rho = rho; }
    s.C[owner] = Cx;
    s.C[1 - owner] = 1;
    return s;
}

// A public value v: S_0 holds v, S_1 holds 0, both with randomness 0.
// So C_0 = g^v and C_1 = 1, which anyone can compute.
inline CShared share_public(const Party& c, const mpz_class& v) {
    mpz_class m = c.F().mod(v);
    return share_of_one(c, 0, m, 0, c.ped.g_pow(m));
}

// ---- Local operations ----

// <x + y>: C_{x+y,b} = C_{x,b} * C_{y,b}
inline CShared add(const Party& c, const CShared& A, const CShared& B) {
    CShared R;
    R.x = c.F().add(A.x, B.x);
    R.rho = c.F().add(A.rho, B.rho);
    for (int b = 0; b < 2; b++) R.C[b] = c.ped.mul(A.C[b], B.C[b]);
    return R;
}

// <x - y>: C_{x-y,b} = C_{x,b} / C_{y,b}
inline CShared sub(const Party& c, const CShared& A, const CShared& B) {
    CShared R;
    R.x = c.F().sub(A.x, B.x);
    R.rho = c.F().sub(A.rho, B.rho);
    for (int b = 0; b < 2; b++) R.C[b] = c.ped.mul(A.C[b], c.ped.inv(B.C[b]));
    return R;
}

// <k * x> for a public k: C_{kx,b} = C_{x,b}^k
inline CShared scale(const Party& c, const CShared& A, const mpz_class& k) {
    CShared R;
    R.x = c.F().mul(A.x, k);
    R.rho = c.F().mul(A.rho, k);
    for (int b = 0; b < 2; b++) R.C[b] = c.ped.pow(A.C[b], k);
    return R;
}

// <x + k> for a public k: only S_0 adds k, so C_{x+k,0} = C_{x,0} * g^k.
inline CShared add_const(const Party& c, const CShared& A, const mpz_class& k) {
    CShared R = A;
    if (c.me == 0) R.x = c.F().add(A.x, k);
    R.C[0] = c.ped.mul(A.C[0], c.ped.g_pow(k));
    return R;
}

// ---- Checked opening, many values in one round ----
// Each server posts (x_b, rho_b) for every item; each checks the peer's pairs
// against the peer's commitments. Returns the opened values x = x_0 + x_1.
inline std::vector<mpz_class> open_many(Party& c, const std::vector<CShared>& items, const std::string& label) {
    std::vector<mpz_class> mine;
    mine.reserve(2 * items.size());
    for (const CShared& s : items) {
        mpz_class v = s.x;
        if (c.cheat(label)) v = c.F().add(v, 1);   // test hook: wrong share
        mine.push_back(v);
        mine.push_back(s.rho);
    }
    Round r = c.post1(label, mine);
    const long rec = r.first(c.peer());
    const std::vector<mpz_class>& th = c.peer_values(r);
    if (th.size() != 2 * items.size()) throw Fault(c.peer(), "malformed opening record: " + label, rec);
    std::vector<char> bad(items.size(), 0);
    parallel_for(items.size(), [&](size_t i) {
        const mpz_class& v = th[2 * i];
        bad[i] = v < 0 || v >= c.P.p() || th[2 * i + 1] < 0 || !c.ped.opens(items[i].C[c.peer()], v, th[2 * i + 1]);
    });
    for (size_t i = 0; i < items.size(); i++)
        if (bad[i]) throw Fault(c.peer(), "opening does not match commitment: " + label + " item " + std::to_string(i), rec);
    std::vector<mpz_class> out(items.size());
    for (size_t i = 0; i < items.size(); i++) out[i] = c.F().add(items[i].x, th[2 * i]);
    return out;
}

inline mpz_class open1(Party& c, const CShared& s, const std::string& label) { return open_many(c, {s}, label)[0]; }

// ---- Beaver multiplication ----
struct Triple { CShared a, b, c; };   // (alpha, beta, gamma) with gamma = alpha * beta

// A multiplication split in two halves, so that the openings of many different
// multiplications (and other openings) can share one round.
// begin: the two values to open, eps = x - alpha and del = y - beta (uniformly random,
//        so they hide x and y).
// end:   z_b = gamma_b + eps*beta_b + del*alpha_b, and S_0 also adds eps*del.
//        The commitments of <z> are derived by both servers, so z = x*y is checked
//        without any proof.
struct MulJob { CShared eps_in, del_in; const Triple* t = nullptr; };

inline MulJob mul_begin(const Party& c, const CShared& X, const CShared& Y, const Triple& t) {
    return {sub(c, X, t.a), sub(c, Y, t.b), &t};
}

inline CShared mul_end(const Party& c, const MulJob& j, const mpz_class& eps, const mpz_class& del) {
    CShared Z = add(c, add(c, j.t->c, scale(c, j.t->b, eps)), scale(c, j.t->a, del));
    return add_const(c, Z, c.F().mul(eps, del));
}

// Collects values to open in one round.
struct Opener {
    std::vector<CShared> items;
    size_t add(const CShared& s) { items.push_back(s); return items.size() - 1; }
    size_t add(const MulJob& j) { items.push_back(j.eps_in); items.push_back(j.del_in); return items.size() - 2; }
    std::vector<mpz_class> run(Party& c, const std::string& label) { return open_many(c, items, label); }
};

// n multiplications in one round.
inline std::vector<CShared> mul_many(Party& c, const std::vector<CShared>& X, const std::vector<CShared>& Y,
                                     const std::vector<Triple>& T, const std::string& label) {
    if (X.size() != Y.size() || X.size() > T.size()) throw std::invalid_argument("mul_many: sizes");
    std::vector<MulJob> jobs(X.size());
    parallel_for(X.size(), [&](size_t i) { jobs[i] = mul_begin(c, X[i], Y[i], T[i]); });
    Opener op;
    for (const MulJob& j : jobs) op.add(j);
    std::vector<mpz_class> v = op.run(c, label);
    std::vector<CShared> Z(X.size());
    parallel_for(X.size(), [&](size_t i) { Z[i] = mul_end(c, jobs[i], v[2 * i], v[2 * i + 1]); });
    return Z;
}

// ---- saving sharings and triples (prep files, state files) ----
inline void write_cs(Writer& w, const CShared& s) { w.mpz(s.x); w.mpz(s.rho); w.mpz(s.C[0]); w.mpz(s.C[1]); }
inline CShared read_cs(Reader& r) {
    CShared s;
    s.x = r.mpz(); s.rho = r.mpz(); s.C[0] = r.mpz(); s.C[1] = r.mpz();
    return s;
}
inline void write_triple(Writer& w, const Triple& t) { write_cs(w, t.a); write_cs(w, t.b); write_cs(w, t.c); }
inline Triple read_triple(Reader& r) {
    Triple t;
    t.a = read_cs(r); t.b = read_cs(r); t.c = read_cs(r);
    return t;
}

}  // namespace parry
