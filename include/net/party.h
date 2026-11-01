// =====================================================================
// party.h -- one server process (S0 or S1) and its copy of the public log
//
// Each server runs in its own process and holds ONLY its own secrets.
// What the paper calls the public log is kept as two identical copies:
// every protocol round, both servers sign their messages, swap them, check
// the peer's signature, and append S0's records and then S1's records.
// So both copies contain the same records in the same order, and any record
// that fails a public check is signed evidence against its sender.
//
// A record is signed together with the session id, the round number and its
// position in the round. So a record cannot be moved to another place, run or
// phase without breaking its signature.
// =====================================================================
#pragma once
#include "net/net.h"
#include "zk/cggmp.h"
#include <deque>
#include <set>

namespace parry {

constexpr unsigned kPaillierBits = 3072;   // Paillier modulus size (N_A)
constexpr unsigned kStatBits = 80;         // kappa: statistical hiding of the shift mask (MtA sizes are in ZkSizes)

struct Record {
    int sender = -1;
    uint64_t round = 0, pos = 0;          // where the record was posted
    std::string label;
    std::vector<mpz_class> values;
    Bytes sig;
};

inline std::string hex(const Bytes& b, size_t max = 64) {
    static const char* d = "0123456789abcdef";
    std::string s;
    for (size_t i = 0; i < b.size() && i < max; i++) { s += d[b[i] >> 4]; s += d[b[i] & 15]; }
    return s;
}

// The bytes a server signs for one record.
inline Bytes record_encoding(const Bytes& session, const Record& r) {
    Writer w;
    w.str("parry-log-v1");
    w.bytes(session);
    w.u64(r.round);
    w.u64(static_cast<uint64_t>(r.sender));
    w.u64(r.pos);
    w.str(r.label);
    w.mpzs(r.values);
    return w.buf;
}

inline bool record_authentic(const Bytes& session, const Record& r, const Bytes& sender_pub) {
    return ed25519_verify(sender_pub, record_encoding(session, r), r.sig);
}

class Log {
public:
    void reset(const Bytes& session) { session_ = session; recs_.clear(); }
    long append(Record r) { recs_.push_back(std::move(r)); return static_cast<long>(recs_.size()) - 1; }
    const Record& at(long id) const { return recs_.at(static_cast<size_t>(id)); }
    size_t size() const { return recs_.size(); }
    const Bytes& session() const { return session_; }

    // Saves the whole log, e.g. as evidence after a blame.
    void save(const std::string& path) const {
        Writer w;
        w.bytes(session_);
        w.u64(recs_.size());
        for (const Record& r : recs_) {
            w.u64(static_cast<uint64_t>(r.sender)); w.u64(r.round); w.u64(r.pos);
            w.str(r.label); w.mpzs(r.values); w.bytes(r.sig);
        }
        save_file(path, w.buf);
    }

private:
    Bytes session_;
    std::deque<Record> recs_;
};

// TEST HOOK ONLY. A server started with --cheat <point> deviates at every step
// whose label ends with <point>. Honest runs never set it. The tests use it to
// check that each kind of cheating is caught and blamed on the right server.
struct Adversary {
    std::set<std::string> points;
    bool at(const std::string& label) const {
        for (const std::string& pt : points)
            if (label.size() >= pt.size() && label.compare(label.size() - pt.size(), pt.size(), pt) == 0) return true;
        return false;
    }
};

// Record ids of one round: ids[b] are the records of S_b, in order.
struct Round {
    std::vector<long> ids[2];
    long first(int b) const { return ids[b].at(0); }
};

class Party {
public:
    const int me;                 // 0 or 1
    const Params& P;
    Pedersen ped;
    SigningKey sig;               // my signing key
    Bytes sig_pub[2];             // both servers' public signing keys
    PaillierKey pal_key;          // my Paillier key (secret)
    PaillierPub pal[2];           // both public Paillier keys; pal[me] = pal_key.pub
    RingPedersen rp[2];           // both ring-Pedersen parameters (modulus = that server's Paillier N)
    mpz_class rp_lam;             // my ring-Pedersen trapdoor lambda = log_t s (secret)
    ZkSizes Z;                    // sizes of the range proofs
    Channel* ch = nullptr;
    Log log;
    Adversary adv;

    Party(int b, const Params& prm) : me(b), P(prm), ped(prm), Z(prm) {}
    Party(const Party&) = delete;
    Party& operator=(const Party&) = delete;

    int peer() const { return 1 - me; }

    // My ring-Pedersen parameters with the trapdoor, for checking the peer's range proofs.
    RpVerifier rp_verifier() const { return RpVerifier{rp[me], rp_lam, &pal_key}; }
    const Field& F() const { return P.F; }
    bool cheat(const std::string& label) const { return adv.at(label); }

    // Unique name for one protocol instance. Both servers create names in the same
    // order, so they get the same names. The session id makes them unique across runs,
    // which also binds every Fiat-Shamir proof to this run.
    std::string fresh_id(const std::string& what) {
        return hex(log.session(), 8) + "/" + what + "#" + std::to_string(ctr_++);
    }

    // Starts a run: both servers agree on a fresh session id and check that the
    // peer is the server they set up with at init time.
    void start_session(const std::string& phase) {
        unsigned char nonce[16];
        random_bytes(nonce, 16);
        Writer w;
        w.str(phase);
        w.raw(nonce, 16);
        w.bytes(sig_pub[me]);
        Bytes theirs = ch->exchange(w.buf);
        Reader r(theirs);
        std::string phase_peer = r.str();
        unsigned char nonce_peer[16];
        r.raw(nonce_peer, 16);
        Bytes pub_peer = r.bytes();
        if (phase_peer != phase)
            throw std::runtime_error("the peer runs '" + phase_peer + "' but we run '" + phase + "'");
        if (!sig_pub[peer()].empty() && pub_peer != sig_pub[peer()])
            throw std::runtime_error("the peer's signing key differs from the one stored at init");
        sig_pub[peer()] = pub_peer;
        Writer s;
        s.str("parry-session");
        s.str(phase);
        const unsigned char* n0 = me == 0 ? nonce : nonce_peer;
        const unsigned char* n1 = me == 0 ? nonce_peer : nonce;
        s.raw(n0, 16);
        s.raw(n1, 16);
        log.reset(sha256(s.buf));
        round_ = 0;
        ctr_ = 0;
    }

    // One round of the public log: I post `mine` (one record per entry), the peer
    // posts its records, both append S0's records and then S1's.
    // If expect_peer is given, the peer must post exactly that many records.
    Round post(const std::string& label, const std::vector<std::vector<mpz_class>>& mine,
               size_t expect_peer = static_cast<size_t>(-1)) {
        const uint64_t rnd = round_++;
        std::vector<Record> my_recs;
        Writer w;
        w.u64(mine.size());
        for (size_t i = 0; i < mine.size(); i++) {
            Record r;
            r.sender = me; r.round = rnd; r.pos = i; r.label = label; r.values = mine[i];
            r.sig = sig.sign(record_encoding(log.session(), r));
            w.str(r.label); w.mpzs(r.values); w.bytes(r.sig);
            my_recs.push_back(std::move(r));
        }
        Bytes theirs = ch->exchange(w.buf);
        std::vector<Record> peer_recs;
        try {
            Reader rd(theirs);
            uint64_t k = rd.u64();
            if (k > theirs.size()) throw std::runtime_error("bad count");
            for (uint64_t i = 0; i < k; i++) {
                Record r;
                r.sender = peer(); r.round = rnd; r.pos = i;
                r.label = rd.str(); r.values = rd.mpzs(); r.sig = rd.bytes();
                peer_recs.push_back(std::move(r));
            }
            if (!rd.done()) throw std::runtime_error("trailing bytes");
        } catch (const std::runtime_error&) {
            throw Fault(peer(), "unreadable message in round '" + label + "'", -1);
        }
        // An unsigned or misplaced message is not evidence against anyone, so it
        // counts as the peer not sending the required message (record -1).
        for (const Record& r : peer_recs) {
            if (!record_authentic(log.session(), r, sig_pub[peer()]))
                throw Fault(peer(), "bad signature in round '" + label + "'", -1);
            if (r.label != label)
                throw Fault(peer(), "sent '" + r.label + "' when '" + label + "' was due", -1);
        }
        if (expect_peer != static_cast<size_t>(-1) && peer_recs.size() != expect_peer)
            throw Fault(peer(), "wrong number of records in round '" + label + "'", -1);
        Round out;
        auto add = [&](std::vector<Record>& recs, int b) {
            for (Record& r : recs) out.ids[b].push_back(log.append(std::move(r)));
        };
        if (me == 0) { add(my_recs, 0); add(peer_recs, 1); }
        else { add(peer_recs, 0); add(my_recs, 1); }
        return out;
    }

    // Shorthand: one record per server.
    Round post1(const std::string& label, const std::vector<mpz_class>& mine) { return post(label, {mine}, 1); }

    const Record& rec(long id) const { return log.at(id); }
    const std::vector<mpz_class>& peer_values(const Round& r, size_t i = 0) const {
        return log.at(r.ids[peer()].at(i)).values;
    }
    const std::vector<mpz_class>& my_values(const Round& r, size_t i = 0) const {
        return log.at(r.ids[me].at(i)).values;
    }

    // Raw bytes swap outside the log (used only for large data that is bound by
    // a signed hash on the log first, e.g. a database share being published).
    Bytes swap_raw(const Bytes& mine) { return ch->exchange(mine); }

private:
    uint64_t round_ = 0, ctr_ = 0;
};

}  // namespace parry
