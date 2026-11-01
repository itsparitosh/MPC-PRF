// =====================================================================
// mpcprf.cpp -- one server of the accountable MPC-PRF. Run S0 and S1 in two
// terminals (S0 listens on --port, S1 connects to --host:--port):
//
//   ./mpcprf p0 init                 ./mpcprf p1 init
//   ./mpcprf p0 run --n 4            ./mpcprf p1 run --n 4
//   ./mpcprf p0 run --n 16 --reps 3  ./mpcprf p1 run --n 16 --reps 3
//
// init  keys, proofs about the Paillier keys (pi^prm, pi^mod, pi^fac), and the
//       PRF key k = K_0 + K_1 (each share with a proof of knowledge).
// run   preprocessing for n evaluations in each mode (one 3-round batch per mode),
//       then n evaluations in shared mode and n in delivery mode, and a timing
//       report: every round of preprocessing and of the online part, total and
//       per call. --reps R repeats the whole run R times and averages.
//       --check opens the key at the end and checks every output (demo only:
//       a real system never opens k).
// A cheating server (--cheat <point>, test hook) is blamed by the honest one,
// which saves the signed log as evidence.
// =====================================================================
#include "init/keys.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>

using namespace parry;
using clk = std::chrono::steady_clock;

namespace {

struct Options {
    int me = -1;
    std::string cmd, dir = "run", host = "127.0.0.1";
    uint16_t port = 9300;
    unsigned deadline = 3600;     // seconds a server waits for the peer in one round
    size_t n = 4, reps = 1;
    std::string mode = "both";
    bool check = false;
    std::set<std::string> cheats;
};

[[noreturn]] void usage() {
    std::cerr << "usage: ./mpcprf <p0|p1> <init|run> [--n N] [--mode both|shared|delivery] [--reps R] [--check]\n"
                 "common: --cores C (default: all cores)  --dir run  --host 127.0.0.1  --port 9300  --deadline S (default 3600)  --cheat <point> (tests only)\n";
    std::exit(2);
}

Options parse(int argc, char** argv) {
    if (argc < 3) usage();
    Options o;
    std::string r = argv[1];
    if (r == "p0") o.me = 0;
    else if (r == "p1") o.me = 1;
    else usage();
    o.cmd = argv[2];
    if (o.cmd != "init" && o.cmd != "run") usage();
    for (int i = 3; i < argc; i++) {
        std::string a = argv[i];
        auto val = [&]() -> std::string { if (i + 1 >= argc) usage(); return argv[++i]; };
        if (a == "--dir") o.dir = val();
        else if (a == "--host") o.host = val();
        else if (a == "--deadline") o.deadline = static_cast<unsigned>(std::stoul(val()));
        else if (a == "--port") o.port = static_cast<uint16_t>(std::stoi(val()));
        else if (a == "--n") o.n = std::max<size_t>(1, std::stoul(val()));
        else if (a == "--reps") o.reps = std::max<size_t>(1, std::stoul(val()));
        else if (a == "--cores" || a == "--core" || a == "--threads") num_threads() = static_cast<unsigned>(std::max(1, std::stoi(val())));
        else if (a == "--cheat") o.cheats.insert(val());
        else if (a == "--check") o.check = true;
        else if (a == "--mode") {
            o.mode = val();
            if (o.mode != "both" && o.mode != "shared" && o.mode != "delivery") usage();
        }
        else {
            std::cerr << "unknown option: " << a << "\n";
            usage();
        }
    }
    return o;
}

double secs_since(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

std::string kib(uint64_t b) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(1) << b / 1024.0 << " KiB";
    return s.str();
}

// ---- the timing report ----
// One row: time on this server, and how many PRF calls used it (for the average).
struct Row {
    std::string name;
    bool heading = false;
    double secs = 0;
    size_t calls = 1;
    bool total = false;
};

using Rows = std::vector<Row>;

void heading(Rows& R, const std::string& name) { R.push_back({name, true, 0, 1, false}); }
void row(Rows& R, const std::string& name, double secs, size_t calls, bool total = false) {
    R.push_back({name, false, secs, calls, total});
}

// Prints rows averaged over `reps` runs: total time and time per PRF call.
void print_report(const Rows& R, size_t n, size_t reps, int me) {
    const double div = static_cast<double>(reps);
    const std::string rule(76, '-');
    std::cout << "\n[S" << me << "] " << n << " PRF call(s) per mode, " << num_threads() << " core(s)"
              << (reps > 1 ? ", average of " + std::to_string(reps) + " runs" : "") << "\n";
    for (const Row& r : R) {
        if (r.heading) {
            if (r.name.empty()) continue;
            if (r.name.find("MODE") != std::string::npos)
                std::cout << rule << "\n" << std::left << std::setw(48) << r.name.substr(0, r.name.find(':')) << std::right
                          << std::setw(14) << "total (ms)" << std::setw(14) << "per call (ms)" << "\n"
                          << "  (" << r.name.substr(r.name.find(':') + 2) << ")\n";
            else
                std::cout << r.name << "\n";
            continue;
        }
        const double t = r.secs / div * 1000;
        std::cout << std::left << std::setw(48) << ((r.total ? "" : "  ") + r.name) << std::right << std::fixed
                  << std::setprecision(1) << std::setw(14) << t << std::setw(14) << t / static_cast<double>(r.calls) << "\n";
    }
    std::cout << rule << "\n"
              << "Times are wall-clock time on this server, with all its cores working. Each line counts\n"
              << "everything of that kind done for all calls (e.g. every pi^dec made and checked), and the\n"
              << "lines add up to the totals. per call = total / " << n << " (the number of PRF calls).\n"
              << "The check of a mask's pi^prod is batched with its pi^pow and counted there.\n";
}

// Measures consecutive steps: each lap() returns the time since the previous one.
struct Meter {
    clk::time_point t0 = clk::now();
    void reset() { t0 = clk::now(); }
    std::pair<double, uint64_t> lap() {
        std::pair<double, uint64_t> r{secs_since(t0), 0};
        reset();
        return r;
    }
};


// One complete measurement: preprocessing for n calls in each mode, then n calls in each mode.
// The rows of one mode's preprocessing batch.
void offline_rows(Rows& R, const PrepStats& ps, size_t n) {
    row(R, "MtA: Paillier encrypt, E', decrypt", ps.get(Op::MtA), n);
    row(R, "Beaver triple: local alpha, beta, psi", ps.get(Op::TripleLocal), n);
    row(R, "mask: local m, m^e, inverse", ps.get(Op::MaskLocal), n);
    row(R, "Pedersen commitments (make, check in G)", ps.get(Op::Commit), n);
    row(R, "proof pi^enc  (make + check)", ps.get(Op::ProofEnc), n);
    row(R, "proof pi^aff  (make + check)", ps.get(Op::ProofAff), n);
    row(R, "proof pi^dec  (make + check)", ps.get(Op::ProofDec), n);
    row(R, "proof pi^prod (make + check)", ps.get(Op::ProofProd), n);
    row(R, "proof pi^pow  (make + check)", ps.get(Op::ProofPow), n);
    row(R, "network: 3 rounds, incl. waiting for the peer", ps.get(Op::Network), n);
    row(R, "OFFLINE TOTAL", ps.total(), n, true);
}

// One complete measurement of the chosen mode(s): preprocessing for n calls, then n calls.
Rows run_once(Party& c, const ServerState& st, size_t n, bool check, bool do_shared, bool do_delivery) {
    const Field& F = c.F();
    Rows R;
    Meter m;
    mpz_class k_open;                    // only for --check
    std::vector<mpz_class> xv, yv, wp(n);

    if (do_shared) {
        // ---- offline: per call 1 write mask (2 MtA runs) + 1 Beaver triple (2 MtA runs) ----
        heading(R, "SHARED MODE (write phase): y = (k + x)^d stays secret-shared");
        heading(R, "OFFLINE (per call: 1 write mask + 1 Beaver triple, 4 MtA runs)");
        PrepStats ps;
        PrepBatch B = gen_prep_batch(c, n, n, 0, 0, &ps);
        offline_rows(R, ps, n);
        // inputs: in Parry they come from the writer; here each server picks random shares (not timed)
        std::vector<mpz_class> mine(n);
        for (mpz_class& x : mine) x = F.rand();
        std::vector<CShared> xs = share_fresh(c, mine, c.fresh_id("inputs"));
        // ---- online ----
        heading(R, "ONLINE (2 rounds)");
        m.reset();
        const std::string id = c.fresh_id("prf");
        std::vector<MulJob> jobs(n);
        for (size_t k = 0; k < n; k++) jobs[k] = prf_begin(c, st.key, xs[k], B.wmasks[k].r, B.triples[k]);   // s = k + x
        Opener op;
        for (const MulJob& j : jobs) op.add(j);
        std::vector<mpz_class> ev = op.run(c, id + ".mul");
        std::vector<CShared> omega(n);
        for (size_t k = 0; k < n; k++) omega[k] = mul_end(c, jobs[k], ev[2 * k], ev[2 * k + 1]);   // omega = r (k + x)
        const double t1 = m.lap().first;
        std::vector<mpz_class> w = open_many(c, omega, id + ".omega");
        std::vector<CShared> y(n);
        for (size_t k = 0; k < n; k++) {
            if (w[k] == 0) throw std::runtime_error("omega = 0 (negligible probability)");
            y[k] = prf_shared_end(c, B.wmasks[k].v, w[k]);                                           // y = omega^d v
        }
        const double t2 = m.lap().first;
        row(R, "round 1: Beaver Mult, open eps and delta", t1, n);
        row(R, "round 2: open omega; omega^d; y = omega^d * v", t2, n);
        row(R, "ONLINE TOTAL", t1 + t2, n, true);
        heading(R, "");
        row(R, "OFFLINE + ONLINE", ps.total() + t1 + t2, n, true);
        if (check) {
            k_open = open1(c, st.key, "demo.check.k");
            xv = open_many(c, xs, "demo.check.x");
            yv = open_many(c, y, "demo.check.y");
        }
    }

    if (do_delivery) {
        // ---- offline: per call 1 setup mask (1 MtA run) + 1 Beaver triple (2 MtA runs) ----
        heading(R, "DELIVERY MODE (setup phase): the writer learns (k + x)^d");
        heading(R, "OFFLINE (per call: 1 setup mask + 1 Beaver triple, 3 MtA runs)");
        PrepStats ps;
        PrepBatch B = gen_prep_batch(c, n, 0, n, 0, &ps);
        offline_rows(R, ps, n);
        // writer keys: in Parry they come from SETUPDB; here S0 makes them and posts them (not timed)
        std::vector<EcKeyPair> writers(c.me == 0 ? n : 0);
        std::vector<mpz_class> pks;
        for (EcKeyPair& wk : writers) { wk = EcKeyPair::generate(); pks.push_back(ec_to_mpz(wk.pk)); }
        Round rk = c.post1("demo.writer_keys", pks);
        const std::vector<mpz_class>& pk_list = c.rec(rk.first(0)).values;
        if (pk_list.size() != n) throw Fault(0, "demo: wrong number of writer keys", rk.first(0));
        // ---- online ----
        heading(R, "ONLINE (2 rounds)");
        m.reset();
        std::vector<MulJob> dj(n);
        for (size_t k = 0; k < n; k++) dj[k] = prf_begin(c, st.key, share_public(c, k), B.smasks[k], B.triples[k]);   // x = k
        Opener od;
        for (const MulJob& j : dj) od.add(j);
        std::vector<mpz_class> dv = od.run(c, "delivery.mul");
        std::vector<CShared> om(n);
        for (size_t k = 0; k < n; k++) om[k] = mul_end(c, dj[k], dv[2 * k], dv[2 * k + 1]);   // omega is NOT opened
        const double t1 = m.lap().first;
        std::vector<std::vector<mpz_class>> recs(n);
        for (size_t k = 0; k < n; k++) {
            EcPoint pk;
            if (!ec_from_mpz(pk_list[k], pk) || !ec_valid(pk)) throw std::runtime_error("demo: invalid writer key");
            HCipher ct = delivery_encrypt(c, pk, k, om[k]);
            recs[k] = {ec_to_mpz(ct.R), ct.c[0], ct.c[1], ct.c[2]};
        }
        Round rd = c.post("delivery.ct", recs, n);
        const double t2 = m.lap().first;
        // the writers (simulated by S0): decrypt both shares, check them, omega^d
        if (c.me == 0) {
            for (size_t k = 0; k < n; k++) {
                HCipher ct[2];
                for (int b = 0; b < 2; b++) {
                    const std::vector<mpz_class>& v = c.rec(rd.ids[b][k]).values;
                    if (v.size() != 4 || !ec_from_mpz(v[0], ct[b].R)) throw Fault(b, "malformed delivery record", rd.ids[b][k]);
                    ct[b].c = {v[1], v[2], v[3]};
                }
                DeliveryResult d = delivery_finish(c.P, c.ped, writers[k], k, ct, om[k].C);
                if (d.bad_server >= 0) {
                    const int b = d.bad_server;
                    const std::string bind = hex(c.log.session()) + "/delivery/" + std::to_string(k);
                    const bool guilty = !ec_valid(ct[b].R) || delivery_judge(c.P, c.ped, bind, writers[k].pk, k, om[k].C[b],
                                                                            ct[b], make_complaint(bind, writers[k], ct[b]));
                    if (guilty) throw Fault(b, "writer's complaint upheld: wrong share in delivery mode", rd.ids[b][k]);
                    throw std::runtime_error("writer's complaint rejected");
                }
                if (!d.ok) throw std::runtime_error("omega = 0 (negligible probability)");
                wp[k] = d.wp;
            }
        }
        const double t3 = m.lap().first;
        row(R, "round 1: Beaver Mult, open eps and delta", t1, n);
        row(R, "round 2: encrypt to the writer, post", t2, n);
        row(R, "ONLINE TOTAL (servers)", t1 + t2, n, true);
        if (c.me == 0) row(R, "writer: decrypt, check both shares, omega^d", t3, n);
        heading(R, "");
        row(R, "OFFLINE + ONLINE (servers)", ps.total() + t1 + t2, n, true);
        if (check && !do_shared) k_open = open1(c, st.key, "demo.check.k");
    }

    if (check) {
        // Demo only: compare every output with (k + x)^d. A real system never opens k.
        bool ok = true;
        for (size_t j = 0; j < n; j++) {
            if (do_shared) ok = ok && yv[j] == F.pow(k_open + xv[j], c.P.d);
            if (do_delivery && c.me == 0) ok = ok && wp[j] == F.pow(k_open + j, c.P.d);
        }
        std::cout << (ok ? "CHECK OK: every output is (k + x)^d" : "CHECK FAILED") << std::endl;
        if (!ok) throw std::runtime_error("check failed");
    }
    return R;
}

}  // namespace

int main(int argc, char** argv) {
    Options o = parse(argc, argv);
    const std::string tag = "[S" + std::to_string(o.me) + "] ";
    std::system(("mkdir -p '" + o.dir + "'").c_str());
    Params P = get_params(file_exists(o.dir + "/params.txt") ? o.dir + "/params.txt" : "params.txt");

    Channel ch;
    Party c(o.me, P);
    c.adv.points = o.cheats;
    ServerState st;
    if (o.cmd == "run") {
        if (!file_exists(state_path(o.dir, o.me))) { std::cerr << tag << "no state: run init first" << std::endl; return 2; }
        st = load_state(o.dir, o.me);
        install(c, st);
    }
    std::cout << tag << o.cmd << ", " << num_threads() << " core(s), "
              << (o.me == 0 ? "listening on port " + std::to_string(o.port) : "connecting to " + o.host + ":" + std::to_string(o.port))
              << std::endl;
    ch = o.me == 0 ? Channel::listen_on(o.port) : Channel::connect_to(o.host, o.port);
    ch.set_deadline(o.deadline);
    c.ch = &ch;
    const auto t_start = clk::now();
    const std::string log_file = o.dir + "/s" + std::to_string(o.me) + "." + o.cmd + ".log";
    try {
        if (o.cmd == "init") {
            InitStats is;
            st = run_init(c, &is);
            save_state(o.dir, st);
            std::cout << tag << "keys created, Paillier keys proven, PRF key registered\n\n"
                      << "[S" << o.me << "] INIT (once per system, not per PRF call)      time (ms)\n";
            auto ln = [](const std::string& what, double sec) {
                std::cout << "  " << std::left << std::setw(46) << what << std::right << std::fixed << std::setprecision(1)
                          << std::setw(12) << sec * 1000 << "\n";
            };
            ln("Paillier key (two safe primes), ring-Pedersen", is.keygen);
            ln("make pi^prm, pi^mod", is.prove_prm_mod);
            ln("check the peer's pi^prm, pi^mod", is.check_prm_mod);
            ln("make pi^fac", is.prove_fac);
            ln("check the peer's pi^fac", is.check_fac);
            ln("PRF key share + proof of knowledge (1 round)", is.key_reg);
            ln("network: 2 rounds, incl. waiting for the peer", is.network);
            ln("INIT TOTAL", is.keygen + is.prove_prm_mod + is.check_prm_mod + is.prove_fac + is.check_fac + is.key_reg + is.network);
            std::cout << "(the safe-prime search is random, so key generation time varies from run to run)" << std::endl;
        } else {
            c.start_session("run");
            if (c.cheat("silent")) {   // test hook: a server that stops sending (it is blamed at the deadline)
                std::this_thread::sleep_for(std::chrono::seconds(o.deadline + 5));
                return 5;
            }
            Rows sum;
            for (size_t r = 0; r < o.reps; r++) {
                Rows R = run_once(c, st, o.n, o.check, o.mode != "delivery", o.mode != "shared");
                if (sum.empty()) sum = R;
                else
                    for (size_t k = 0; k < R.size(); k++) sum[k].secs += R[k].secs;
                if (o.reps > 1) std::cout << tag << "run " << r + 1 << " of " << o.reps << " done" << std::endl;
            }
            print_report(sum, o.n, o.reps, o.me);
        }
    } catch (const Fault& e) {
        ch.send_abort(e.what());
        c.log.save(log_file);
        std::cerr << tag << "BLAME S" << e.server << ": " << e.what()
                  << (e.record >= 0 ? " (log record " + std::to_string(e.record) + ")" : "") << std::endl
                  << tag << "evidence saved in " << log_file << std::endl;
        return 3;
    } catch (const PeerGone& e) {
        c.log.save(log_file);
        if (e.deadline)
            std::cerr << tag << "S" << c.peer() << " sent nothing within the round deadline (" << o.deadline
                      << " s). A missing message counts as cheating: BLAME S" << c.peer() << "." << std::endl;
        else if (e.reason.empty())
            std::cerr << tag << "S" << c.peer() << " stopped answering. Not sending a required message counts as "
                      << "cheating: BLAME S" << c.peer() << "." << std::endl;
        else
            std::cerr << tag << "S" << c.peer() << " stopped the run: " << e.reason << std::endl;
        return 4;
    } catch (const std::exception& e) {
        ch.send_abort(std::string("error: ") + e.what());
        std::cerr << tag << "error: " << e.what() << std::endl;
        return 1;
    }
    c.log.save(log_file);
    std::cout << tag << "total " << std::fixed << std::setprecision(2) << secs_since(t_start) << " s, sent "
              << kib(ch.stats.bytes_sent) << ", received " << kib(ch.stats.bytes_recv) << ", " << ch.stats.messages
              << " rounds" << std::endl;
    return 0;
}
