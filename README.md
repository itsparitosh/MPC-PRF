# Parry: the accountable MPC-PRF

This is the MPC-PRF of Parry on its own, with everything it is built on:
commitments, committed sharing, checked opening and Beaver multiplication,
the zero-knowledge proofs, preprocessing (triples, masks, ζ), and the PRF
evaluation in shared mode and in delivery mode.

    F_k(x) = (k + x)^d mod p,   k = K_0 + K_1,   server S_b holds K_b

Two servers run as two programs. Every message a server sends is signed and
kept in a log that both servers hold. If a server cheats, the other one names
it (`BLAME S1: ...`) and saves the signed record that proves it.

`docs/IMPLEMENTATION.md` explains where each part is, where every proof comes
from, and why each optimization is safe.

## Build

You need g++ with C++17, GMP (`libgmp-dev`), OpenSSL (`libssl-dev`) and
libsodium (`libsodium-dev`).

No sudo, so libsodium cannot be installed system-wide? Then build it inside
this folder once (about 30 s; the official source release is included in
`third_party/` and its SHA-256 is checked):

```
make deps     # builds libsodium into ./deps
make          # uses ./deps automatically
```

`params.txt` holds the public parameters (p, d, Q, g, h). If it is missing, it
is generated once (about a minute).

## Run in two terminals

Open two terminals in this folder. S0 waits for S1 on port 9300; start them in
any order.

```
# terminal 1                          # terminal 2
make                                  
./mpcprf p0 init                      ./mpcprf p1 init        # once: keys and the PRF key
./mpcprf p0 run --n 8                 ./mpcprf p1 run --n 8   # n PRF calls in each mode
```

* `init` (once): signing keys; Paillier keys with π^prm, π^mod,
  π^fac; the PRF key shares with proofs of knowledge. Saved in
  `run/s<b>.state`. It prints its own timing table (INIT); this is a one-time
  cost of the system and is not part of any PRF call.
* `run --n N`: N PRF calls in shared mode and N in delivery mode, each with
  its own preprocessing, and one timing table per mode. Both servers must use
  the same options.
* `--cores C`: how many cores this server may use (default: all of them).
  Each server takes its own value, e.g. `./mpcprf p0 run --n 16 --cores 4`
  and `./mpcprf p1 run --n 16 --cores 4`. The report shows the number used.
  If both servers run on the same machine, they share its cores, so give each
  at most half of them.

What uses several cores: every preprocessing step runs over all calls of
the batch in parallel (and the ~380 product proofs inside a power proof run
in parallel too), and so do the local steps of the online part. The network rounds are
the same whatever the number of cores, so results are identical; only the
time changes.

Each server then prints one table per mode (S0, `--n 1 --cores 1`, both
servers on one 2-core VM):

```
[S0] 1 PRF call(s) per mode, 1 core(s)
----------------------------------------------------------------------------
SHARED MODE (write phase)                           total (ms) per call (ms)
  (y = (k + x)^d stays secret-shared)
OFFLINE (per call: 1 write mask + 1 Beaver triple, 4 MtA runs)
  MtA: Paillier encrypt, E', decrypt                     128.7         128.7
  Beaver triple: local alpha, beta, psi                    0.0           0.0
  mask: local m, m^e, inverse                              0.1           0.1
  Pedersen commitments (make, check in G)                 13.7          13.7
  proof pi^enc  (make + check)                           155.1         155.1
  proof pi^aff  (make + check)                           338.9         338.9
  proof pi^dec  (make + check)                           361.6         361.6
  proof pi^prod (make + check)                             8.8           8.8
  proof pi^pow  (make + check)                          1132.3        1132.3
  network: 3 rounds, incl. waiting for the peer          141.9         141.9
OFFLINE TOTAL                                           2281.1        2281.1
ONLINE (2 rounds)
  round 1: Beaver Mult, open eps and delta                 9.0           9.0
  round 2: open omega; omega^d; y = omega^d * v            4.6           4.6
ONLINE TOTAL                                              13.6          13.6
OFFLINE + ONLINE                                        2294.7        2294.7
----------------------------------------------------------------------------
DELIVERY MODE (setup phase)                         total (ms) per call (ms)
  (the writer learns (k + x)^d)
OFFLINE (per call: 1 setup mask + 1 Beaver triple, 3 MtA runs)
  MtA: Paillier encrypt, E', decrypt                      89.1          89.1
  Beaver triple: local alpha, beta, psi                    0.0           0.0
  mask: local m, m^e, inverse                              0.0           0.0
  Pedersen commitments (make, check in G)                 11.8          11.8
  proof pi^enc  (make + check)                           126.6         126.6
  proof pi^aff  (make + check)                           200.4         200.4
  proof pi^dec  (make + check)                           216.2         216.2
  proof pi^prod (make + check)                             8.5           8.5
  proof pi^pow  (make + check)                           462.0         462.0
  network: 3 rounds, incl. waiting for the peer          189.0         189.0
OFFLINE TOTAL                                           1303.7        1303.7
ONLINE (2 rounds)
  round 1: Beaver Mult, open eps and delta                12.8          12.8
  round 2: encrypt to the writer, post                     0.4           0.4
ONLINE TOTAL (servers)                                    13.2          13.2
  writer: decrypt, check both shares, omega^d              0.9           0.9
OFFLINE + ONLINE (servers)                              1316.9        1316.9
----------------------------------------------------------------------------
Times are wall-clock time on this server, with all its cores working. Each line counts
everything of that kind done for all calls (e.g. every pi^dec made and checked), and the
lines add up to the totals. per call = total / 1 (the number of PRF calls).
The check of a mask's pi^prod is batched with its pi^pow and counted there.
```

How to read it:
* The two modes are different PRF calls in Parry: shared mode is used in the
  write phase, delivery mode in the setup phase. Each has its own table.
* Every time is **wall-clock time on this server**, with all its cores
  working. Each line counts everything of that kind done for all N calls: for
  example "pi^dec" is every pi^dec this server made plus every pi^dec of the
  other server it checked. The lines add up to OFFLINE TOTAL and ONLINE TOTAL.
* **per call** = total / N, the average time of one PRF call. With `--n 1`
  the two columns are equal.
* "network" includes waiting while the other server is still computing.
* The writer's own work in delivery mode is shown below ONLINE TOTAL, since
  the writer is not a server (S0 plays the writers, so only S0 prints it).
* For the cost of one call on one core, run with `--cores 1`. With C cores
  and a large N, per call shows the amortized time with all cores working.

## Code map

| folder | phase | what it holds |
|---|---|---|
| `include/base/` | all | randomness, F_p, bytes, SHA-256 transcripts, Ed25519, public parameters |
| `include/crypto/` | all | Pedersen commitments, Paillier, Ristretto255 (libsodium), hashed ElGamal with complaints |
| `include/net/` | all | TCP link, one server with its signed public log |
| `include/zk/` | init + prep | opening, product, power proofs; π^prm, π^mod, π^fac, π^enc, π^aff (CGGMP20) and π^dec |
| `include/init/` | init | keys, proofs about the Paillier keys, PRF key registration, state file |
| `include/mpc/` | online | committed sharing ⟨x⟩, checked opening, Beaver Mult |
| `include/prep/` | preprocessing | MtA, triples, write masks, setup masks and ζ in one 3-round batch |
| `include/prf/` | online | the MPC-PRF: shared mode and delivery mode |
| `src/mpcprf.cpp` | | the server program |
| `tests/` | | unit tests (`test_prf.cpp`, `harness.h`), program tests (`run_tests.sh`) |
