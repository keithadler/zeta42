"""Cube-and-conquer driver: split one hard CNF across all M3 cores.

march_cu partitions the search space into cubes (partial assignments chosen
by lookahead); each cube plus the original CNF is an independent, much easier
instance for kissat.  All cubes UNSAT -> UNSAT.  Any cube SAT -> SAT (and we
stop early).  This is the Heule pipeline (Pythagorean triples, Schur 5) —
its home turf is exactly our forced-pair/colouring UNSAT instances.

    ./.venv/bin/python cnc.py instance.cnf [workers] [cutoff]

Cubes are solved in shuffled order with a process pool; progress is flushed
per batch so an interrupted run shows how far it got.
"""

import os
import random
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
MARCH = os.path.join(HERE, "external", "CnC", "march_cu", "march_cu")
KISSAT = "kissat"


def make_cubes(cnf, cutoff=None, log=print):
    cubes_path = cnf + ".cubes"
    cmd = [MARCH, cnf, "-o", cubes_path]
    if cutoff:
        cmd += ["-d", str(cutoff)]      # cut_depth: <= 2^cutoff cubes
    t0 = time.time()
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    cubes = []
    with open(cubes_path) as fh:
        for line in fh:
            if line.startswith("a "):
                lits = [int(x) for x in line.split()[1:-1]]
                cubes.append(lits)
    log(f"march_cu: {len(cubes)} cubes  [{time.time()-t0:.0f}s]")
    return cubes


def solve_cube(args):
    cnf, lits, idx = args
    # append cube literals as unit clauses; header clause-count patched
    with open(cnf) as fh:
        header = fh.readline().split()
        body = fh.read()
    nv, nc = int(header[2]), int(header[3])
    piece = f"/tmp/cnc-{os.getpid()}-{idx}.cnf"
    with open(piece, "w") as fh:
        fh.write(f"p cnf {nv} {nc + len(lits)}\n")
        fh.write(body)
        for l in lits:
            fh.write(f"{l} 0\n")
    t0 = time.time()
    r = subprocess.run([KISSAT, "-q", piece], capture_output=True, text=True)
    os.unlink(piece)
    verdict = ("SAT" if "s SATISFIABLE" in r.stdout
               else "UNSAT" if "s UNSATISFIABLE" in r.stdout else "?")
    return idx, verdict, time.time() - t0


def run(cnf, workers=6, cutoff=None):
    print(f"cube-and-conquer: {cnf}, {workers} workers", flush=True)
    t0 = time.time()
    cubes = make_cubes(cnf, cutoff)
    order = list(range(len(cubes)))
    random.Random(0).shuffle(order)

    done = unsat = 0
    sat_found = None
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(solve_cube, (cnf, cubes[i], i)): i for i in order}
        for fut in as_completed(futs):
            idx, verdict, dt = fut.result()
            done += 1
            if verdict == "UNSAT":
                unsat += 1
            elif verdict == "SAT":
                sat_found = idx
                print(f"cube {idx}: SAT — instance is SAT, stopping", flush=True)
                for f in futs:
                    f.cancel()
                break
            else:
                print(f"cube {idx}: verdict '?' — investigate", flush=True)
            if done % 25 == 0 or done == len(cubes):
                el = time.time() - t0
                eta = el / done * (len(cubes) - done)
                print(f"  {done}/{len(cubes)} cubes, {unsat} UNSAT  "
                      f"[{el:.0f}s, ~{eta:.0f}s left]", flush=True)

    total = time.time() - t0
    if sat_found is not None:
        print(f"RESULT: SAT  [{total:.0f}s]", flush=True)
    elif unsat == len(cubes):
        print(f"RESULT: UNSAT — all {len(cubes)} cubes refuted  "
              f"[{total:.0f}s total]", flush=True)
    else:
        print(f"RESULT: incomplete ({unsat}/{len(cubes)})  [{total:.0f}s]",
              flush=True)


if __name__ == "__main__":
    cnf = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    cutoff = int(sys.argv[3]) if len(sys.argv) > 3 else None
    run(cnf, workers, cutoff)
