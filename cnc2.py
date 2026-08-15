"""Cube-and-conquer, incremental edition.

cnc.py paid kissat startup + preprocessing per cube (~5x per-core overhead,
measured).  Here each worker is ONE persistent iglucose process fed an .icnf
file -- "p inccnf" header, the full CNF once, then its share of cubes -- so
preprocessing is paid once per worker and cubes are solved incrementally
under assumptions.  This is the pipeline march_cu was designed for
(cube-glucose.sh in Heule's CnC repo), parallelised by interleaved cube
striping across workers.

    ./.venv/bin/python cnc2.py instance.cnf [workers] [cube_depth]

Every worker must report UNSAT for the instance to be UNSAT; any SAT line
means SAT (workers poll-killed on first SAT).
"""

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MARCH = os.path.join(HERE, "external", "CnC", "march_cu", "march_cu")
IGLUCOSE = os.path.join(HERE, "external", "CnC", "iglucose", "core",
                        "iglucose_release")


def run(cnf, workers=6, depth=12):
    print(f"cnc2 (incremental): {cnf}, {workers} workers, depth {depth}",
          flush=True)
    t0 = time.time()

    cubes_path = cnf + ".cubes"
    subprocess.run([MARCH, cnf, "-o", cubes_path, "-d", str(depth)],
                   capture_output=True, check=True)
    cubes = [l for l in open(cubes_path) if l.startswith("a ")]
    print(f"march_cu: {len(cubes)} cubes  [{time.time()-t0:.0f}s]", flush=True)

    with open(cnf) as fh:
        header = fh.readline()
        body = fh.read()

    procs = []
    for w in range(workers):
        share = cubes[w::workers]                 # interleaved striping
        path = f"/tmp/cnc2-{os.getpid()}-{w}.icnf"
        with open(path, "w") as fh:
            fh.write("p inccnf\n")
            fh.write(body)
            fh.writelines(share)
        p = subprocess.Popen([IGLUCOSE, path, "-verb=0"],
                             stdout=subprocess.PIPE, text=True)
        procs.append((w, path, p, len(share)))
    print(f"{workers} iglucose workers launched", flush=True)

    verdicts = {}
    try:
        while len(verdicts) < workers:
            for w, path, p, n in procs:
                if w in verdicts or p.poll() is None:
                    continue
                out = p.stdout.read()
                sat = "SATISFIABLE" in out.replace("UNSATISFIABLE", "")
                unsat_all = "UNSATISFIABLE" in out
                verdicts[w] = "SAT" if sat else ("UNSAT" if unsat_all else "?")
                print(f"  worker {w}: {verdicts[w]} ({n} cubes)  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
                if verdicts[w] == "SAT":
                    raise KeyboardInterrupt
            time.sleep(5)
    except KeyboardInterrupt:
        for _, _, p, _ in procs:
            p.kill()
        print(f"RESULT: SAT  [{time.time()-t0:.0f}s]", flush=True)
        return
    finally:
        for _, path, _, _ in procs:
            try:
                os.unlink(path)
            except OSError:
                pass

    total = time.time() - t0
    if all(v == "UNSAT" for v in verdicts.values()):
        print(f"RESULT: UNSAT — all workers refuted their stripes  "
              f"[{total:.0f}s]", flush=True)
    else:
        print(f"RESULT: mixed/incomplete {verdicts}  [{total:.0f}s]", flush=True)


if __name__ == "__main__":
    run(sys.argv[1],
        int(sys.argv[2]) if len(sys.argv) > 2 else 6,
        int(sys.argv[3]) if len(sys.argv) > 3 else 12)
