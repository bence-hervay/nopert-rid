"""Verify every record of the certificate in a uniformly random order (seed 20260928),
12 worker processes. The first 20,000 records of that order form the uniform random
sample; the run continues through all records. Results go to results.jsonl."""
from pathlib import Path
PROOF = Path(__file__).resolve().parent.parent / "proof"
import json
import os
import random
import sys
import time
from multiprocessing import Pool

CERT = str(PROOF / "results" / "full" / "search.cert")


def init():
    global verify_record
    from verify import verify_record as vr
    verify_record = vr


def work(item):
    idx, txt = item
    try:
        r = verify_record(txt)
    except Exception as e:  # any exception is a refusal, reported
        r = {"ok": False, "error": repr(e), "path": json.loads(txt)["path"], "component": json.loads(txt)["component"]}
    r["line"] = idx
    return r


if __name__ == "__main__":
    import verify  # noqa: F401  (a missing dependency fails here, not in every worker)
    lines = open(CERT).read().split("\n")[1:-1]
    items = [(i + 2, l[9:]) for i, l in enumerate(lines)]  # file line numbers (header = line 1)
    order = list(range(len(items)))
    random.Random(20260928).shuffle(order)
    todo = [items[k] for k in order]
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else len(todo)
    todo = todo[:limit]
    out = open("results.jsonl", "w")
    t0 = time.time()
    n = bad = 0
    with Pool(12, initializer=init) as pool:
        for r in pool.imap(work, todo, chunksize=64):
            out.write(json.dumps(r, separators=(",", ":")) + "\n")
            n += 1
            if not r["ok"]:
                bad += 1
                print("REFUSED", r, flush=True)
            if n % 10000 == 0:
                out.flush()
                print(f"{n} done, {bad} refused, {time.time() - t0:.0f} s", flush=True)
    out.close()
    print(f"FINISHED {n} records, {bad} refused, {time.time() - t0:.0f} s", flush=True)
