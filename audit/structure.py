"""Independent structural check of the certificate: line framing, SHA-256
prefixes, canonical JSON spelling, record types and ranges, path order,
prefix-freeness and complete coverage of B0 (three independent arguments)."""
from pathlib import Path
PROOF = Path(__file__).resolve().parent.parent / "proof"
import hashlib
import json
import re
import sys
from collections import Counter
from flint import fmpq
import model as M

CERT = str(PROOF / "results" / "full" / "search.cert")
HEADER = re.compile(r'\{"format":"rid-certificate/1","policy":"[0-9a-f]{64}","root":"","max_depth":4096\}')
EXC = ["square", "pentagon", "arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"]
problems = []


def bad(i, msg):
    problems.append((i, msg))
    if len(problems) < 30:
        print("PROBLEM line", i, msg)


raw = open(CERT, "rb").read()
print("bytes", len(raw))
check_ascii = all(b < 128 for b in raw)
print("pure ASCII:", check_ascii)
if not raw.endswith(b"\n"):
    bad(-1, "file does not end with newline")
lines = raw.split(b"\n")[:-1]
print("lines", len(lines))
LINE = re.compile(rb"^[0-9a-f]{8} (\{.*\})$")
paths = []
comps = Counter()
records = []
for i, ln in enumerate(lines):
    m = LINE.match(ln)
    if not m:
        bad(i, "bad frame")
        continue
    payload = m.group(1)
    h = hashlib.sha256(payload).hexdigest()[:8].encode()
    if h != ln[:8]:
        bad(i, "hash mismatch")
    txt = payload.decode()
    if i == 0:
        if not HEADER.fullmatch(txt):
            bad(i, "header differs")
        continue
    # parse with duplicate-key detection
    def hook(pairs):
        keys = [k for k, _ in pairs]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate key")
        return dict(pairs)
    try:
        obj = json.loads(txt, object_pairs_hook=hook)
    except Exception as e:
        bad(i, f"json {e}")
        continue
    # canonical re-serialisation (compact, key order preserved)
    if json.dumps(obj, separators=(",", ":")) != txt:
        bad(i, "not canonical compact JSON")
    keys = list(obj)
    comp = obj.get("component")
    p = obj.get("path")
    if not isinstance(p, str) or not re.fullmatch(r"[01]*", p) or len(p) > 4096:
        bad(i, "bad path")
        continue
    if comp == "Domain":
        ok = keys == ["path", "component", "inequality"] and type(obj["inequality"]) is int and 0 <= obj["inequality"] < 73
    elif comp == "Global":
        w = obj
        ok = (keys == ["path", "component", "edge", "vertex"]
              and isinstance(w["edge"], list) and len(w["edge"]) == 2 and all(type(x) is int for x in w["edge"])
              and type(w["vertex"]) is int and 0 <= w["vertex"] < 60
              and tuple(sorted(w["edge"])) in M.EDGESET)
    elif comp == "Local":
        ok = keys == ["path", "component", "cover"] and type(obj["cover"]) is int and 0 <= obj["cover"] < 30
    elif comp == "Exotic":
        ok = keys == ["path", "component", "cover"] and obj["cover"] in EXC
    else:
        ok = False
    if not ok:
        bad(i, f"bad record shape {txt}")
    comps[comp] += 1
    paths.append(p)
    records.append((i, txt))

print("records", len(paths), dict(comps))

# order: strictly increasing (length, lexical)
order_ok = all((len(paths[k]), paths[k]) < (len(paths[k + 1]), paths[k + 1]) for k in range(len(paths) - 1))
print("strict search order:", order_ok)
if not order_ok:
    bad(-1, "order")

# prefix-free (no path equals or is a proper prefix of another)
pset = set(paths)
dup = len(pset) != len(paths)
pref = [p for p in paths if any(p[:k] in pset for k in range(len(p)))]
print("duplicates:", dup, " paths having a recorded proper prefix:", len(pref))
if dup or pref:
    bad(-1, "not prefix-free")

# coverage 1: Kraft sum exactly 1
Lmax = max(len(p) for p in paths)
kraft = sum(1 << (Lmax - len(p)) for p in paths)
print("max depth", Lmax, " Kraft sum == 1:", kraft == (1 << Lmax))

# coverage 2: full binary tree - every internal node (proper prefix) has both children present
internal = set()
for p in paths:
    for k in range(len(p)):
        internal.add(p[:k])
full = all(((q + "0") in internal or (q + "0") in pset) and ((q + "1") in internal or (q + "1") in pset) for q in internal)
root_ok = "" in internal or "" in pset
disjoint = not (internal & pset)
print("internal nodes", len(internal), " every internal node has both children:", full, " root present:", root_ok,
      " leaves = records only:", disjoint)

# coverage 3: exact volume of the boxes built independently from the paths
vol0 = fmpq(1)
for lo, hi in M.B0:
    vol0 *= hi - lo
tot = fmpq(0)
for p in paths:
    b = M.box_of_path(p)
    v = fmpq(1)
    for lo, hi in b:
        v *= hi - lo
    tot += v
print("sum of box volumes == vol(B0):", tot == vol0)

# coverage 4: point probes - random exact points of B0 (and points on B0's faces/corners) lie in a record box
import random
rng = random.Random(7)
by_len = {}
for p in paths:
    by_len.setdefault(len(p), set()).add(p)


def locate(x):
    """All records whose closed box contains x, found by walking both halves at ties."""
    found = []
    stack = [("", [list(iv) for iv in M.B0])]
    while stack:
        q, box = stack.pop()
        if q in pset:
            found.append(q)
            continue
        if q not in internal:
            return None  # hole in the tree
        ax = len(q) % 5
        lo, hi = box[ax]
        mid = (lo + hi) / 2
        if x[ax] <= mid:
            b = [list(iv) for iv in box]
            b[ax] = [lo, mid]
            stack.append((q + "0", b))
        if x[ax] >= mid:
            b = [list(iv) for iv in box]
            b[ax] = [mid, hi]
            stack.append((q + "1", b))
    return found


probe_ok = True
nprobe = 0
for _ in range(20000):
    x = []
    for lo, hi in M.B0:
        c = rng.random()
        if c < 0.05:
            x.append(lo)
        elif c < 0.10:
            x.append(hi)
        elif c < 0.2:
            k = rng.randint(0, 64)
            x.append(lo + (hi - lo) * fmpq(k, 64))  # on bisection hyperplanes
        else:
            x.append(lo + (hi - lo) * fmpq(rng.randint(0, 10**9), 10**9))
    f = locate(x)
    nprobe += 1
    if not f:
        probe_ok = False
        bad(-1, f"point not covered {x}")
    else:
        for q in f:
            b = M.box_of_path(q)
            if not all(b[m][0] <= x[m] <= b[m][1] for m in range(5)):
                probe_ok = False
print(f"point probes: {nprobe} exact points of B0 (incl. faces and bisection planes) all covered:", probe_ok)

print("STRUCTURE PROBLEMS:", len(problems))
json.dump({"records": len(paths), "components": dict(comps), "problems": problems[:100], "max_depth": Lmax},
          open("structure-result.json", "w"), indent=1)
