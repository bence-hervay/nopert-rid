"""Floating-point semantic oracle (evidence only, not proof): at sample points of record
boxes, is the configuration a fit (plug shadow strictly inside hole shadow), and is it in D?
A sound certificate never has a sampled point that is a fit AND in D. For Global records the
sampled points must not even be weakly contained (plug vertex p beyond the hole shadow)."""
from pathlib import Path
PROOF = Path(__file__).resolve().parent.parent / "proof"
import json
import random
import numpy as np
from scipy.spatial import ConvexHull
import model as M

rng = random.Random(5)
V = np.array([[float(c) for c in v] for v in M.VERTS])
DPf = None


def rot(r):
    r = np.asarray(r, float)
    n2 = r @ r
    Rh = (1 - n2) * np.eye(3) + 2 * np.outer(r, r) + 2 * np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
    return Rh / (1 + n2)


def shadow(X, s, t):
    return np.stack([X[:, 0] - s * X[:, 2], X[:, 1] - t * X[:, 2]], axis=1)


def outside_margin(s, t, r):
    """max over plug vertices of the signed distance beyond the hole-shadow hull (>0: poke)."""
    H = shadow(V, s, t)
    hull = ConvexHull(H)
    P = shadow(V @ rot(r).T, s, t)
    eq = hull.equations  # a x + b y + c <= 0 inside, (a,b) unit
    return (P @ eq[:, :2].T + eq[:, 2]).max(axis=1)  # per plug vertex


# D in float: evaluate the exact polynomials at float points via their Q5 coefficients
DP = M.domain_polynomials()
DPF = [[(e, float(c)) for e, c in p.items()] for p in DP]


def in_D(x):
    vals = []
    for p in DPF:
        v = 0.0
        for e, c in p:
            v += c * x[0] ** e[0] * x[1] ** e[1] * x[2] ** e[2] * x[3] ** e[3] * x[4] ** e[4]
        vals.append(v)
    return max(vals)


lines = open(str(PROOF / "results" / "full" / "search.cert")).read().split("\n")[1:-1]
recs = [json.loads(l[9:]) for l in lines]
by = {}
for o in recs:
    by.setdefault(o["component"], []).append(o)
plan = {"Global": 4000, "Domain": 2000, "Local": 2000, "Exotic": len(by["Exotic"])}
summary = {}
for comp, k in plan.items():
    fits_in_D = 0
    weak_global = 0
    npts = 0
    inD = 0
    for o in rng.sample(by[comp], k):
        box = M.box_of_path(o["path"])
        for j in range(8):
            if j == 0:
                x = [float((lo + hi) / 2) for lo, hi in box]
            else:
                x = [float(lo) + (float(hi) - float(lo)) * rng.random() for lo, hi in box]
            npts += 1
            m = outside_margin(x[0], x[1], x[2:])
            d = in_D(x)
            if d <= 1e-12:
                inD += 1
            if m.max() < -1e-12 and d <= 1e-12:
                fits_in_D += 1
                print("FIT IN D?", comp, o, x)
            if comp == "Global":
                p = o["vertex"]
                if m[p] <= 1e-12:
                    # plug vertex p not beyond the hole shadow: witness claim would be false here
                    weak_global += 1
                    print("GLOBAL VERTEX NOT OUTSIDE", o, x, m[p])
            if comp == "Domain":
                if not d > -1e-12:
                    print("DOMAIN POINT INSIDE D?", o, x, d)
    summary[comp] = {"records": k, "points": npts, "points_in_D(float)": inD, "fits_in_D": fits_in_D,
                     "global_vertex_not_outside": weak_global}
    print(comp, summary[comp], flush=True)
json.dump(summary, open("oracle-result.json", "w"), indent=1)
print("ORACLE DONE")
