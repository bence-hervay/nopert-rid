"""Negative and consistency controls for the independent checker."""
from pathlib import Path
PROOF = Path(__file__).resolve().parent.parent / "proof"
import json
import random
import sys
from itertools import product
from flint import fmpq
import numpy as np
import model as M
import verify as Vf
from model import Q5, q5

rng = random.Random(99)
V = M.VERTS
report = {}


def rq(lo, hi, den=2**20):
    return lo + (hi - lo) * fmpq(rng.randint(0, den), den)


# ---------------------------------------------------------------- 1. Global: boxes containing r = 0
# At r = 0 the plug is the hole, so every gap of a valid support and every maximum form is >= 0:
# a sound checker must refuse every witness on every box containing an aligned configuration.
acc = tot = 0
for trial in range(120):
    s0 = rq(fmpq(0), fmpq(2, 3) - fmpq(1, 64))
    t0 = rq(fmpq(0), fmpq(2, 5) - fmpq(1, 64))
    w = fmpq(1, 2 ** rng.randint(6, 14))
    box = [(s0, s0 + w), (t0, t0 + w)]
    for k in range(3):
        a = fmpq(rng.randint(0, 4), 2 ** rng.randint(8, 16))
        b = fmpq(rng.randint(0, 4), 2 ** rng.randint(8, 16))
        if rng.random() < 0.3:
            a = fmpq(0)
        box.append((-a, b))
    for _ in range(40):
        e = rng.choice(M.EDGES)
        if rng.random() < 0.5:
            e = (e[1], e[0])
        p = rng.randrange(60)
        if rng.random() < 0.5:
            p = e[0]  # the plug copy of the contact vertex: the most dangerous candidate
        r = Vf.verify_global(box, e[0], e[1], p)
        tot += 1
        acc += r["ok"]
print(f"[1] Global witnesses on boxes containing r=0: {tot} tested, {acc} accepted (must be 0)")
report["global_r0"] = (tot, acc)

# ---------------------------------------------------------------- 2. Global: boxes around touching
# configurations where the plug is a symmetric copy of the hole (r = Cayley vector of a symmetry in B0)
acc = tot = 0
cay = []
for g in M.SYMS:
    if g[0].iszero():
        continue
    r = [c / g[0] for c in g[1]]
    if all(abs(float(c)) <= 0.4 for c in r):
        cay.append(r)
print("   symmetry Cayley vectors inside [-2/5,2/5]^3:", len(cay))
for r in cay:
    for trial in range(20):
        s0 = rq(fmpq(0), fmpq(1, 2))
        t0 = rq(fmpq(0), fmpq(1, 4))
        w = fmpq(1, 2 ** rng.randint(6, 12))
        box = [(s0, s0 + w), (t0, t0 + w)]
        for c in r:
            # rational interval containing the irrational c
            f = float(c)
            lo = fmpq(int(np.floor(f * 2**20)) - rng.randint(1, 5), 2**20)
            hi = fmpq(int(np.ceil(f * 2**20)) + rng.randint(1, 5), 2**20)
            assert q5(lo) < c < q5(hi)
            box.append((lo, hi))
        for _ in range(20):
            e = rng.choice(M.EDGES)
            if rng.random() < 0.5:
                e = (e[1], e[0])
            p = rng.randrange(60)
            tot += 1
            acc += Vf.verify_global(box, e[0], e[1], p)["ok"]
print(f"[2] Global witnesses on boxes containing a symmetry rotation: {tot} tested, {acc} accepted (must be 0)")
report["global_sym"] = (tot, acc)

# ---------------------------------------------------------------- 3. Global: real record boxes with foreign witnesses
lines = open(str(PROOF / "results" / "full" / "search.cert")).read().split("\n")[1:-1]
glob = [json.loads(l[9:]) for l in lines if '"Global"' in l]
sample = rng.sample(glob, 400)
acc_other = acc_rev = acc_parent = 0
for rec in sample:
    box = M.box_of_path(rec["path"])
    a, b = rec["edge"]
    p = rec["vertex"]
    other = rng.choice(glob)
    acc_other += Vf.verify_global(box, other["edge"][0], other["edge"][1], other["vertex"])["ok"]
    acc_rev += Vf.verify_global(box, b, a, p)["ok"]
    parent = M.box_of_path(rec["path"][:-6])
    acc_parent += Vf.verify_global(parent, a, b, p)["ok"]
print(f"[3] 400 Global record boxes: foreign witness accepted {acc_other}, reversed edge accepted {acc_rev}, "
      f"witness on the ancestor 6 levels up accepted {acc_parent}")
report["global_mutations"] = (acc_other, acc_rev, acc_parent)

# ---------------------------------------------------------------- 4. Domain: boxes containing points of D
DP = M.domain_polynomials()
acc = tot = npts = 0
while npts < 150:
    x = [rq(lo, hi) for lo, hi in M.B0]
    x[2] = x[2] / 4
    x[3] = x[3] / 4
    x[4] = x[4] / 4
    if not all(M.peval(p, x).sign() <= 0 for p in DP):
        continue
    npts += 1
    w = fmpq(1, 2 ** rng.randint(4, 30))
    box = []
    for m in range(5):
        lo = x[m] - fmpq(rng.randint(0, 100), 100) * w
        hi = x[m] + fmpq(rng.randint(0, 100), 100) * w
        box.append((lo, hi))
    for k in range(73):
        r = Vf.verify_domain(box, k)
        tot += 1
        acc += r["bernstein"] or r["corner"]
print(f"[4] Domain claims on boxes containing a point of D: {tot} tested, {acc} accepted (must be 0)")
report["domain_inD"] = (tot, acc)
# also degenerate boxes on the boundary: a point with c_0 = 0 exactly is impossible rationally; use fold/affine walls
# r_1 wall: eps r_i + delta(phi-1) r_{i+1} + phi - 2 = 0 has no rational points in B0 (irrational); skip.

# ---------------------------------------------------------------- 5. covers: exact containment vs float sampling
def z_float(cv, s, t, r):
    x = [float(v) for v in cv.view_coords(s, t)]
    rot = cv.rot_affine()
    for j in range(3):
        c, lin = rot[j]
        x.append(float(c) + sum(float(lin[k]) * float(r[k]) for k in range(3)))
    z = []
    for i in range(5):
        z.append(sum(float(cv.Minv[i][j]) * (x[j] - float(cv.c0[j])) for j in range(5)))
    return z


def inside_float(cv, z, tol=1e-12):
    for i in range(5):
        L, H = float(cv.iv[i][0]), float(cv.iv[i][1])
        if z[i] < L - tol:
            return False
        if cv.beyond != i and z[i] > H + tol:
            return False
    return True


names = list(range(30)) + ["square", "pentagon", "arc+", "arc-", "endpoint+", "endpoint-", "crossing+", "crossing-"]
cert_cov = {}
for l in lines:
    if '"Local"' in l or '"Exotic"' in l:
        o = json.loads(l[9:])
        cert_cov.setdefault(o["cover"], []).append(o["path"])
mismatch = 0
yes = no = 0
for nm in names:
    cv = Vf.cover(nm)
    paths = cert_cov.get(nm, [])
    for trial in range(300):
        if paths and trial < 200:
            # perturb a real record box: ancestors, siblings, random moves
            p = rng.choice(paths)
            k = rng.randint(0, min(12, len(p)))
            q = p[: len(p) - k]
            if q and rng.random() < 0.5:
                q = q[:-1] + ("1" if q[-1] == "0" else "0")
            box = M.box_of_path(q)
        else:
            box = M.box_of_path("".join(rng.choice("01") for _ in range(rng.randint(10, 40))))
        ex = cv.contains(box)
        corners = [[box[m][c[m]] for m in range(5)] for c in product((0, 1), repeat=5)]
        fl = all(inside_float(cv, z_float(cv, x[0], x[1], x[2:])) for x in corners)
        # interior random points must agree with a positive answer
        if ex:
            yes += 1
            for _ in range(5):
                x = [rq(box[m][0], box[m][1]) for m in range(5)]
                if not inside_float(cv, z_float(cv, x[0], x[1], x[2:])):
                    fl = False
        else:
            no += 1
        if ex != fl:
            # tolerate exact-boundary ties only when float is within 1e-12
            mismatch += 1
            print("   cover mismatch", nm, q if paths else "", ex, fl)
print(f"[5] cover containment: {yes} contained / {no} not contained test boxes, {mismatch} disagreements with float corner sampling")
report["cover_float"] = (yes, no, mismatch)

# mutated cover records: record box with another cover of the same component
acc = tot = 0
for nm, paths in cert_cov.items():
    for p in rng.sample(paths, min(40, len(paths))):
        box = M.box_of_path(p)
        others = [o for o in names if o != nm and isinstance(o, type(nm))]
        for o in rng.sample(others, min(5, len(others))):
            tot += 1
            acc += Vf.cover(o).contains(box)
print(f"[5b] cover records re-labelled with another cover of the same component: {tot} tested, {acc} accepted "
      "(covers overlap, so some acceptances are legitimate)")
report["cover_relabel"] = (tot, acc)

# ---------------------------------------------------------------- 6. pentagon coordinates and 'beyond'
pc = Vf.cover("pentagon")
phi = M.PHI
# z0 should be s - t/phi - (s_P - t_P/phi), z1 = s + phi t - 1/phi  (event coordinates)
ok = True
sP, tP = pc.c0[0], pc.c0[1]
eta_c = [Q5(1), -(phi.inv())]
xi_c = [Q5(1), phi]
ok &= pc.Minv[0][0] == eta_c[0] and pc.Minv[0][1] == eta_c[1]
ok &= pc.Minv[1][0] == xi_c[0] and pc.Minv[1][1] == xi_c[1]
# at the centre: eta = 0, xi = 0 means sP - tP/phi = 0 and sP + phi tP = 1/phi
ok &= (sP - tP / phi).iszero() and (sP + phi * tP - phi.inv()).iszero()
# D's inequality 0 as a function of z: c_0 = phi s + phi^2 t - 1 = kappa z1
# c_0(c0 + M z) : coefficients
coef_z = [phi * pc.M[0][j] + phi * phi * pc.M[1][j] for j in range(2)]
const = phi * sP + phi * phi * tP - 1
ok &= const.iszero() and coef_z[0].iszero() and coef_z[1].sign() > 0 and pc.beyond == 1 and pc.beyond_ineq == 0
print("[6] pentagon cover coordinates are (eta, xi) = (s - t/phi, s + phi t - 1/phi); c_0 = kappa*xi with kappa =",
      coef_z[1], "> 0 ; beyond axis 1 uses inequality 0:", ok)
report["pentagon_beyond"] = bool(ok)

# covered sets as intervals (for the report)
for nm in names:
    cv = Vf.cover(nm)
    print("   cover", nm, "sigma", cv.sigma, "z-intervals", [(str(a), str(b)) for a, b in cv.iv], "beyond", cv.beyond)

json.dump(report, open("controls-result.json", "w"), indent=1, default=str)
print("CONTROLS DONE")
