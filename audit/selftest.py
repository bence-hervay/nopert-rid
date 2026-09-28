"""Self-tests of the independent model. Every check is exact unless stated."""
import random
from itertools import product
from flint import fmpq
import numpy as np
import model as M
from model import Q5, q5

rng = random.Random(20260928)
ok = True


def check(cond, msg):
    global ok
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        ok = False


# ---------------- the solid
V = M.VERTS
check(len(V) == 60, "60 vertices")
check(all(M.dot(v, v) == Q5(11, 4) for v in V), "all |v|^2 = 11 + 4 sqrt5")
check(M.VKEYS[0] == ((-5, -1), (0, 0), (-3, -1)), "vertex 0 = (-(5+r5)/2, 0, -(3+r5)/2)")
check(all(M.vsub([Q5(0)] * 3, V[i]) == V[59 - i] for i in range(60)) or
      all(any(M.vsub([Q5(0)] * 3, V[i]) == V[j] for j in range(60)) for i in range(60)),
      "central symmetry K = -K")
check([M.vsub([Q5(0)] * 3, V[0])[k] == V[59][k] for k in range(3)] == [True] * 3, "vertex 59 = -vertex 0")
check(len(M.EDGES) == 120, "120 edges (pairs at distance 2)")
deg = [0] * 60
for a, b in M.EDGES:
    deg[a] += 1
    deg[b] += 1
check(all(d == 4 for d in deg), "every vertex has degree 4")
# minimal distance is 2
mind = min(float(M.dot(M.vsub(V[i], V[j]), M.vsub(V[i], V[j]))) for i in range(60) for j in range(i + 1, 60))
check(abs(mind - 4) < 1e-12, "minimal squared distance 4 (float)")

# ---------------- the group
G = M.SYMS
check(len(G) == 60, "60 symmetries")
check(all((g[0] * g[0] + M.dot(g[1], g[1])) == Q5(1) for g in G), "unit quaternions")


def qkey(q):
    return tuple((x.a, x.b) for x in [q[0]] + list(q[1]))


allq = set()
for g in G:
    allq.add(qkey(g))
    allq.add(qkey((-g[0], [-x for x in g[1]])))
check(len(allq) == 120, "+-G has 120 elements")
closed = True
for g in G:
    for h in G:
        if qkey(M.qmul(g, h)) not in allq:
            closed = False
check(closed, "+-G closed under products")
vset = set(tuple((c.a, c.b) for c in v) for v in V)
perms = set()
maps = True
for g in G:
    img = []
    for v in V:
        w = M.qrotate(g, v)
        k = tuple((c.a, c.b) for c in w)
        if k not in vset:
            maps = False
        img.append(k)
    perms.add(tuple(img))
check(maps, "every symmetry maps the vertex set onto itself")
check(len(perms) == 60, "60 distinct permutations")
# 15 elements with scalar part 0 (F2)
check(sum(1 for g in G if g[0].iszero()) == 15, "15 half-turns (F2)")
# F4: the twelve quaternions (phi/2, eps/2 e_i + delta (phi-1)/2 e_{i+1}) lie in +-G
f4 = True
for i in range(3):
    for eps in (-1, 1):
        for delta in (-1, 1):
            v = [Q5(0)] * 3
            v[i] = Q5(fmpq(eps, 2))
            v[(i + 1) % 3] = (M.PHI - 1) * fmpq(delta, 2)
            if qkey((M.PHI * fmpq(1, 2), v)) not in allq:
                f4 = False
check(f4, "F4 quaternions are in the group")

# ---------------- Cayley rotation
def rand_q(den=64):
    return fmpq(rng.randint(-den, den), rng.randint(1, den))


good = True
for trial in range(30):
    r = [rand_q() for _ in range(3)]
    for vi in rng.sample(range(60), 5):
        p = V[vi]
        poly = M.cayley_hat_apply(p)
        val = [M.peval(poly[i], [0, 0, r[0], r[1], r[2]]) for i in range(3)]
        ref = M.qrotate((Q5(1), [q5(x) for x in r]), p)
        if not all(val[i] == ref[i] for i in range(3)):
            good = False
check(good, "Rhat(r)p equals (1,r)(0,p)conj(1,r) at random rational r")
# Rhat Rhat^T = (1+|r|^2)^2 I at random r (via columns)
good = True
for trial in range(20):
    r = [rand_q() for _ in range(3)]
    cols = []
    for k in range(3):
        e = [Q5(0)] * 3
        e[k] = Q5(1)
        poly = M.cayley_hat_apply(e)
        cols.append([M.peval(poly[i], [0, 0, r[0], r[1], r[2]]) for i in range(3)])
    n2 = 1 + r[0] ** 2 + r[1] ** 2 + r[2] ** 2
    for a in range(3):
        for b in range(3):
            if M.dot(cols[a], cols[b]) != (q5(n2 * n2) if a == b else Q5(0)):
                good = False
    # determinant positive: det = (1+|r|^2)^3
    c0, c1, c2 = cols
    det = M.dot(c0, M.cross(c1, c2))
    if det != q5(n2 ** 3):
        good = False
check(good, "Rhat orthogonal with factor (1+|r|^2)^2, det (1+|r|^2)^3")

# ---------------- Bernstein
def bern_value(bA, bB, degs, y):
    """sum_k b_k prod B_{k_m,n_m}(y_m) for one polynomial (exact)."""
    from math import comb
    tot = Q5(0)
    for idx in product(*[range(d + 1) for d in degs]):
        w = fmpq(1)
        for m in range(5):
            k, n = idx[m], degs[m]
            w *= comb(n, k) * y[m] ** k * (1 - y[m]) ** (n - k)
        tot = tot + Q5(bA[idx], bB[idx]) * w
    return tot


good_c = good_i = True
for trial in range(25):
    degs = [rng.randint(0, 2) for _ in range(5)]
    poly = {}
    for e in product(*[range(d + 1) for d in degs]):
        if rng.random() < 0.6:
            poly[e] = Q5(rand_q(9), rand_q(9))
    box = []
    for m in range(5):
        lo = rand_q(32)
        box.append((lo, lo + fmpq(rng.randint(0, 20), rng.randint(1, 64))))
    A, B, K = M.poly_to_int_tensors([poly], degs)
    bA, bB, sc = M.bernstein_int(A, B, degs, box)
    bA, bB = bA[0], bB[0]
    fac = fmpq(1, sc * K)
    # corners
    for c in product((0, 1), repeat=5):
        idx = tuple(degs[m] * c[m] for m in range(5))
        x = [box[m][c[m]] for m in range(5)]
        if Q5(bA[idx] * fac, bB[idx] * fac) != M.peval(poly, x):
            good_c = False
    # interior points
    for _ in range(3):
        y = [fmpq(rng.randint(0, 50), 50) for _ in range(5)]
        x = [box[m][0] + (box[m][1] - box[m][0]) * y[m] for m in range(5)]
        bv = bern_value(bA, bB, degs, y)
        if Q5(bv.a * fac, bv.b * fac) != M.peval(poly, x):
            good_i = False
check(good_c, "Bernstein corner coefficients = values at corners (random polys/boxes)")
check(good_i, "Bernstein expansion reproduces values at random interior points")

# sign test
good = True
for _ in range(20000):
    a = rng.randint(-10**6, 10**6)
    b = rng.randint(-10**6, 10**6)
    s = M.sign_q5(a, b)
    f = a + b * 5 ** 0.5
    if abs(f) > 1e-3 and (s > 0) != (f > 0):
        good = False
# near-cancellation: Pell-like pairs
for (a, b) in [(9, -4), (-9, 4), (161, -72), (-161, 72), (2889, -1292), (-2889, 1292), (0, 0), (5, 0), (0, -3)]:
    f = a + b * 5 ** 0.5
    s = M.sign_q5(a, b)
    if (a, b) == (0, 0):
        good = good and s == 0
    else:
        good = good and s == (1 if f > 0 else -1)
check(good, "exact sign of a + b sqrt5 (random and near-cancelling)")

# ---------------- documented examples
# Figure 3: box, witness edge 6->13, plug 13; coefficients in [-0.113, -0.067], edge valid at corner views
def gap_poly(a, b, p, w=None):
    d = M.vsub(V[b], V[a])
    h = V[a] if w is None else V[w]
    n = M.pcross(M.U, [M.pconst(c) for c in d])
    nh = M.pdot(n, [M.pconst(c) for c in h])
    Rp = M.cayley_hat_apply(V[p])
    return M.psub(M.pmul(M.padd(M.ONE, M.RR), nh), M.pdot(n, Rp))


def valid_at(a, b, s, t):
    d = M.vsub(V[b], V[a])
    n = M.cross([q5(s), q5(t), Q5(1)], d)
    return all(M.dot(n, M.vsub(V[a], w)).sign() >= 0 for w in V)


box3 = [(fmpq(3, 16), fmpq(7, 32)), (fmpq(3, 32), fmpq(7, 64)), (fmpq(1, 8), fmpq(9, 64)),
        (fmpq(1, 16), fmpq(5, 64)), (fmpq(-1, 32), fmpq(-1, 64))]
g = gap_poly(6, 13, 13)
check(len(g) == 26, f"Figure 3 gap has 26 terms (got {len(g)})")
check(M.poly_degrees(g) == [1, 1, 2, 2, 2], "Figure 3 gap degrees (1,1,2,2,2)")
A, B, K = M.poly_to_int_tensors([g], [1, 1, 2, 2, 2])
bA, bB, sc = M.bernstein_int(A, B, [1, 1, 2, 2, 2], box3)
vals = [float(Q5(fmpq(int(x), sc * K), fmpq(int(y), sc * K))) for x, y in zip(bA.ravel(), bB.ravel())]
print("   Figure 3 coefficient range:", round(min(vals), 4), round(max(vals), 4), len(vals))
check(abs(min(vals) + 0.113) < 0.0006 and abs(max(vals) + 0.067) < 0.0006, "Figure 3 range -0.113 .. -0.067")
check(all(valid_at(6, 13, s, t) for s in box3[0] for t in box3[1]), "Figure 3 edge valid at the four corner views")
bigbox = [(fmpq(0), fmpq(1, 3)), (fmpq(0), fmpq(1, 5)), (fmpq(1, 16), fmpq(3, 16)), (fmpq(0), fmpq(1, 8)), (fmpq(-1, 16), fmpq(0))]
bA, bB, sc = M.bernstein_int(A, B, [1, 1, 2, 2, 2], bigbox)
vals = [float(Q5(fmpq(int(x), sc * K), fmpq(int(y), sc * K))) for x, y in zip(bA.ravel(), bB.ravel())]
print("   Figure 3 larger box range:", round(min(vals), 3), round(max(vals), 3))
check(abs(min(vals) + 0.34) < 0.006 and abs(max(vals) - 0.32) < 0.006, "larger box range -0.34 .. 0.32")
bad = []
for s in bigbox[0]:
    for t in bigbox[1]:
        d = M.vsub(V[13], V[6])
        n = M.cross([q5(s), q5(t), Q5(1)], d)
        for wi, w in enumerate(V):
            if M.dot(n, M.vsub(V[6], w)).sign() < 0:
                bad.append(wi)
check(17 in bad, f"vertex 17 beyond the line at a corner view of the larger box (violators {sorted(set(bad))})")

# Domain: second arc theta=-e lies outside D: c_13 (sign +) and c_24 (sign -) equal 5e^2+5e^4 in homogeneous form
a1 = Q5(fmpq(-1, 2), fmpq(3, 10))
a3 = Q5(fmpq(-1, 2), fmpq(1, 10))
DP = M.domain_polynomials()
good = True
for sigma, idx in ((1, 13), (-1, 24)):
    for e in [fmpq(1, 7), fmpq(1, 5), fmpq(3, 100), fmpq(1, 9)]:
        th = -e
        s = (1 - 2 * e) / (2 + e)
        r = [sigma * a1 + th * sigma * a3, q5(th), sigma * a3 - th * sigma * a1]
        val = M.peval(DP[idx], [q5(s), Q5(0)] + r)
        # homogeneous degree 2 in u: value at u=(1-2e,0,2+e) = (2+e)^2 * affine value
        hom = val * (2 + e) ** 2
        if hom != q5(5 * e ** 2 + 5 * e ** 4):
            good = False
            print("   ", sigma, e, hom)
check(good, "folds 13 and 24 equal 5e^2+5e^4 on the second arcs (documented identity; checks group order)")
# the arc planes are no-fit sets: max_w e2.R(sigma r_A) w = 2 + sqrt5 = max_w e2.w
good = True
for sigma in (1, -1):
    rA = [sigma * a1, Q5(0), sigma * a3]
    n2 = 1 + M.dot(rA, rA)
    mx = None
    for w in V:
        y = M.qrotate((Q5(1), rA), w)[1] / n2
        mx = y if mx is None or y > mx else mx
    good = good and mx == Q5(2, 1) and max((w[1] for w in V), key=lambda x: float(x)) == Q5(2, 1)
check(good, "arc planes: plug reach in e2 equals hole reach 2+sqrt5 (both signs)")
# D contains the aligned configuration (1/5, 1/10, 0,0,0)
x = [q5(fmpq(1, 5)), q5(fmpq(1, 10)), Q5(0), Q5(0), Q5(0)]
check(all(M.peval(p, x).sign() <= 0 for p in DP), "aligned configuration (1/5,1/10,0) is in D")
# the three examples of documented Domain inequality index for r = (-t, s, 0) near square: fold 29 limit
s0, t0 = fmpq(1, 100), fmpq(3, 500)
x = [q5(s0), q5(t0), q5(-t0), q5(s0), Q5(0)]
check(M.peval(DP[29], x).sign() > 0, "fold 29 positive at r=(-t,s,0) near the square view (documented)")

print("SELFTEST", "OK" if ok else "FAILED")
