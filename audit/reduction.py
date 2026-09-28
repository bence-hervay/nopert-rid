"""Floating-point check (evidence, not proof) that the domain D, as built here from the
documented inequalities, contains a representative of random configurations: executes the
documented reduction (step 3: the view maximising c.u over the orbit, c = (1,2,9); step 4: the
quaternion of largest scalar part among +-q b, +-j q b) and checks the representative lies in
D and B0 and has an isometric pair of shadows (same signed outside-margin)."""
import numpy as np
from scipy.spatial import ConvexHull
import model as M

rng = np.random.default_rng(11)
V = np.array([[float(c) for c in v] for v in M.VERTS])
G = np.array([[float(g[0])] + [float(c) for c in g[1]] for g in M.SYMS])
GT = np.concatenate([G, -G])  # 120 unit quaternions
DPF = [[(e, float(c)) for e, c in p.items()] for p in M.domain_polynomials()]


def qm(p, q):
    w1, v1 = p[0], p[1:]
    w2, v2 = q[0], q[1:]
    return np.concatenate([[w1 * w2 - v1 @ v2], w1 * v2 + w2 * v1 + np.cross(v1, v2)])


def qmat(q):
    w, x, y, z = q
    return np.array([[w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z]]) / (q @ q)


def margin(u, R):
    """signed max distance of plug shadow outside hole shadow, orthogonal projection on u-perp."""
    u = u / np.linalg.norm(u)
    a = np.cross(u, [1.0, 0, 0] if abs(u[0]) < 0.9 else [0, 1.0, 0])
    a /= np.linalg.norm(a)
    b = np.cross(u, a)
    H = np.stack([V @ a, V @ b], 1)
    P = V @ R.T
    P = np.stack([P @ a, P @ b], 1)
    hull = ConvexHull(H)
    return (P @ hull.equations[:, :2].T + hull.equations[:, 2]).max()


def inD(x):
    worst = -1e9
    for p in DPF:
        v = 0.0
        for e, c in p:
            v += c * x[0] ** e[0] * x[1] ** e[1] * x[2] ** e[2] * x[3] ** e[3] * x[4] ** e[4]
        worst = max(worst, v)
    return worst


c = np.array([1.0, 2.0, 9.0])
fails = 0
mism = 0
N = 3000
for it in range(N):
    e = rng.normal(size=3)
    e /= np.linalg.norm(e)
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    m0 = margin(e, qmat(q))
    # step 3: u* = sigma g e maximising c.u*
    best = None
    for g in G:
        ge = qmat(g) @ e
        for sg in (1, -1):
            val = c @ (sg * ge)
            if best is None or val > best[0]:
                best = (val, sg * ge, g)
    _, us, g = best
    # relative rotation g R g^-1
    gi = g * np.array([1, -1, -1, -1])
    q1 = qm(qm(g, q), gi)
    u = us / us[2]
    j = np.concatenate([[0.0], u / np.linalg.norm(u)])
    S = []
    for b in GT:
        qb = qm(q1, b)
        jqb = qm(j, qb)
        S += [qb, -qb, jqb, -jqb]
    S = np.array(S)
    qs = S[np.argmax(S[:, 0])]
    r = qs[1:] / qs[0]
    x = [u[0], u[1], r[0], r[1], r[2]]
    inb0 = (0 <= x[0] <= 2 / 3 + 1e-12) and (0 <= x[1] <= 0.4 + 1e-12) and all(abs(v) <= 0.4 + 1e-12 for v in r)
    d = inD(x)
    if not (inb0 and d <= 1e-9):
        fails += 1
        print("NOT IN D", x, d)
    m1 = margin(np.array([u[0], u[1], 1.0]), qmat(np.concatenate([[1.0], r])))
    if abs(m1 - m0) > 1e-9:
        mism += 1
        print("MARGIN MISMATCH", m0, m1)
print(f"reduction: {N} random configurations, {fails} representatives outside D/B0, {mism} shadow-margin mismatches")
