"""Independent exact model of the RID Rupert problem (written from the documented
mathematical definitions only; no code of the crate is imported or copied).

Numbers of Q(sqrt5) are Q5(a, b) = a + b*sqrt5 with a, b flint.fmpq.
Polynomials are sparse dicts {exponent 5-tuple over (s,t,r1,r2,r3): Q5}.
"""
from itertools import product
from math import comb
from flint import fmpq, fmpz
import numpy as np


# ---------------------------------------------------------------- Q(sqrt5)
class Q5:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, fmpq) else fmpq(a)
        self.b = b if isinstance(b, fmpq) else fmpq(b)

    def __add__(self, o):
        o = q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __sub__(self, o):
        o = q5(o)
        return Q5(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return q5(o) - self

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __mul__(self, o):
        o = q5(o)
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def conj(self):
        return Q5(self.a, -self.b)

    def norm(self):
        return self.a * self.a - 5 * self.b * self.b

    def inv(self):
        n = self.norm()
        assert n != 0
        c = self.conj()
        return Q5(c.a / n, c.b / n)

    def __truediv__(self, o):
        return self * q5(o).inv()

    def __eq__(self, o):
        o = q5(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((int(self.a.p), int(self.a.q), int(self.b.p), int(self.b.q)))

    def iszero(self):
        return self.a == 0 and self.b == 0

    def sign(self):
        return sign_q5(self.a, self.b)

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __gt__(self, o):
        return (self - o).sign() > 0

    def __ge__(self, o):
        return (self - o).sign() >= 0

    def __float__(self):
        return float(self.a) + float(self.b) * 5 ** 0.5

    def __repr__(self):
        return f"({self.a} + {self.b}*r5)"


def q5(x):
    if isinstance(x, Q5):
        return x
    return Q5(x, 0)


def sign_q5(a, b):
    """Exact sign of a + b*sqrt5 (a, b rational or integer)."""
    if a >= 0 and b >= 0:
        return 0 if (a == 0 and b == 0) else 1
    if a <= 0 and b <= 0:
        return -1
    # opposite strict signs
    if a > 0:  # b < 0 : compare a with |b| sqrt5
        return 1 if a * a > 5 * b * b else -1
    # a < 0, b > 0
    return 1 if 5 * b * b > a * a else -1


R5 = Q5(0, 1)
PHI = Q5(fmpq(1, 2), fmpq(1, 2))


def half(A, B):
    """(A + B sqrt5)/2 for integers A, B."""
    return Q5(fmpq(A, 2), fmpq(B, 2))


# ---------------------------------------------------------------- vectors
def vadd(x, y):
    return [x[i] + y[i] for i in range(3)]


def vsub(x, y):
    return [x[i] - y[i] for i in range(3)]


def vscale(c, x):
    return [c * x[i] for i in range(3)]


def dot(x, y):
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]


def cross(x, y):
    return [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]]


# ---------------------------------------------------------------- the RID
def rid_vertices():
    """The 60 vertices: all sign changes and cyclic permutations of the three
    generating points (edge length 2), duplicates removed, sorted by the
    documented key: each coordinate (A + B sqrt5)/2, lexicographic order of
    the triple of integer pairs (A, B)."""
    gens = [
        (Q5(1), Q5(1), Q5(2, 1)),  # (1, 1, 2+sqrt5)
        (half(3, 1), half(1, 1), Q5(1, 1)),  # ((3+r5)/2, (1+r5)/2, 1+r5)
        (half(5, 1), Q5(0), half(3, 1)),  # ((5+r5)/2, 0, (3+r5)/2)
    ]
    pts = set()
    for g in gens:
        for signs in product((1, -1), repeat=3):
            p = tuple(g[i] * signs[i] for i in range(3))
            for k in range(3):
                q = (p[k % 3], p[(k + 1) % 3], p[(k + 2) % 3])
                pts.add(q)

    def key(p):
        out = []
        for c in p:
            A, B = c.a * 2, c.b * 2
            assert A.q == 1 and B.q == 1
            out.append((int(A.p), int(B.p)))
        return tuple(out)

    vs = sorted(pts, key=key)
    return [list(v) for v in vs], [key(v) for v in vs]


VERTS, VKEYS = rid_vertices()


def rid_edges():
    E = []
    for i in range(60):
        for j in range(i + 1, 60):
            d = vsub(VERTS[i], VERTS[j])
            if dot(d, d) == Q5(4):
                E.append((i, j))
    return E


EDGES = rid_edges()
EDGESET = set(EDGES)


# ---------------------------------------------------------------- quaternions
def qmul(p, q):
    """Hamilton product of quaternions given as (w, [x,y,z])."""
    pw, pv = p
    qw, qv = q
    w = pw * qw - dot(pv, qv)
    v = vadd(vadd(vscale(pw, qv), vscale(qw, pv)), cross(pv, qv))
    return (w, v)


def qconj(q):
    return (q[0], [-q[1][0], -q[1][1], -q[1][2]])


def qrotate(q, p):
    """q (0,p) conj(q)  = |q|^2 R(q) p."""
    return qmul(qmul(q, (Q5(0), list(p))), qconj(q))[1]


def symmetries():
    """The 60 rotation quaternions (one of each pair +-g, first nonzero coordinate
    positive), sorted by the documented key: each coordinate (a + b sqrt5)/4,
    lexicographic order of the integer pairs."""
    q4 = lambda a, b: Q5(fmpq(a, 4), fmpq(b, 4))
    cands = []
    # 8 coordinate units
    for k in range(4):
        for sg in (1, -1):
            c = [Q5(0)] * 4
            c[k] = Q5(sg)
            cands.append(c)
    # 16 (+-1,+-1,+-1,+-1)/2
    for sg in product((1, -1), repeat=4):
        cands.append([Q5(fmpq(x, 2)) for x in sg])
    # 96 even permutations of (0, +-1/2, +-phi/2, +-(phi-1)/2)
    base = [Q5(0), q4(2, 0), q4(1, 1), q4(-1, 1)]

    def parity(perm):
        perm = list(perm)
        s = 0
        for i in range(len(perm)):
            for j in range(i + 1, len(perm)):
                if perm[i] > perm[j]:
                    s += 1
        return s % 2

    from itertools import permutations

    for perm in permutations(range(4)):
        if parity(perm):
            continue
        for sg in product((1, -1), repeat=3):
            vals = [base[0], base[1] * sg[0], base[2] * sg[1], base[3] * sg[2]]
            c = [None] * 4
            # position perm[i] receives vals[i]
            for i in range(4):
                c[perm[i]] = vals[i]
            cands.append(c)
    assert len(cands) == 120
    reps = {}
    for c in cands:
        first = next(x for x in c if not x.iszero())
        if first.sign() < 0:
            c = [-x for x in c]
        key = []
        for x in c:
            A, B = x.a * 4, x.b * 4
            assert A.q == 1 and B.q == 1
            key.append((int(A.p), int(B.p)))
        reps[tuple(key)] = c
    assert len(reps) == 60, len(reps)
    keys = sorted(reps)
    return [(reps[k][0], reps[k][1:]) for k in keys]


SYMS = symmetries()


# ---------------------------------------------------------------- sparse polynomials
NV = 5  # s, t, r1, r2, r3


def pconst(c):
    c = q5(c)
    return {} if c.iszero() else {(0,) * NV: c}


def pvar(i):
    e = [0] * NV
    e[i] = 1
    return {tuple(e): Q5(1)}


def padd(p, q):
    out = dict(p)
    for e, c in q.items():
        if e in out:
            v = out[e] + c
            if v.iszero():
                del out[e]
            else:
                out[e] = v
        else:
            out[e] = c
    return out


def pneg(p):
    return {e: -c for e, c in p.items()}


def psub(p, q):
    return padd(p, pneg(q))


def pmul(p, q):
    out = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(e1[i] + e2[i] for i in range(NV))
            v = c1 * c2
            if e in out:
                v = out[e] + v
                if v.iszero():
                    del out[e]
                else:
                    out[e] = v
            elif not v.iszero():
                out[e] = v
    return out


def pscale(c, p):
    c = q5(c)
    if c.iszero():
        return {}
    return {e: c * v for e, v in p.items()}


def peval(p, x):
    """Exact evaluation at a point x of 5 Q5 (or fmpq) values."""
    x = [q5(v) for v in x]
    tot = Q5(0)
    for e, c in p.items():
        term = c
        for i in range(NV):
            for _ in range(e[i]):
                term = term * x[i]
        tot = tot + term
    return tot


S_, T_, R1_, R2_, R3_ = (pvar(i) for i in range(5))
ONE = pconst(1)
U = [S_, T_, ONE]  # affine view (s, t, 1)
RV = [R1_, R2_, R3_]


def pdot(x, y):
    return padd(padd(pmul(x[0], y[0]), pmul(x[1], y[1])), pmul(x[2], y[2]))


def pcross(x, y):
    return [
        psub(pmul(x[1], y[2]), pmul(x[2], y[1])),
        psub(pmul(x[2], y[0]), pmul(x[0], y[2])),
        psub(pmul(x[0], y[1]), pmul(x[1], y[0])),
    ]


RR = pdot(RV, RV)  # |r|^2


def cayley_hat_apply(p):
    """Rhat(r) p = (1-|r|^2) p + 2 (r.p) r + 2 r x p, as three polynomials, for a
    constant vector p (the rotation of the quaternion (1, r) times 1+|r|^2)."""
    P = [pconst(c) for c in p]
    one_minus = psub(ONE, RR)
    rp = pdot(RV, P)
    rxp = pcross(RV, P)
    return [padd(padd(pmul(one_minus, P[i]), pscale(2, pmul(rp, RV[i]))), pscale(2, rxp[i])) for i in range(3)]


# ---------------------------------------------------------------- the domain D
def domain_polynomials():
    """The 73 polynomials c_i (D = {x in B0 : c_i(x) <= 0}) in (s,t,r1,r2,r3),
    with the affine view u = (s, t, 1)."""
    polys = []
    # 0: phi u1 + phi^2 u2 - u3
    polys.append(padd(padd(pscale(PHI, S_), pscale(PHI * PHI, T_)), pconst(-1)))
    # 1-12: eps r_i + delta (phi-1) r_{i+1} + phi - 2 ; order i, eps, delta; -1 before +1
    for i in range(3):
        for eps in (-1, 1):
            for delta in (-1, 1):
                p = padd(pscale(eps, RV[i]), pscale((PHI - 1) * delta, RV[(i + 1) % 3]))
                p = padd(p, pconst(PHI - 2))
                polys.append(p)
    # 13-72: L_g^2 - |u|^2, L_g = scalar((0,u)(1,r) g)
    for g in SYMS:
        L = poly_quat_scalar_0u_1r_g(g)
        polys.append(psub(pmul(L, L), pdot(U, U)))
    assert len(polys) == 73
    return polys


def poly_quat_scalar_0u_1r_g(g):
    """scalar part of (0,u)(1,r)g computed by polynomial Hamilton products."""
    # quaternion with polynomial entries
    def pqmul(p, q):
        pw, pv = p
        qw, qv = q
        w = psub(pmul(pw, qw), pdot(pv, qv))
        v = [padd(padd(pmul(pw, qv[i]), pmul(qw, pv[i])), pcross(pv, qv)[i]) for i in range(3)]
        return (w, v)

    a = ({}, U)
    b = (ONE, RV)
    gq = (pconst(g[0]), [pconst(c) for c in g[1]])
    return pqmul(pqmul(a, b), gq)[0]


# ---------------------------------------------------------------- boxes
B0 = [(fmpq(0), fmpq(2, 3)), (fmpq(0), fmpq(2, 5))] + [(fmpq(-2, 5), fmpq(2, 5))] * 3


def box_of_path(path):
    box = [list(iv) for iv in B0]
    for d, ch in enumerate(path):
        ax = d % 5
        lo, hi = box[ax]
        mid = (lo + hi) / 2
        if ch == "0":
            box[ax] = [lo, mid]
        elif ch == "1":
            box[ax] = [mid, hi]
        else:
            raise ValueError("bad path byte")
    return [tuple(iv) for iv in box]


# ---------------------------------------------------------------- Bernstein
def bernstein_matrix(lo, hi, n):
    """T[k][a]: the degree-n Bernstein coefficient k on [lo,hi] of the monomial x^a
    (a <= n): x = lo + w y, y^j = sum_{k>=j} C(k,j)/C(n,j) B_{k,n}(y)."""
    w = hi - lo
    T = [[fmpq(0)] * (n + 1) for _ in range(n + 1)]
    for k in range(n + 1):
        for a in range(n + 1):
            v = fmpq(0)
            for j in range(0, min(a, k) + 1):
                v += comb(a, j) * lo ** (a - j) * w ** j * fmpq(comb(k, j), comb(n, j))
            T[k][a] = v
    return T


def lcm(a, b):
    from math import gcd
    return a // gcd(a, b) * b


def to_int_matrix(T):
    """Scale a rational matrix by the positive lcm of its denominators."""
    L = 1
    for row in T:
        for v in row:
            L = lcm(L, int(v.q))
    M = np.empty((len(T), len(T[0])), dtype=object)
    for i, row in enumerate(T):
        for j, v in enumerate(row):
            x = v * L
            assert x.q == 1
            M[i, j] = int(x.p)
    return M, L


def poly_to_int_tensors(polys, degs):
    """Dense integer tensors (A, B) of shape (len(polys), degs+1...) with
    sum_e (A[e] + B[e] sqrt5) x^e = K * poly for one common positive K."""
    L = 1
    for p in polys:
        for e, c in p.items():
            for i in range(NV):
                if e[i] > degs[i]:
                    raise ValueError("degree exceeds requested")
            L = lcm(L, int(c.a.q))
            L = lcm(L, int(c.b.q))
    shape = (len(polys),) + tuple(d + 1 for d in degs)
    A = np.zeros(shape, dtype=object)
    B = np.zeros(shape, dtype=object)
    A[...] = 0
    B[...] = 0
    for k, p in enumerate(polys):
        for e, c in p.items():
            a, b = c.a * L, c.b * L
            A[(k,) + e] = int(a.p)
            B[(k,) + e] = int(b.p)
    return A, B, L


def bernstein_int(A, B, degs, box):
    """Tensor Bernstein coefficients (up to one common positive factor) of the
    polynomials encoded by the integer tensors A, B (leading axis = polynomial
    index) on the box, with degree degs[m] in variable m.
    Returns (A', B', scale): true coefficient = (A' + B' sqrt5) / (scale * K)."""
    scale = 1
    for m in range(NV):
        T, L = to_int_matrix(bernstein_matrix(box[m][0], box[m][1], degs[m]))
        scale *= L
        A = np.moveaxis(np.tensordot(A, T, axes=([1 + m], [1])), -1, 1 + m)
        B = np.moveaxis(np.tensordot(B, T, axes=([1 + m], [1])), -1, 1 + m)
    return A, B, scale


def neg_mask(A, B):
    """Elementwise exact test A + B sqrt5 < 0 for integer object arrays."""
    A = np.asarray(A, dtype=object)
    B = np.asarray(B, dtype=object)
    Af = A.ravel()
    Bf = B.ravel()
    out = np.empty(Af.shape, dtype=bool)
    for i in range(Af.size):
        out[i] = sign_q5(Af[i], Bf[i]) < 0
    return out.reshape(A.shape)


def all_negative(A, B):
    Af = np.asarray(A, dtype=object).ravel()
    Bf = np.asarray(B, dtype=object).ravel()
    for i in range(Af.size):
        if sign_q5(Af[i], Bf[i]) >= 0:
            return False
    return True


def poly_degrees(p):
    d = [0] * NV
    for e in p:
        for i in range(NV):
            d[i] = max(d[i], e[i])
    return d
