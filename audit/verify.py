"""Exact re-verification of certificate records with the independent model.

Domain  {"inequality": k}: c_k > 0 on the whole closed box, by
        (A) tensor Bernstein coefficients of c_k (natural degrees) all > 0, and
        (B) an exact corner argument (affine: minimum at a corner; folds: the
            concavity rule proved in the report). Accepted if either proves it.
Global  {"edge":[a,b],"vertex":p}: with n(u) = u x (v_b - v_a),
        (S) shortcut: n(u).(v_a - w) >= 0 for all 60 w at the 4 corner views, and
            all 108 Bernstein coefficients (degrees (1,1,2,2,2)) of
            g = (1+|r|^2) n.v_a - n.Rhat(r) v_p strictly negative;
        (M) maximum form: for every hole vertex w, all 108 Bernstein coefficients of
            g_w = (1+|r|^2) n.w - n.Rhat(r) v_p strictly negative.
        Accepted if (S) or (M).
Local / Exotic {"cover": c}: the box lies in the covered set computed from the
        cover file's parameters (data/), exactly.
"""
from pathlib import Path
PROOF = Path(__file__).resolve().parent.parent / "proof"
import json
import sys
import time
from functools import lru_cache
from flint import fmpq
import numpy as np
import model as M
from model import Q5, q5

DATA = str(PROOF / "data")
V = M.VERTS
DP = M.domain_polynomials()
DEGS_G = [1, 1, 2, 2, 2]


# ------------------------------------------------------------------ Domain
@lru_cache(maxsize=None)
def domain_tensor(k):
    p = DP[k]
    degs = M.poly_degrees(p)
    A, B, K = M.poly_to_int_tensors([p], degs)
    return A, B, degs


def domain_bernstein(box, k):
    A, B, degs = domain_tensor(k)
    bA, bB, _ = M.bernstein_int(A, B, degs, box)
    # c_k > 0 everywhere iff (sufficient) all coefficients > 0 ; test -c_k < 0
    return M.all_negative(-bA, -bB)


def affine_min(p, box):
    """Exact minimum over the box of an affine polynomial (coordinatewise choice)."""
    const = Q5(0)
    mn = Q5(0)
    for e, c in p.items():
        if sum(e) == 0:
            mn = mn + c
            continue
        assert sum(e) == 1
        i = e.index(1)
        lo, hi = box[i]
        mn = mn + (c * lo if c.sign() >= 0 else c * hi)
    return mn


@lru_cache(maxsize=None)
def fold_L(k):
    g = M.SYMS[k - 13]
    return M.poly_quat_scalar_0u_1r_g(g)


def domain_corner(box, k):
    p = DP[k]
    if k <= 12:
        return affine_min(p, box).sign() > 0
    L = fold_L(k)
    sigma = None
    for s in box[0]:
        for t in box[1]:
            # L at fixed view is affine in r: substitute s, t
            Lr = {}
            for e, c in L.items():
                v = c * (s ** e[0]) * (t ** e[1])
                e2 = (0, 0) + e[2:]
                Lr[e2] = Lr.get(e2, Q5(0)) + v
            lo = affine_min(Lr, box)
            hi = -affine_min({e: -c for e, c in Lr.items()}, box)
            if lo.sign() > 0:
                sg, m = 1, lo
            elif hi.sign() < 0:
                sg, m = -1, -hi
            else:
                return False
            if sigma is None:
                sigma = sg
            elif sigma != sg:
                return False
            N = 1 + s * s + t * t
            if not (m * m - N).sign() > 0:
                return False
    return True


def verify_domain(box, k):
    a = domain_bernstein(box, k)
    b = domain_corner(box, k)
    return {"bernstein": a, "corner": b, "ok": a or b}


# ------------------------------------------------------------------ Global
@lru_cache(maxsize=4096)
def global_tensor(a, b, p):
    d = M.vsub(V[b], V[a])
    n = M.pcross(M.U, [M.pconst(c) for c in d])
    nRp = M.pdot(n, M.cayley_hat_apply(V[p]))
    E = M.padd(M.ONE, M.RR)
    members = []
    for w in V:
        nw = M.pdot(n, [M.pconst(c) for c in w])
        members.append(M.psub(M.pmul(E, nw), nRp))
    A, B, K = M.poly_to_int_tensors(members, DEGS_G)
    return A, B


def valid_corners(box, a, b):
    d = M.vsub(V[b], V[a])
    for s in box[0]:
        for t in box[1]:
            n = M.cross([q5(s), q5(t), Q5(1)], d)
            for w in V:
                if M.dot(n, M.vsub(V[a], w)).sign() < 0:
                    return False
    return True


def verify_global(box, a, b, p, want_both=True):
    A, B = global_tensor(a, b, p)
    bA, bB, _ = M.bernstein_int(A, B, DEGS_G, box)
    gap_neg = M.all_negative(bA[a], bB[a])  # member v_a is the edge's gap
    valid = valid_corners(box, a, b) if gap_neg or want_both else None
    shortcut = bool(gap_neg and valid)
    maximum = None
    if want_both or not shortcut:
        maximum = gap_neg and M.all_negative(bA, bB)  # all 60 members
        if not gap_neg:
            maximum = False
    return {"shortcut": shortcut, "valid": valid, "gap_negative": gap_neg, "maximum": maximum,
            "ok": bool(shortcut or maximum)}


# ------------------------------------------------------------------ covers
def num(x):
    """["a","b"] -> a + b sqrt5 ; "a" -> rational."""
    if isinstance(x, list):
        return Q5(fmpq(x[0]) if "/" not in x[0] else fmpq(*map(int, x[0].split("/"))),
                  fmpq(x[1]) if "/" not in x[1] else fmpq(*map(int, x[1].split("/"))))
    if isinstance(x, str):
        return q5(fmpq(*map(int, x.split("/"))) if "/" in x else fmpq(int(x)))
    return q5(fmpq(x))


def inverse(Mt):
    n = len(Mt)
    A = [[Mt[i][j] for j in range(n)] + [Q5(1) if i == j else Q5(0) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if not A[r][col].iszero())
        A[col], A[piv] = A[piv], A[col]
        iv = A[col][col].inv()
        A[col] = [x * iv for x in A[col]]
        for r in range(n):
            if r != col and not A[r][col].iszero():
                f = A[r][col]
                A[r] = [A[r][j] - f * A[col][j] for j in range(2 * n)]
    return [row[n:] for row in A]


A1 = Q5(fmpq(-1, 2), fmpq(3, 10))  # (-5+3 sqrt5)/10
A3 = Q5(fmpq(-1, 2), fmpq(1, 10))  # (-5+sqrt5)/10


class Cover:
    def __init__(self, path):
        d = json.load(open(path))
        self.name = d["name"]
        P = d["parameters"]
        co = P["coordinates"]
        if co == "configuration":
            self.sigma = None
        else:
            self.sigma = {"plus": 1, "minus": -1}[co["arc-plane"]]
        self.c0 = [num(x) for x in P["centre"]]
        self.M = [[num(x) for x in row] for row in P["map"]]
        self.Minv = inverse(self.M)
        # check inverse
        for i in range(5):
            for j in range(5):
                s = Q5(0)
                for k in range(5):
                    s = s + self.M[i][k] * self.Minv[k][j]
                assert s == (Q5(1) if i == j else Q5(0))
        sh = P["shape"]
        iv = []
        if "tube" in sh:
            T = sh["tube"]
            for lo, hi in T["base"]:
                iv.append([num(lo), num(hi)])
            eps = num(T["radius"])
            for lo, hi in T["offset"]:
                iv.append([eps * lo, eps * hi])
        else:
            Pt = sh["point"]
            radii = [num(x) for x in Pt["radii"]]
            faces = []
            for j, (lo, hi) in enumerate(Pt["base"]):
                assert lo in (-1, 0) and hi in (0, 1) and lo < hi
                if lo == -1:
                    faces.append((j, -1))
                if hi == 1:
                    faces.append((j, 1))
            assert len(faces) == len(radii)
            rad = dict(zip(faces, radii))
            for j, (lo, hi) in enumerate(Pt["base"]):
                iv.append([rad[(j, -1)] * lo if lo == -1 else Q5(0), rad[(j, 1)] * hi if hi == 1 else Q5(0)])
            eps = num(Pt["radius"])
            for lo, hi in Pt["offset"]:
                iv.append([eps * lo, eps * hi])
        assert len(iv) == 5
        self.iv = iv
        self.beyond = None
        if P["beyond"] is not None:
            self.beyond = P["beyond"]["axis"]
            self.beyond_ineq = P["beyond"]["inequality"]

    # coordinates x(s,t,r) of the cover's coordinate system
    def view_coords(self, s, t):
        if self.sigma is None:
            return [q5(s), q5(t)]
        return [q5((1 - 2 * s) / (2 + s)), q5(5 * t / (2 + s))]

    def rot_affine(self):
        """x2,x3,x4 as affine functions of r: list of (const, [c1,c2,c3])."""
        if self.sigma is None:
            return [(Q5(0), [Q5(1 if k == j else 0) for k in range(3)]) for j in range(3)]
        sg = self.sigma
        return [
            (Q5(0), [Q5(0), Q5(1), Q5(0)]),                       # theta = r2
            (-(sg * A1), [Q5(1), -(sg * A3), Q5(0)]),            # zeta1 = r1 - sg(a1 + a3 r2)
            (-(sg * A3), [Q5(0), sg * A1, Q5(1)]),               # zeta3 = r3 - sg(a3 - a1 r2)
        ]

    def z_ranges(self, box):
        corners = [self.view_coords(s, t) for s in box[0] for t in box[1]]
        rot = self.rot_affine()
        out = []
        for i in range(5):
            row = self.Minv[i]
            shift = Q5(0)
            for j in range(5):
                shift = shift + row[j] * self.c0[j]
            # view part: extremes at the four (s,t) corners (monotone in s and in t separately)
            vals = [row[0] * c[0] + row[1] * c[1] for c in corners]
            vlo = vals[0]
            vhi = vals[0]
            for x in vals[1:]:
                if x < vlo:
                    vlo = x
                if x > vhi:
                    vhi = x
            # rotation part: const + sum c_k r_k, independent r_k
            const = Q5(0)
            coef = [Q5(0)] * 3
            for j in range(3):
                cj, lin = rot[j]
                const = const + row[2 + j] * cj
                for k in range(3):
                    coef[k] = coef[k] + row[2 + j] * lin[k]
            rlo = const
            rhi = const
            for k in range(3):
                lo, hi = box[2 + k]
                a, b = coef[k] * lo, coef[k] * hi
                if a <= b:
                    rlo, rhi = rlo + a, rhi + b
                else:
                    rlo, rhi = rlo + b, rhi + a
            out.append((vlo + rlo - shift, vhi + rhi - shift))
        return out

    def contains(self, box):
        rng = self.z_ranges(box)
        for i in range(5):
            lo, hi = rng[i]
            L, H = self.iv[i]
            if lo < L:
                return False
            if self.beyond == i:
                continue
            if hi > H:
                return False
        return True


COVERS = {}


def cover(name):
    if name not in COVERS:
        if isinstance(name, int):
            COVERS[name] = Cover(f"{DATA}/local/{name}.json")
        else:
            COVERS[name] = Cover(f"{DATA}/exotic/{name}.json")
        assert str(COVERS[name].name) == str(name)
    return COVERS[name]


# ------------------------------------------------------------------ driver
def verify_record(txt):
    obj = json.loads(txt)
    box = M.box_of_path(obj["path"])
    c = obj["component"]
    t0 = time.time()
    if c == "Domain":
        res = verify_domain(box, obj["inequality"])
    elif c == "Global":
        w = obj.get("witness", obj)  # the flat record of 2026-09-28: {"edge":[a,b],"vertex":p}
        res = verify_global(box, w["edge"][0], w["edge"][1], w["vertex"])
    else:
        res = {"ok": cover(obj["cover"]).contains(box)}
    res["path"] = obj["path"]
    res["component"] = c
    res["sec"] = round(time.time() - t0, 4)
    return res


if __name__ == "__main__":
    import random
    lines = open(str(PROOF / "results" / "full" / "search.cert")).read().split("\n")[1:-1]
    rng = random.Random(1)
    for comp in ["Domain", "Global", "Local", "Exotic"]:
        sel = [l for l in lines if f'"component":"{comp}"' in l]
        for l in rng.sample(sel, 5):
            txt = l[9:]
            print(verify_record(txt))
