# The rhombicosidodecahedron is not Rupert

Can a solid pass through a straight tunnel drilled through an identical copy
of itself? Prince Rupert of the Rhine wagered in the seventeenth century that
a cube can, and it was long believed that every convex polyhedron might share
this property. This repository contains a proof that the
rhombicosidodecahedron does not: however the two copies are turned, and
whichever direction the tunnel takes, the shadow of one never fits strictly
inside the shadow of the other.

The proof is computer-assisted. It consists of

- an article, [`paper/article.pdf`](paper/article.pdf), with its source in
  [`paper/`](paper/), which explains the mathematics;
- a program, [`proof/`](proof/), written in Rust, which searches for and
  checks a certificate using exact arithmetic;
- the final certificate, [`proof/results/full/`](proof/results/full/);
- an independent checker, [`audit/`](audit/), written in Python.

## How the proof works

Symmetry reduces the question to a bounded five-dimensional set of
configurations: a viewing direction and a relative rotation. This set is cut
into boxes, and each box is eliminated by one of four components, each of
which certifies with exact arithmetic in Q(√5) that no configuration in the
box is a fit:

- **Domain**: the box lies outside the symmetry-reduced domain;
- **Global**: some vertex of the rotated copy sticks out of the other's
  shadow throughout the box;
- **Local** and **Exotic**: the box lies in a precomputed region around
  configurations where the two shadows touch.

The difficulty is the touching configurations, around which no strict
inequality can hold. They are handled by one lemma: near a set of
configurations known to contain no fit, a necessary inequality is divided by
the distance from that set, and its sign is then certified on a whole box by
Bernstein coefficients, touching configurations included. The article
develops this step by step.

## Checking the proof

The program builds with a recent Rust toolchain; `Cargo.lock` pins the
exact versions of its few dependencies:

```sh
cd proof
cargo build --release
(cd results/full && sha256sum -c SHA256SUMS)
target/release/rid check results/full/check.json
```

The check rebuilds and verifies every zoom cover, re-verifies every record of
the certificate exactly, and confirms that the records cover the whole
search domain; it prints `"complete":true` on success. It takes about half a
minute on 16 threads. Running `target/release/rid search
results/full/search.json` from a copy without the certificate reproduces it
byte for byte in about a minute.

| | |
| --- | --- |
| Records | 192,696 |
| Size | 18,195,036 bytes |
| SHA-256 | `df1f8dbab68b561f22fa9daeca48153d71bde7b6fba3345d9f6b7b8077df4b5f` |
