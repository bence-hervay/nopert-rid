# The program

`rid` searches for and checks a certificate that no configuration of the
rhombicosidodecahedron (RID) is a fit. All decisions use exact arithmetic in
Q(√5); floating point is used only to choose which witness to try, never to
accept one. The mathematics is explained in the article in [`../paper`](../paper).

## Commands

```sh
cargo build --release
target/release/rid prepare examples/prepare.json   # load and verify the zoom covers
target/release/rid search  examples/search.json    # search the whole domain, writing rid.cert
target/release/rid check   examples/check.json     # re-verify a certificate completely
target/release/rid generate <config.json>          # recompute zoom cover files
```

Every configuration is a JSON file whose fields are all required. A search
can be interrupted at any time and resumed with the same command: the
certificate is written as it goes and is the only saved state. The
certificate does not depend on the number of threads.

The final certificate and the configurations that produced and checked it
are in [`results/full/`](results/full/).

## Source

The source is organised in five areas, each depending only on the ones
before it:

| Area | Contents |
| --- | --- |
| `arithmetic` | exact numbers `a + b√5` with rational `a, b`; polynomials in five variables; tensor Bernstein coefficients and the sign test |
| `problem` | the RID's vertices and symmetries; configurations and search boxes; the 73 inequalities of the reduced domain `D` |
| `elimination` | witnesses (support gaps, the all-vertex form, inequalities of `D`); zooms, zoom cells and zoom covers; the zoom lemma; loading and checking cover files; the cover generator |
| `components` | the Domain, Exotic, Local and Global components and their records |
| `search` | the certificate format, the ordered parallel search, the configuration files and the commands |

A component proposes candidate witnesses for a box and accepts one only when
its exact test holds; a box that no component eliminates is halved.

## Zoom covers

[`data/`](data/) holds the 38 zoom covers used by the Local component (30,
around aligned configurations) and the Exotic component (8: the square
and pentagon configurations, the arcs, their endpoints and crossings). Each
file lists the cells of every zoom and the witness of each cell; the
parameters of each cover are fixed in `src/elimination/proof/catalogue`.
Every cell is checked with the zoom lemma at every start of the program, so
the files are not trusted. `rid generate` recomputes them.

## Certificate format

The first line is a header recording the format, a fingerprint of the
program's source code, data and dependencies (a certificate is only read by
the program that wrote it), the root box and the depth limit. Every further
line records one eliminated box: an 8-digit checksum (the start of the
SHA-256 hash of the rest of the line), then JSON with the box's path of
halvings and the component's data, for example

```text
651fc54e {"path":"11","component":"Domain","inequality":0}
```

## Tests

`cargo test --release` runs the unit, integration and adversarial
tests in a few minutes. Heavier tests are marked `#[ignore]`.

## Companion packages

- [`validation/`](validation/): experiments on the components, such as
  whether every configuration of a list of special configurations has a
  neighbourhood that the collection eliminates, and whether each component
  and cover is needed.
- [`figures/`](figures/): the tool that draws the article's figures from
  JSON descriptions, using the program's exact geometry and components.
