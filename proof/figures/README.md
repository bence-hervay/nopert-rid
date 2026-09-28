# Figures

`rid-figures` draws the article's figures from JSON descriptions, using the
exact geometry and the components of `rid`:

- the shadows of the two copies at a configuration, with enlarged insets;
- envelopes of every vertex and edge over a box of configurations;
- slices of the configuration space, each piece coloured by the component
  that eliminates it;
- a support at a configuration: the supporting line of the hole's shadow,
  the silhouette on it, the outward normal and a vertex of the plug beyond
  it;
- the graph of a polynomial with its Bernstein coefficients;
- explicit exact segments, polygons and dots, for guide lines and diagrams;
- a small label in the corner of a panel, such as an enlargement's
  magnification.

A slice names the parts of the collection it is subdivided with (some of
the components and Exotic covers), so a figure can show what a partial
collection leaves unresolved.

```sh
cargo build --release
target/release/rid-figures <description.json> <output-directory>
```

The output directory receives `figure.svg` and a PNG preview. The
article's figures are written, all at once, by `paper-figures` (see
[`../../paper/figures`](../../paper/figures)). Every field of
a description is required; numbers are exact rationals or pairs `[a, b]` for
`a + b√5`. [`examples/`](examples/) has one description of each kind of
scene, and the article's own descriptions are in
[`../../paper/figures/source`](../../paper/figures/source).
