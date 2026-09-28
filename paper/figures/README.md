# Figures of the article

The figures of the configuration space, the shadows and the zooms are drawn
by the figure tool of the program (`proof/figures`) from the JSON
descriptions in [`source/`](source/), as SVG, and converted to PDF with
`rsvg-convert` for LaTeX. The drawing is exact: positions come from the
crate's exact arithmetic, Bernstein coefficients from its own conversion,
slices from its components and the crossing figure from the cover file.
Floating point only chooses what to show: the examples' enlargement windows
and magnifications. The only text is the magnification in the corner of an
enlargement; the captions explain the colours.

All of these figures are written by one command, from `proof/figures`
(with `rsvg-convert` on `PATH`):

```sh
cargo build --release
target/release/paper-figures ../../paper/figures            # every figure
target/release/paper-figures ../../paper/figures global     # or some of them
```

It writes each description to `source/`, subdivides the slices a figure
needs into `source/partitions/` (a saved partition is reused while its slice
and the crate's policy fingerprint are unchanged), and renders `<name>.svg`
and `<name>.pdf`. `RID_FIGURES_PREVIEW=<directory>` also writes PNG
previews there.

| File | Section | Content |
| --- | --- | --- |
| `hierarchy` | 4 | two slices (`t = 0`, `r₁ = r₂ = 0` in `(s, r₃)`; the plane Π₊ in `(s, θ)`) subdivided with Domain and Global, then with Local, then with the whole collection |
| `bernstein` | 5.2 | `−1 + 3x − 3x²` on `[0, 1]`, `[0, 1/2]`, `[1/2, 1]` |
| `global` | 6 | midpoints of the Global records `00101000111101011` (edge `[1,10]`, plug vertex 11) and `01000101010110001111` (edge `[0,9]`, plug vertex 10), with their supports |
| `global-refinement` | 6 | Domain and Global alone around the aligned configurations `(s, t) = (1/4, 1/10)` in `(s, r₃)` and `(1/5, 1/8)` in `(t, r₁)`, in windows of half-width 1/8, 1/32, 1/128 |
| `local` | 7 | `[1/4, 1/20, (1/1000, 1/2000, −1/1000)]` and `[2/5, 1/10, (−1/1000, 1/1000, 1/2000)]` |
| `local-refinement` | 7 | the slices of `global-refinement` with Local added |
| `exotic-directions` | 7 | the plane `r = 0` with Domain, Global and Local, and its part around the pentagon view |
| `square-pentagon` | 8 | the half-turn configuration `[1/16, 1/16, (−1/16, 1/16, 0)]`; `[1/64, 1/128, (1/1000, −1/1000, −1/1000)]`; the pentagon views `(η, ξ) = (±1/200, −1/400)` with `r = (1/1000, 0, 0)` |
| `toy-zoom` | 8 | the model `p = y(x + y)`: cells of a face zoom and of the point zooms with `ρ₀ = 1`, coloured by the smallest value of the quotient |
| `pentagon-support` | 8 | the pentagon cover's views in `(η, ξ)` and the two pentagon examples with a support on each side of `η = 0` |
| `arcs-remaining` | 8 | the plane `r = 0` and the plane Π₊ with the square and pentagon covers added |
| `arcs` | 9 | first-arc configurations of Π₊ at `e = 1/10` and `3/20` (exact, in Q(√5)) |
| `arcs-residual` | 9 | the plane Π₊ with the arc covers added, and its parts around the endpoint and the crossing |
| `slice-arcs` | 9 | slice through Π₊ in `(s, θ)` with the whole collection |
| `endpoint-crossing` | 10 | Π₊ at `e = e*`, `e = e* + 1/256`, `e = 0`, and `(e, θ) = (1/200, −1/200)` |
| `crossing-plane` | 10 | the crossing cover `crossing+` where it meets Π₊, in `(e, θ)`, and the zoom through which the second arc runs, in `(20μ, ρ)` |
| `slice-aligned` | 11 | slice `t = 1/10`, `r₁ = r₂ = 0` with the whole collection |
| `local-tiling` | B | the views of the Local covers' rows, the square cover and the pentagon cover |
