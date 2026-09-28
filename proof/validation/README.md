# Validation experiments

Experiments on the components of `rid`. They are evidence about the
implementation, not part of the proof: they record what the components do on
chosen configurations.

- **completeness**: for each of 683 exact configurations near the special
  configurations of the article (`catalogue/points.json`), the largest tested
  radius `2^-k` at which the collection eliminates every box around it;
- **irredundancy**: for one target configuration per component and per zoom
  cover (`catalogue/targets.json`), the scales at which only that component or
  cover eliminates the boxes around it;
- **pilot**: bounded searches below chosen boxes, to estimate the cost of a
  full search;
- **check** and **compare**: consistency checks of one transcript, and
  changes between two.

From the `proof/` directory:

```sh
cargo build --release --manifest-path validation/Cargo.toml
validation/target/release/rid-validation completeness validation/examples/completeness.json
validation/target/release/rid-validation irredundancy validation/examples/irredundancy.json
```

Each command prints a JSON summary and writes a transcript, one JSON line per
configuration. For the final collection, every one of the 683 configurations
is resolved, and each of the 45 targets is needed.
