# An independent checker

A second implementation of the certificate check, in Python, written from the
mathematical definitions without using the program's code. It re-verifies
every record of the final certificate exactly:

- Domain records, by exact Bernstein coefficients of the named inequality
  (and, independently, by the corner criterion for affine and fold
  inequalities);
- Global records, by the sixty gap polynomials of the named edge and plug
  vertex;
- Local and Exotic records, by exact inclusion of the box in the named
  zoom cover's covered set, read from the cover files.

It also checks the certificate's framing, checksums and canonical JSON, that
the paths are prefix-free and cover the whole search domain (`structure.py`),
and that deliberately wrong records are refused (`controls.py`).

It is a [uv](https://docs.astral.sh/uv/) project; `uv.lock` pins `numpy` and
`python-flint`:

```sh
uv run python runall.py      # every record, on 12 processes; writes results.jsonl
uv run python structure.py   # framing, order and coverage
uv run python controls.py    # wrong records are refused
uv run python selftest.py    # the checker's own consistency tests
```

The cell proofs of the zoom covers are checked by the program at start-up
and are not re-checked here.
