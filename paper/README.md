# The article

*Non-Rupertness of the Rhombicosidodecahedron*:
[`article.pdf`](article.pdf). To build it, with a TeX distribution that
provides `latexmk` (which also copies the result from `build/main.pdf` to
`article.pdf`):

```sh
latexmk main.tex
./arxiv.sh        # then: arxiv.tar.gz, the source as submitted to arXiv
```

The figures in [`figures/`](figures/) are drawn by the figure tool of the
program; [`figures/README.md`](figures/README.md) explains how each was made.
