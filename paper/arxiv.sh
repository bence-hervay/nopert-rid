#!/bin/sh
# Writes arxiv.tar.gz, the article's source as arXiv takes it: main.tex with
# its compiled bibliography, the sections, the appendices and the figures as
# PDF. Run after `latexmk main.tex`; the archive is the same for the same
# sources.
set -e
cd "$(dirname "$0")"
rm -rf build/arxiv
mkdir -p build/arxiv/sections build/arxiv/appendices build/arxiv/figures
cp main.tex build/main.bbl build/arxiv/
cp sections/*.tex build/arxiv/sections/
cp appendices/*.tex build/arxiv/appendices/
cp figures/*.pdf build/arxiv/figures/
(cd build/arxiv && tar --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner -cf - *) | gzip -n > arxiv.tar.gz
