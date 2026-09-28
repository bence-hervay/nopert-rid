$pdf_mode = 1;
$out_dir = 'build';
$pdflatex = 'pdflatex -synctex=1 -interaction=nonstopmode -file-line-error -halt-on-error %O %S';
# The article, built, next to its source.
$success_cmd = 'cp build/main.pdf article.pdf';
