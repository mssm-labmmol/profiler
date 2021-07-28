all: doc.pdf

doc.pdf: doc.tex doc.bib
	pdflatex doc.tex
	bibtex doc
	pdflatex doc.tex
	pdflatex doc.tex
