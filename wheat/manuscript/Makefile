default : frontiers_wheatsenescence2016.pdf

tgz : wheatsenescence.tgz

frontiers_wheatsenescence2016.pdf : frontiers_wheatsenescence2016.tex bioinfo.bib
	pdflatex frontiers_wheatsenescence2016
	bibtex frontiers_wheatsenescence2016
	pdflatex frontiers_wheatsenescence2016
	pdflatex frontiers_wheatsenescence2016
	pdflatex frontiers_wheatsenescence2016

wheatsenescence.tgz : frontiers_wheatsenescence2016.pdf
	rm -rf wheatsenescence
	mkdir wheatsenescence
	cp *.tex *.eps *.jpg *.cls *.bib Makefile frontiers_wheatsenescence2016.pdf wheatsenescence/
	tar -zcvf wheatsenescence.tgz wheatsenescencet

.PHONY : tgz clean
