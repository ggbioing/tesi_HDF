FILE=tesi

all : 
	make pdf
	make clean

pdf :
	latex $(FILE).tex
	bibtex $(FILE).aux
	latex $(FILE).tex
	latex $(FILE).tex
	dvipdf $(FILE).dvi $(FILE).pdf

clean :
	rm -f *.aux *.bbl *.blg *.dvi *.log *.toc
