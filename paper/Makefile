all: *.tex
	#convert graphics/delauny_3.svg graphics/delauny_3.png
	pdflatex ase_flux.tex
	#bibtex ase_flux 
	#pdflatex ase_flux.tex
	pdflatex ase_flux.tex

clean:
	rm -f *.aux
	rm -f *.toc
	rm -f *.bbl
	rm -f *.blg
	rm -f *.out
	rm -f *.bcf
	rm -f *.log
	rm -f *.run.xml

graphics/%.png: graphics/%.svg
	convert $< -format "png" -
