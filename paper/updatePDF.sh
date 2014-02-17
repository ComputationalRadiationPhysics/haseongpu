#!/bin/bash

while [ 1 ] ; do

	make

	echo -e "\n\n\n\n\n\nWarnings:\n"
	grep  Warning ase_flux.log 
	echo -e "\n\n\nErrors:\n"
	grep  Error ase_flux.log
	echo -e "\n\n\n\n"

	notify-send -t 3000 "pdflatex" "compilation done"

	inotifywait -e modify *.tex #*.bib
done

