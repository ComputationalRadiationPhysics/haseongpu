#! /usr/bin/env gnuplot
clear
reset
set key opaque
set key Left top left

set title "Emission and absorption spectrum"
set ylabel "intensity[cm²]"
set format y "%.1t · 10^{%T}"
set xtics 910,50
set grid ytics xtics
set xlabel "wavelength[nm]"

# Output
set output "sigma.png"
set term pngcairo enhanced font 'Serif,14'

plot\
"sigma.dat" u 1:3 t "Emission" w lines lw 5 lt 1 lc rgb "green",\
"sigma.dat" u 1:2 t "Absorption" w lines lw 5 lt 1 lc rgb "red"

# Output as pdf
set terminal pdfcairo enhanced font "Serif, 14"
set output "sigma.pdf"
replot

