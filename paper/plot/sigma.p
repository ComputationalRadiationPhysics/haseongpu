#! /usr/bin/env gnuplot
clear
reset
set key opaque
set key Left top left

set title "Emission and absorption spectrum"
set ylabel "intensity[?]"
set format y "%.1t · 10^{%T}"
set xtics 905,45,1095
set grid ytics xtics
set xlabel "wavelength[nm]"

plot\
"sigma.dat" u 1:3 t "Emission" w lines lw 5 lt 1 lc rgb "green",\
"sigma.dat" u 1:2 t "Absorption" w lines lw 5 lt 1 lc rgb "red"

# Output
set term pngcairo enhanced font 'Serif,14'
set output "sigma.png"
replot
set term x11
