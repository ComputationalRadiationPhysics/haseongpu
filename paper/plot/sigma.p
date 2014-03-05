#! /usr/bin/env gnuplot
#set key opaque
set key top left
set title "Emission and absorption spectrum"
set ylabel "Intensity[?]"
#set ytics nomirror
#set xtics nomirror
set xtics 905,45,1095
#set grid #noxtics noytics noztics front
set grid ytics xtics
set xlabel "Wavelength[nm]"
#set format x "%.0te%+03T";

plot\
"sigma.dat" u 1:3 t "Emission" w lines lw 3,\
"sigma.dat" u 1:2 t "Absorption" w lines lw 3
# Output
set term pngcairo enhanced font 'Serif,14'
set output "sigma.png"
replot
set term postscript
set output "runtime.ps"
replot
set term x11
