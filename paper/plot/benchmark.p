#! /usr/bin/env gnuplot
clear
reset
set term wxt enhanced font "Serif, 14"
set key on
set border 3
set grid

set key left
set key opaque
set xlabel "t"
set ylabel "gain"

plot\
"benchmark.dat" u (($1*$1) * 1.026)t "benchmark" with linespoints

# Output
set term pngcairo enhanced font "Serif, 14"
set output "benchmark.png"
replot
set term postscript
set output "benchmark.ps"
replot
set term x11