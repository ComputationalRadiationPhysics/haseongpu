#! /usr/bin/env gnuplot
clear
reset
set term wxt enhanced font "Serif, 14"
set key on
set border 3
set grid

set key out
set key center top
#set key opaque
set xlabel "time[ns]"
set ylabel "gain"
#set xtics 0,0.1,1
set xtics nomirror
set ytics nomirror

plot\
"benchmark_100k.dat" u 1:(($2*$2) * 1.0263) t "monochromatic, no reflection" with linespoints,\
"benchmark_refl.dat" u 1:(($2*$2) * 1.0263) t "monochromatic, with reflection" with linespoints,\
"benchmark_polychromatic.dat" u 1:(($2*$2) * 1.0263) t "polychromatic, no reflection" with linespoints,\
"benchmark_polychromatic_refl.dat" u 1:(($2*$2) * 1.0263) t "polychromatic, with reflection" with linespoints

# Output
set term pngcairo enhanced font "Serif, 14"
set output "benchmark.png"
replot
set term postscript
set output "benchmark.ps"
replot
set term x11