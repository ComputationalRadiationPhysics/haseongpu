#! /usr/bin/env gnuplot
clear
reset
set key on
set border 3
set grid

set key out
set key top center
#set key width -6
set key reverse
set xlabel "time[µs]"
set ylabel "gain"
set xtics 500
set ytics 0,1,5
set yrange [0:5]
set xrange [0:2000]
set xtics nomirror
set ytics nomirror

plot\
"exp_crystal.dat" u 1:2 t "experiment measurement crystal" with lines ,\
"exp_ceramics.dat" u 1:2 t "experiment measurement ceramics" with lines ,\
"benchmark_100k.dat" u 1:(($2*$2) * 1.0263) t "monochromatic, no reflection" with lines,\
"benchmark_polychromatic_refl.dat" u 1:(($2*$2) * 1.0263) t "polychromatic, with reflection" with lines,\
"benchmark_refl.dat" u 1:(($2*$2) * 1.0263) t "monochromatic, with reflection" with lines,\
"benchmark_polychromatic.dat" u 1:(($2*$2) * 1.0263) t "polychromatic, no reflection" with lines

# Output
set term pngcairo enhanced font "Serif, 14"
set output "benchmark.png"
replot
set term postscript
set output "benchmark.ps"
replot
set term x11