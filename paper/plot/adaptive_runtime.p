#! /usr/bin/env gnuplot
set xlabel "1 / max. MSE"
set grid xtics
set grid ytics
set ylabel "runtime[s]"
set yrange [0:500]
set xrange [1:510]
set logscale x
unset logscale y

plot \
"adaptive_runtime.dat" u 4:3 w linespoints t "runtime adaptive" lw 3 ps 2 pt 57,\
"adaptive_runtime.dat" u 4:2 w linespoints t "runtime non adaptive" lw 3 ps 2 pt 37

# Output
set term png
set output "adaptive_runtime.png"
replot
set term postscript
set output "adaptive_runtime.ps"
replot
set term x11
