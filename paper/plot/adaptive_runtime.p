#! /usr/bin/env gnuplot
clear
reset
set ylabel "max. MSE"
set grid xtics ytics
set xtics nomirror
set ytics nomirror
set xlabel "runtime[s]"
set xrange [15:3084]
set logscale x
set logscale y

# Output
set term png enhanced font "Serif, 14"
set output "adaptive_runtime.png"

plot \
"adaptive_runtime.dat" u 2:1 t "IS" w linespoints lc rgb "blue" lw 3 ps 2 pt 37,\
"adaptive_runtime.dat" u 4:3 t "IS + AS + RS" w linespoints lc rgb "orange" lw 3 ps 2 pt 57


