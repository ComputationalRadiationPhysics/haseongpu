#! /usr/bin/env gnuplot
clear
reset
set term wxt enhanced font "Serif, 14"
set ylabel "max. MSE"
set grid xtics ytics
set xtics nomirror
set ytics nomirror
set xlabel "runtime[s]"
#set yrange [0.001:10]
set xrange [15:3084]
set logscale x
set logscale y

plot \
"adaptive_runtime.dat" u 2:1 t "IS" w linespoints lc rgb "blue" lw 3 ps 2 pt 37,\
"adaptive_runtime.dat" u 4:3 t "IS + AS + RS" w linespoints lc rgb "orange" lw 3 ps 2 pt 57



# Output
set term png enhanced font "Serif, 14"
set output "adaptive_runtime.png"
replot
set term postscript
set output "adaptive_runtime.ps"
replot
set term wxt
