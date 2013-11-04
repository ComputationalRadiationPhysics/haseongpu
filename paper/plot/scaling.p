#! /usr/bin/env gnuplot
set key top left
set grid xtics
set grid ytics
set xlabel "GPUs"
set ylabel "speedup"
#set ylabel "efficency"
set xtics 1,1,4
set ytics 1,0.5,4
set yrange [1:4]
#set yrange [0:1]

plot \
"scaling.dat" u 1:2 w linespoints t "scaling adaptive" lw 3 ps 2 pt 57,\
"scaling.dat" u 1:3 w linespoints t "scaling non adaptive" lw 3 ps 2 pt 37


# Output
set term png
set output "scaling.png"
replot
set term postscript
set output "scaling.ps"
replot
set term x11