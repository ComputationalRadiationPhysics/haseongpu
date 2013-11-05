#! /usr/bin/env gnuplot
set key top right
set grid xtics
set grid ytics
set xlabel "GPUs"
#set ylabel "speedup"
set ylabel "efficency"
set xtics 8,8,64
#set ytics 0,4,64
#set yrange [1:64]
set yrange [0:1]
set xrange [1:64]
#set yrange [0:1]


plot \
"scaling.dat" u 1:((472/$3)/$1) w linespoints t "efficency adaptive" lw 3 ps 2 pt 57,\
"scaling.dat" u 1:((1057/$2)/$1) w linespoints t "efficency non adaptive" lw 3 ps 2 pt 37


# Output
set term png
set output "scaling.png"
replot
set term postscript
set output "scaling.ps"
replot
set term x11