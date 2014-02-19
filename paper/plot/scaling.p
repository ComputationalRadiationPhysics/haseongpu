#! /usr/bin/env gnuplot
set term wxt enhanced font "Serif, 14"
set key top left
set grid xtics
set grid ytics
set xlabel "GPUs"
set ylabel "speedup"
#set ylabel "efficency"
#set xtics 0,10,70
#set ytics 0,8,70
#set ytics 0,4,64
set yrange [0:70]
#set yrange [0:1]
set xrange [0:70]
#set yrange [0:1]


f(x)=x


plot \
"scaling.dat" u 1:(7858/$3) w linespoints t "IS"  lc rgb "blue" lw 3 ps 2 pt 57,\
"scaling.dat" u 1:(13619/$4) w linespoints t "IS + AS + RS"  lc rgb "orange" lw 3 ps 2 pt 57,\
f(x) t "" lc rgb "black"

# Output
set term pngcairo enhanced font "Serif, 14"
set output "scaling.png"
replot
set term postscript
set output "scaling.ps"
replot
set term x11