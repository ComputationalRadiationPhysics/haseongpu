#! /usr/bin/env gnuplot
clear
reset
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

# Output
set term pngcairo enhanced font "Serif, 14"
set output "scaling_poster.png"

plot \
"scaling.dat" u 1:(13619/$4) w linespoints t "extended ASE simulation"  lc rgb "orange" lw 4 ps 1.5 pt 57,\
f(x) t "" lw 4 lc rgb "black"

# Output as pdf
set terminal pdfcairo enhanced font "Serif, 14"
set output "scaling_poster.pdf"
replot
