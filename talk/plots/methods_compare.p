#! /usr/bin/env gnuplot
clear
reset
set key on
set border 3
set grid noxtics
set grid 
set logscale y
set ytics nomirror
set xtics nomirror
unset border
set border 2
unset xtics
set boxwidth 0.2
set style fill solid border -1

# Output
set term pngcairo enhanced font "Serif, 14"
set output "methods_compare.png"

plot\
"methods_compare.dat" t "compare methods" with boxes
