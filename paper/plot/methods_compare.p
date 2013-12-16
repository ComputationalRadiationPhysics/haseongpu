#! /usr/bin/env gnuplot
clear
reset
set term wxt enhanced font "Serif, 14"
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


plot\
"methods_compare.dat" t "compare methods" with boxes

# Output
set term pngcairo enhanced font "Serif, 14"
set output "methods_compare.png"
replot
set term postscript
set output "methods_compare.ps"
replot
set term svg
set output "methods_compare.svg"
replot
set term x11