#! /usr/bin/env gnuplot
clear
reset
set term wxt enhanced font "Serif, 14"
set key on

set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set grid noxtics noytics
set xtics nomirror
set ytics nomirror
#set format x "%.0te%+03T";

plot\
"repetitive_sampling.dat" u 3:xtic(1) t "IS+AS",\
"repetitive_sampling.dat" u 2:xtic(1) t "IS+AS+RS" lt rgb "gold" 
   

# Output
#set term pngcairo enhanced font "Serif, 14"
#set output "repetitive_sampling.png"
#replot
#set term postscript
#set output "repetitive_sampling.ps"
#replot
#set term svg
#set output "repetitive_sampling.svg"
#replot
#set term x11