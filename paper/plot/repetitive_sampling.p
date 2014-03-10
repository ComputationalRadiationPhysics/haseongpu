#! /usr/bin/env gnuplot
clear
reset
set key on
set title "MSE threshold = 0.001 with 3210 sample points"
set xlabel "rays per sample point"
set ylabel "sample points"

set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set grid noxtics noytics
set xtics nomirror
set ytics nomirror

# Output
set term pngcairo enhanced font "Serif, 14"
set output "repetitive_sampling.png"

plot\
"repetitive_sampling.dat" u 5:xtic(1) t "IS+AS in 1160s",\
"repetitive_sampling.dat" u 4:xtic(1) t "IS+AS+RS in 600s" lt rgb "gold" 
