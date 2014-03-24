#! /usr/bin/env gnuplot
clear
reset
set key opaque
set key Left top left

set xtics nomirror
set logscale x
set logscale y
set ylabel "runtime[s]"
set grid ytics xtics
set xlabel "rays per sample"
set format y "10^{%T}"
set format x "10^{%T}"
set xrange [1e4:1e8]
set yrange [1:1e6]

# Output
set term pngcairo enhanced font 'Serif,14'
set output "runtime.png"

plot \
"runtime_cylindrical.dat" u 1:5 w linespoints axes x1y1 t "  1 x CPU IS "  lw 5 ps 1.5 pt 5 lc rgb "blue" lt 1,\
"runtime_cylindrical.dat" u 1:4 w linespoints axes x1y1 t "  1 x GPU IS"  lw 5 ps 1.5 pt 7 lc rgb "black" lt 1,\
"runtime_cylindrical.dat" u 1:3 w linespoints axes x1y1 t "  4 x GPU IS"  lw 5 ps 1.5 pt 11 lc rgb "violet" lt 1,\
"runtime_cylindrical.dat" u 1:2 w linespoints axes x1y1 t "47 x GPU IS"  lw 5 ps 1.5 pt 13 lc rgb "cyan" lt 1

# Output as pdf
set terminal pdfcairo enhanced font "Serif, 14"
set output "runtime.pdf"
replot
