#! /usr/bin/env gnuplot
#set key opaque
set key top left
set logscale x
set logscale y
set ylabel "runtime[s]"
#set ytics nomirror
#set xtics nomirror
#set grid #noxtics noytics noztics front
set grid ytics xtics
set xlabel "rays per sample"
set format x "%.0te%+03T";
set xrange [10000:100000000]

plot \
"runtime.dat" u 1:2 w linespoints axes x1y1 t "1 x CPU runtime " lw 3 ps 2 pt 57 lt 3,\
"runtime.dat" u 1:4 w linespoints axes x1y1 t "1 x GPU runtime " lw 3 ps 2 lt 7,\
"runtime.dat" u 1:3 w linespoints axes x1y1 t "4 x GPU runtime " lw 3 ps 2 pt 37 lt 5,\
"runtime.dat" u 1:5 w linespoints axes x1y1 t "48 x GPU runtime " lw 3 ps 2 lt 4



# Output
set term pngcairo enhanced font 'Serif,14'
set output "runtime.png"
replot
set term postscript
set output "runtime.ps"
replot
set term x11
