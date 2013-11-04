#! /usr/bin/env gnuplot
set ylabel "max. MSE"
set grid xtics ytics
set xtics nomirror
set ytics nomirror
set xlabel "runtime[s]"
set yrange [0.001:10]
set xrange [4:131]
set logscale x
set logscale y

plot \
"adaptive_runtime.dat" u 3:1 t "runtime adaptive" w linespoints lw 3 ps 2 pt 57,\
"adaptive_runtime.dat" u 2:1 t "runtime non adaptive" w linespoints lw 3 ps 2 pt 37


# Output
set term png
set output "adaptive_runtime.png"
replot
set term postscript
set output "adaptive_runtime.ps"
replot
set term wxt
