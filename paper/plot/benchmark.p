#! /usr/bin/env gnuplot
clear
reset
set key on
set border 3
set grid

set key Left top right
#set key reverse
set xlabel "time[µs]"
set ylabel "gain"
set xtics 500
set ytics 0,1,5
set yrange [0:5]
set xrange [0:1500]
set xtics nomirror
set ytics nomirror

# Output
set term pngcairo enhanced font "Serif, 14"
set output "benchmark.png"

plot\
"benchmark_polychromatic.dat" u 1:(($2*$2) * 1.0263) t "(a)" with lines lw 2,\
"exp_ceramics.dat" u 1:2 t "(b)" with lines lw 2,\
"benchmark_polychromatic_refl.dat" u 1:(($2*$2) * 1.0263) t "(c)" with lines lw 2,\
"benchmark_100k.dat" u 1:(($2*$2) * 1.0263) t "(d)" with lines lw 2,\
"benchmark_refl.dat" u 1:(($2*$2) * 1.0263) t "(e)" with lines lw 2

#a  polychromatic, no reflection
#b  experimental measurement
#c  polychromatic, with reflection
#d  monochromatic, no reflection
#e  monochromatic, with reflection

   
