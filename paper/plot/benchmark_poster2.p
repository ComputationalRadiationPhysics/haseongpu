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
set ytics 0,1,7
set yrange [0:7]
set xrange [0:1500]
#set xrange [0:3000]
set xtics nomirror
set ytics nomirror

# Output
set term pngcairo enhanced font "Serif, 14"
set output "benchmark_poster2.png"

plot\
"exp_ceramics.dat" u 1:2 t "experimental measurement" with lines lw 4 lt rgb "green",\
"benchmark_no_ase.dat" u 1:(($2*$2) * 1.0263) t "no ASE" with lines lw 4 lt rgb "blue",\
"benchmark_100k.dat" u 1:(($2*$2) * 1.0263) t "simplified ASE simulation" with lines lw 4 lt rgb "violet"

#"benchmark_cladding.dat" u 1:2 t "(c)" with lines lw 4,\
#"benchmark_polychromatic_refl.dat" u 1:(($2*$2) * 1.0263) t "(c)" with lines lw 4,\


set terminal pdfcairo enhanced font "Serif, 14"
set output "benchmark_poster2.pdf"
replot

#a  polychromatic, no reflection
#b  experimental measurement
#c  polychromatic, with reflection
#d  monochromatic, no reflection
#e  monochromatic, with reflection

   
