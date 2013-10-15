set key top left
set grid xtics
set grid ytics
set xlabel "GPUs"
set ylabel "speedup"
set xtics 1,1,4
set ytics 1,0.5,4
set yrange [1:4]

plot \
"scaling.dat" u 1:2 w linespoints t "scaling adaptive" pt 57,\
"scaling.dat" u 1:3 w linespoints t "scaling not adaptive" pt 37

# Output
set term png
set output "scaling.png"
replot
set term postscript
set output "scaling.ps"
replot
set term x11