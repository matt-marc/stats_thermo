#!/usr/bin/gnuplot -persist
unset key
set xlabel "Distance (x)"
set ylabel "Desnity"
set title "Density of paricles for given x"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Inital density of paricles"
plot "test.dat" u 1:2
