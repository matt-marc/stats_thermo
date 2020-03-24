#!/usr/bin/gnuplot -persist
unset key
set xlabel "Distance (x)"
set ylabel "Density"
set title "Density of paricles for given x"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Inital density of paricles"
plot "test.dat" u 1:2

set title "Inital velocity of paricles"
plot "test.dat" u 1:3

set title "Inital pressure of paricles"
plot "test.dat" u 1:4
