#!/usr/bin/gnuplot

set terminal pngcairo  enhanced font "arial,14" fontscale 1.0 size 900, 600

unset key
set xlabel "Velocity X (m/s)"
set ylabel "Number distribution"
set title "Density of paricles for given x velocity"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Inital density of paricles for given x velocity right"
set output "right_init_conditions.png"
plot "right.dat" u 1:2

set title "Inital density of paricles for given x velocity left"
set output "left_init_conditions.png"
plot "left.dat" u 1:2
