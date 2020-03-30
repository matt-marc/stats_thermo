#!/usr/bin/gnuplot

set terminal pngcairo  enhanced font "arial,14" fontscale 1.0 size 900, 600

unset key
set xlabel "Velocity X (m/s)"
set ylabel "Number distribution"
set title "Density of paricles for given x velocity"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Velocity space number distribution at x = 0.105 after 0.6ms"
set output "left_f.png"
plot "left_final.dat" u 1:2

set title "Velocity space number distribution at x = 0.355 after 0.6ms"
set output "center_shock.png"
plot "shock.dat" u 1:2
