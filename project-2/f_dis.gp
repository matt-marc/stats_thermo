#!/usr/bin/gnuplot

set terminal pngcairo  enhanced font "arial,14" fontscale 1.0 size 900, 600

unset key
set xlabel "Velocity X (m/s)"
set ylabel "Number distribution"
set title "Density of paricles for given x velocity"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Velocity space number distribution at x = 8.025m after 6ms"
set output "right_f.png"
plot "right_final.dat" u 1:2

set title "Velocity space number distribution at x = 3.025m after 6ms"
set output "left_f.png"
plot "left_final.dat" u 1:2

set title "Velocity space number distribution at x = 4.525m after 6ms"
set output "center_shock.png"
plot "shock.dat" u 1:2

set title "Initial Velocity space number distribution at x = 8.425m"
set output "right_init_conditions.png"
plot "right_init.dat" u 1:2

set title "Initial Velocity space number distribution at x = 0.525m"
set output "left_init_conditions.png"
plot "left_init.dat" u 1:2