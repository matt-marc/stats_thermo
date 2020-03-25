#!/usr/bin/gnuplot

set terminal pngcairo  enhanced font "arial,14" fontscale 1.0 size 900, 600

unset key
set xlabel "Distance x [m]"

set title "Density of paricles for given x velocity"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Initial Density of Particles"
set ylabel "Density"
set output "init-rho.png"
plot "initial_con.dat" u 1:2

set title "Initial Velocity of Particles"
set ylabel "Velocity"
set output "init-u.png"
plot "initial_con.dat" u 1:3

set title "Initial Pressure of Particles"
set ylabel "Pressure"
set output "init-p.png"
plot "initial_con.dat" u 1:4

set title "Initial Heat Flux of Particles"
set ylabel "Heat Flux"
set output "init-q.png"
plot "initial_con.dat" u 1:5
