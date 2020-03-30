#!/usr/bin/gnuplot -persist
unset key
set title "Initial conditions of domain"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set title "density of particles"
#set xlable "Distance x [m]"
#set ylable "Density"
plot "ms6.dat" u 1:2 w l
# --- GRAPH b
set title "velocity of particles"
#set xlable "Distance x [m]"
#set ylable "Velocity [m/s]"
plot "ms6.dat" u 1:3 w l
# --- GRAPH c
set title "pressure of particles"
#set xlable "Distance x [m]"
#set ylable "Pressure [Pa]"
plot "ms6.dat" u 1:4 w l
# --- GRAPH d
set title "heat flux"
#set xlable "Distance x [m]"
plot "ms6.dat" u 1:5 w l
unset multiplot
