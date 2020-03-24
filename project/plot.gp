#!/usr/bin/gnuplot -persist
unset key
set title "Inital conditions of domain"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set title "Inital density of paricles"
plot "test.dat" u 1:2

set title "Inital velocity of paricles"
plot "test.dat" u 1:3

set title "Inital pressure of paricles"
plot "test.dat" u 1:4



set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set title "Inital density of paricles"
#set xlable "Distance x [m]"
#set ylable "Density"
plot "test.dat" u 1:2
# --- GRAPH b
set title "Inital velocity of paricles"
#set xlable "Distance x [m]"
#set ylable "Velocity [m/s]"
plot "test.dat" u 1:3
# --- GRAPH c
set title "Inital pressure of paricles"
#set xlable "Distance x [m]"
#set ylable "Pressure [Pa]"
plot "test.dat" u 1:4
# --- GRAPH d
set title "Inital pressure of paricles"
#set xlable "Distance x [m]"
plot "test.dat" u 1:4
unset multiplot
