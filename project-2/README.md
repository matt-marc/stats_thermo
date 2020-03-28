# Thermo project

## Structure

scr directory contains all code

- DVS.h and cpp contain all the time marching code
- q1.cpp and q2.cpp sets initial conditions for question 1 and 2

all .gp files were used to generate graphs found in the report


## Running

makefile in main directory will build project located in the
src directory this will copy over q1 and q2 that can be run in the main
directory

Running q1 and q2 will generate inital conditiion .dat files and a final .dat
file ms6.dat that contain the properties along the x dir after 6ms