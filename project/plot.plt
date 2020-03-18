set terminal pngcairo  enhanced font "arial,14" fontscale 1.0 size 900, 600

unset key
#unset colorbox
set xlabel "position"
set ylabel "U"
set title "Average vel over domain"

set yrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback



set output "domain.png"
plot "ms6.dat" u 1:3
