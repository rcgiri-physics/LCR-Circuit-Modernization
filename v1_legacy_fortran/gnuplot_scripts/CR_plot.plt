set xzeroaxis linetype -1
#set xtics auto

#set yrange [-0.001:0.004]
#set xrange [0:2]
set autoscale

set xl 'Time'  
set yl 'Current (A)'  

plot "data_complete_ud.dat" u 2:7 w l lw 1 lc "blue" title 'R = 10 Ohm',\
"data_complete_cd.dat" u 2:7 w l lw 1 lc 'red' title 'R = 40 Ohm',\
"data_complete_od.dat" u 2:7 w l lw 1 lc 'green' title 'R = 60 Ohm'