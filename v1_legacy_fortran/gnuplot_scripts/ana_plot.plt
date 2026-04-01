#set title "Graph of Natural response" 

set xzeroaxis linetype -1
set xtics auto

set xrange [0:2]
#set autoscale

set xl 'Time'  
set yl 'Current (A)'  
 
plot 'data_NR_A.dat' u 2:5 w l lw 1 lc 'blue' title 'i(t)'

