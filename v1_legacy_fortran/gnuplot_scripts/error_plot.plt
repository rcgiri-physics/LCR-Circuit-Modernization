#set xzeroaxis linetype -1
set xtics auto

set xrange [0:2]
#set autoscale

set xl 'Time'  
set yl 'Absolute Error (A)'  
 
plot 'data_NR_A.dat' u 2:6 w l lw 0.5 lc 'blue' title 'Exact Current - Calculated Current'