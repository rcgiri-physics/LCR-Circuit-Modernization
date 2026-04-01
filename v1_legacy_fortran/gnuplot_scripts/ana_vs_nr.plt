set xzeroaxis linetype -1
#set xtics auto

set xrange [0:2]
#set autoscale

set xl 'Time'  
set yl 'Current (A)'  

plot "data_NR_A.dat" u 2:5 w l lw 0.5 lc "blue" title 'Analytic Current',\
"data_NR_A.dat" u 2:3 w p lw 0.25 lc 'red' title 'Numerical method i(t)'