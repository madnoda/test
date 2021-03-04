set terminal postscript eps enhanced color solid 14
set output "solar_earth_L2_002_01XYd.eps"
set grid
set ytics nomirror
plot"solar_earth_L2_002_01.dat" using 8:9 title "1" with line lw 5,"solar_earth_L2_002_01.dat" using 14:15 title "2" with line lw 5,"solar_earth_L2_002_01.dat" using 20:21 title "3" with line lw 5,"solar_earth_L2_002_01.dat" using 26:27 title "4" with line lw 5,"solar_earth_L2_002_01.dat" using 32:33 title "5" with line lw 5,"solar_earth_L2_002_01.dat" using 38:39 title "6" with line lw 5
