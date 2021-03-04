all:	solar_earth_L2_002_01XYd.eps

solar_earth_L2_002_01:	solar_earth_L2_002.cpp vector3.cpp vector3.h
	g++ -o solar_earth_L2_002_01 solar_earth_L2_002.cpp vector3.cpp -lm -static-libstdc++

solar_earth_L2_002_01.dat:	solar_earth_L2_002_01
	./solar_earth_L2_002_01 200000000.0000 solar_earth_L2_002_01.dat

solar_earth_L2_002_01XYd.eps:	solar_earth_L2_002_01XYd.plt solar_earth_L2_002_01.dat
	gnuplot solar_earth_L2_002_01XYd.plt

