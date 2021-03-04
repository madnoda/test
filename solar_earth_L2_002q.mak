all:	solar_earth_L2_002q_01XYd.eps

solar_earth_L2_002q_01:	solar_earth_L2_002q.cpp vector3q.cpp vector3q.h
	g++ -o solar_earth_L2_002q_01 solar_earth_L2_002q.cpp vector3q.cpp -lm -static-libstdc++ -lquadmath

solar_earth_L2_002q_01.dat:	solar_earth_L2_002q_01
	./solar_earth_L2_002q_01 300000000.0000 solar_earth_L2_002q_01.dat

solar_earth_L2_002q_01XYd.eps:	solar_earth_L2_002q_01XYd.plt solar_earth_L2_002q_01.dat
	gnuplot solar_earth_L2_002q_01XYd.plt
