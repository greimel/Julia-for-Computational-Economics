# this makes someting
test.out : constants.o testfunctions.o rootfinding.o test.o
	gfortran -o test.out test.o rootfinding.o testfunctions.o constants.o

test.o : test.f90
	gfortran -c test.f90

rootfinding.o : rootfinding.f90
	gfortran -c rootfinding.f90

testfunctions.o : testfunctions.f90
	gfortran -c testfunctions.f90

constants.o : constants.f90
	gfortran -c constants.f90
