FORT=gfortran-4.8


# this makes someting
lib_test_roots.dylib : constants.o testfunctions.o rootfinding.o test_roots.o
	$(FORT) -shared test_roots.o rootfinding.o testfunctions.o constants.o -o lib_test_roots.dylib

#test_roots.out : constants.o testfunctions.o rootfinding.o test_roots.o
#	$(FORT) -o test_roots.out test_roots.o rootfinding.o testfunctions.o constants.o

test_roots.o : test_roots.f90
	$(FORT) -c test_roots.f90

rootfinding.o : rootfinding.f90
	$(FORT) -c rootfinding.f90

testfunctions.o : testfunctions.f90
	$(FORT) -c testfunctions.f90

constants.o : constants.f90
	$(FORT) -c constants.f90
