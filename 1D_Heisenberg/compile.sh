gfortran -c xxz2d.f90
gfortran -c upd.f90
gfortran -c vtx.f90
gfortran -o sse xxz2d.o upd.o vtx.o

#ifort -llapack -g -check all -fpe0 -warn -traceback -debug extended -c xxz2d.f90
#ifort -llapack -g -check all -fpe0 -warn -traceback -debug extended -c upd.f90
#ifort -llapack -g -check all -fpe0 -warn -traceback -debug extended -c vtx.f90
#ifort -llapack -g -check all -fpe0 -warn -traceback -debug extended -o sse xxz2d.o upd.o vtx.o

rm *genmod*
rm *.o 
rm *.mod
