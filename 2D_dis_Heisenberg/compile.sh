mpif90 -cpp -c vars_mods.f90
mpif90 -c map.f90
mpif90 -c misc.f90
mpif90 -c meas.f90
mpif90 -c qmc.f90
mpif90 -c vtx.f90
mpif90 -c upd.f90
mpif90 -c main.f90

mpif90 -llapack -o sse vars_mods.o misc.o vtx.o upd.o main.o map.o qmc.o meas.o

rm *.o *.mod
#ifort -c xxz2d.f90
