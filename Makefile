OBJS = main.o plot3d.o input.o boundaryConditions.o eos.o primVars.o inviscidFlux.o procBlock.o viscousFlux.o output.o matrix.o parallel.o slices.o turbulence.o
CC = mpic++
DEBUG = -ggdb -pg 
OPTIM = -O3 -march=native
PROF = -O3 -march=native -pg
CODENAME = main
CFLAGS = -std=c++11 -Wall -c $(OPTIM)
LFLAGS = -std=c++11 -Wall $(OPTIM) -o $(CODENAME)

$(CODENAME) : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

plot3d.o : plot3d.cpp plot3d.h vector3d.h
	$(CC) $(CFLAGS) plot3d.cpp

main.o : main.cpp plot3d.h vector3d.h input.h procBlock.h eos.h primVars.h boundaryConditions.h inviscidFlux.h tensor.h viscousFlux.h output.h matrix.h parallel.h turbulence.h
	$(CC) $(CFLAGS) main.cpp

input.o : input.cpp input.h boundaryConditions.h
	$(CC) $(CFLAGS) input.cpp

primVars.o : primVars.cpp primVars.h vector3d.h eos.h inviscidFlux.h boundaryConditions.h input.h
	$(CC) $(CFLAGS) primVars.cpp

procBlock.o : procBlock.cpp procBlock.h vector3d.h plot3d.h eos.h primVars.h inviscidFlux.h input.h matrix.h viscousFlux.h boundaryConditions.h macros.h turbulence.h
	$(CC) $(CFLAGS) procBlock.cpp

inviscidFlux.o : inviscidFlux.cpp vector3d.h eos.h primVars.h inviscidFlux.h input.h
	$(CC) $(CFLAGS) inviscidFlux.cpp

boundaryConditions.o : boundaryConditions.cpp boundaryConditions.h plot3d.h vector3d.h
	$(CC) $(CFLAGS) boundaryConditions.cpp

eos.o : eos.cpp eos.h vector3d.h
	$(CC) $(CFLAGS) eos.cpp

slices.o : slices.cpp procBlock.h
	$(CC) $(CFLAGS) slices.cpp

viscousFlux.o : viscousFlux.cpp vector3d.h tensor.h eos.h primVars.h viscousFlux.h input.h
	$(CC) $(CFLAGS) viscousFlux.cpp

output.o : output.cpp output.h procBlock.h tensor.h vector3d.h plot3d.h eos.h primVars.h inviscidFlux.h input.h turbulence.h
	$(CC) $(CFLAGS) output.cpp

parallel.o : parallel.cpp parallel.h primVars.h procBlock.h vector3d.h plot3d.h boundaryConditions.h matrix.h
	$(CC) $(CFLAGS) parallel.cpp

matrix.o : matrix.cpp matrix.h plot3d.h macros.h
	$(CC) $(CFLAGS) matrix.cpp

turbulence.o : turbulence.cpp turbulence.h
	$(CC) $(CFLAGS) turbulence.cpp

clean:
	rm *.o *~ $(CODENAME)
