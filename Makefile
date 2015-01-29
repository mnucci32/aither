OBJS = main.o plot3d.o input.o boundaryConditions.o eos.o primVars.o inviscidFlux.o procBlock.o viscousFlux.o output.o matrix.o geomSlice.o
CC = g++
DEBUG = -ggdb -pg 
OPTIM = -O3 -march=native
CODENAME = main
CFLAGS = -std=c++0x -Wall -c $(OPTIM)
LFLAGS = -std=c++0x -Wall $(OPTIM) -o $(CODENAME)

$(CODENAME) : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

plot3d.o : plot3d.cpp plot3d.h vector3d.h
	$(CC) $(CFLAGS) plot3d.cpp

main.o : main.cpp plot3d.h vector3d.h input.h procBlock.h eos.h primVars.h boundaryConditions.h inviscidFlux.h tensor.h viscousFlux.h output.h matrix.h
	$(CC) $(CFLAGS) main.cpp

input.o : input.cpp input.h boundaryConditions.h
	$(CC) $(CFLAGS) input.cpp

primVars.o : primVars.cpp primVars.h vector3d.h eos.h inviscidFlux.h boundaryConditions.h input.h
	$(CC) $(CFLAGS) primVars.cpp

procBlock.o : procBlock.cpp procBlock.h vector3d.h plot3d.h eos.h primVars.h inviscidFlux.h input.h matrix.h viscousFlux.h boundaryConditions.h geomSlice.h
	$(CC) $(CFLAGS) procBlock.cpp

inviscidFlux.o : inviscidFlux.cpp vector3d.h eos.h primVars.h inviscidFlux.h input.h
	$(CC) $(CFLAGS) inviscidFlux.cpp

boundaryConditions.o : boundaryConditions.cpp boundaryConditions.h plot3d.h vector3d.h
	$(CC) $(CFLAGS) boundaryConditions.cpp

eos.o : eos.cpp eos.h vector3d.h
	$(CC) $(CFLAGS) eos.cpp

viscousFlux.o : viscousFlux.cpp vector3d.h tensor.h eos.h primVars.h viscousFlux.h input.h
	$(CC) $(CFLAGS) viscousFlux.cpp

output.o : output.cpp output.h procBlock.h tensor.h vector3d.h plot3d.h eos.h primVars.h inviscidFlux.h input.h
	$(CC) $(CFLAGS) output.cpp

matrix.o : matrix.cpp matrix.h plot3d.h
	$(CC) $(CFLAGS) matrix.cpp

geomSlice.o : geomSlice.cpp geomSlice.h vector3d.h 
	$(CC) $(CFLAGS) geomSlice.cpp

testMain.o : testMain.cpp
	$(CC) $(CFLAGS) testMain.cpp

testMain : testMain.o
	$(CC) $(LFLAGS) testMain.o

clean:
	rm *.o *~ $(CODENAME)
