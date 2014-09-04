OBJS = main.o plot3d.o input.o boundaryConditions.o eos.o primVars.o inviscidFlux.o blockVars.o viscousFlux.o viscBlockVars.o output.o matrix.o
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

main.o : main.cpp plot3d.h vector3d.h input.h blockVars.h eos.h primVars.h boundaryConditions.h inviscidFlux.h tensor.h viscousFlux.h viscBlockVars.h output.h matrix.h
	$(CC) $(CFLAGS) main.cpp

input.o : input.cpp input.h boundaryConditions.h
	$(CC) $(CFLAGS) input.cpp

primVars.o : primVars.cpp primVars.h vector3d.h eos.h inviscidFlux.h boundaryConditions.h input.h
	$(CC) $(CFLAGS) primVars.cpp

blockVars.o : blockVars.cpp blockVars.h vector3d.h plot3d.h eos.h primVars.h inviscidFlux.h input.h matrix.h
	$(CC) $(CFLAGS) blockVars.cpp

inviscidFlux.o : inviscidFlux.cpp vector3d.h eos.h primVars.h inviscidFlux.h boundaryConditions.h input.h
	$(CC) $(CFLAGS) inviscidFlux.cpp

boundaryConditions.o : boundaryConditions.cpp boundaryConditions.h
	$(CC) $(CFLAGS) boundaryConditions.cpp

eos.o : eos.cpp eos.h vector3d.h
	$(CC) $(CFLAGS) eos.cpp

viscousFlux.o : viscousFlux.cpp vector3d.h tensor.h eos.h primVars.h viscousFlux.h boundaryConditions.h input.h
	$(CC) $(CFLAGS) viscousFlux.cpp

viscBlockVars.o : viscBlockVars.cpp viscBlockVars.h vector3d.h tensor.h plot3d.h eos.h input.h viscousFlux.h blockVars.h primVars.h
	$(CC) $(CFLAGS) viscBlockVars.cpp

output.o : output.cpp output.h blockVars.h tensor.h vector3d.h plot3d.h eos.h primVars.h inviscidFlux.h input.h viscBlockVars.h
	$(CC) $(CFLAGS) output.cpp

matrix.o : matrix.cpp matrix.h plot3d.h
	$(CC) $(CFLAGS) matrix.cpp

testMain.o : testMain.cpp matrix.h
	$(CC) $(CFLAGS) testMain.cpp

testMain : testMain.o matrix.o
	$(CC) $(LFLAGS) testMain.o matrix.o

clean:
	rm *.o *~ $(CODENAME)
