OBJS = main.o plot3d.o input.o boundaryConditions.o eos.o primVars.o procBlock.o output.o matrix.o parallel.o slices.o turbulence.o gradients.o inviscidFlux.o viscousFlux.o source.o
CC = mpic++
DEBUG = -ggdb -pg 
OPTIM = -O3 -march=native
PROF = -O3 -march=native -pg
CODENAME = main
CFLAGS = -std=c++11 -Wall -c $(OPTIM)
LFLAGS = -std=c++11 -Wall $(OPTIM) -o $(CODENAME)

$(CODENAME) : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

plot3d.o : plot3d.cpp plot3d.hpp vector3d.hpp
	$(CC) $(CFLAGS) plot3d.cpp

main.o : main.cpp plot3d.hpp vector3d.hpp input.hpp procBlock.hpp eos.hpp primVars.hpp boundaryConditions.hpp inviscidFlux.hpp tensor.hpp viscousFlux.hpp output.hpp matrix.hpp parallel.hpp turbulence.hpp gradients.hpp
	$(CC) $(CFLAGS) main.cpp

input.o : input.cpp input.hpp boundaryConditions.hpp
	$(CC) $(CFLAGS) input.cpp

primVars.o : primVars.cpp primVars.hpp vector3d.hpp eos.hpp inviscidFlux.hpp boundaryConditions.hpp input.hpp macros.hpp
	$(CC) $(CFLAGS) primVars.cpp

procBlock.o : procBlock.cpp procBlock.hpp vector3d.hpp plot3d.hpp eos.hpp primVars.hpp inviscidFlux.hpp input.hpp matrix.hpp viscousFlux.hpp boundaryConditions.hpp macros.hpp turbulence.hpp
	$(CC) $(CFLAGS) procBlock.cpp

inviscidFlux.o : inviscidFlux.cpp vector3d.hpp eos.hpp primVars.hpp inviscidFlux.hpp input.hpp macros.hpp matrix.hpp
	$(CC) $(CFLAGS) inviscidFlux.cpp

boundaryConditions.o : boundaryConditions.cpp boundaryConditions.hpp plot3d.hpp vector3d.hpp
	$(CC) $(CFLAGS) boundaryConditions.cpp

eos.o : eos.cpp eos.hpp vector3d.hpp
	$(CC) $(CFLAGS) eos.cpp

slices.o : slices.cpp slices.hpp procBlock.hpp
	$(CC) $(CFLAGS) slices.cpp

viscousFlux.o : viscousFlux.cpp vector3d.hpp tensor.hpp eos.hpp primVars.hpp viscousFlux.hpp input.hpp turbulence.hpp macros.hpp
	$(CC) $(CFLAGS) viscousFlux.cpp

output.o : output.cpp output.hpp procBlock.hpp tensor.hpp vector3d.hpp plot3d.hpp eos.hpp primVars.hpp inviscidFlux.hpp input.hpp turbulence.hpp
	$(CC) $(CFLAGS) output.cpp

parallel.o : parallel.cpp parallel.hpp primVars.hpp procBlock.hpp vector3d.hpp plot3d.hpp boundaryConditions.hpp matrix.hpp
	$(CC) $(CFLAGS) parallel.cpp

matrix.o : matrix.cpp matrix.hpp plot3d.hpp macros.hpp
	$(CC) $(CFLAGS) matrix.cpp

turbulence.o : turbulence.cpp turbulence.hpp
	$(CC) $(CFLAGS) turbulence.cpp

source.o : source.cpp source.hpp macros.hpp turbulence.hpp primVars.hpp gradients.hpp
	$(CC) $(CFLAGS) source.cpp

gradients.o : gradients.cpp gradients.hpp primVars.hpp vector3d.hpp tensor.hpp procBlock.hpp
	$(CC) $(CFLAGS) gradients.cpp

clean:
	rm *.o *~ $(CODENAME)
