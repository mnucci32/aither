OBJS = main.o plot3d.o input.o boundaryConditions.o eos.o primVars.o procBlock.o output.o matrix.o parallel.o slices.o turbulence.o inviscidFlux.o viscousFlux.o source.o resid.o kdtree.o genArray.o fluxJacobian.o uncoupledScalar.o
CC = mpic++
DEBUG = -ggdb -pg
OPTIM = -O3 -march=native
PROF = -O3 -march=native -pg
CODENAME = aither
CFLAGS = -std=c++14 -Wall -c $(OPTIM)
LFLAGS = -std=c++14 -Wall $(OPTIM) -o $(CODENAME)

$(CODENAME) : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

plot3d.o : plot3d.cpp plot3d.hpp vector3d.hpp multiArray3d.hpp
	$(CC) $(CFLAGS) plot3d.cpp

main.o : main.cpp plot3d.hpp vector3d.hpp input.hpp procBlock.hpp eos.hpp primVars.hpp boundaryConditions.hpp inviscidFlux.hpp tensor.hpp viscousFlux.hpp output.hpp parallel.hpp turbulence.hpp resid.hpp multiArray3d.hpp genArray.hpp fluxJacobian.hpp
	$(CC) $(CFLAGS) main.cpp

input.o : input.cpp input.hpp boundaryConditions.hpp
	$(CC) $(CFLAGS) input.cpp

primVars.o : primVars.cpp primVars.hpp vector3d.hpp eos.hpp inviscidFlux.hpp boundaryConditions.hpp input.hpp macros.hpp genArray.hpp
	$(CC) $(CFLAGS) primVars.cpp

procBlock.o : procBlock.cpp procBlock.hpp vector3d.hpp plot3d.hpp eos.hpp primVars.hpp inviscidFlux.hpp input.hpp genArray.hpp viscousFlux.hpp boundaryConditions.hpp macros.hpp turbulence.hpp kdtree.hpp uncoupledScalar.hpp fluxJacobian.hpp matrix.hpp
	$(CC) $(CFLAGS) procBlock.cpp

inviscidFlux.o : inviscidFlux.cpp vector3d.hpp eos.hpp primVars.hpp inviscidFlux.hpp input.hpp macros.hpp genArray.hpp turbulence.hpp matrix.hpp
	$(CC) $(CFLAGS) inviscidFlux.cpp

boundaryConditions.o : boundaryConditions.cpp boundaryConditions.hpp plot3d.hpp vector3d.hpp
	$(CC) $(CFLAGS) boundaryConditions.cpp

eos.o : eos.cpp eos.hpp vector3d.hpp
	$(CC) $(CFLAGS) eos.cpp

slices.o : slices.cpp slices.hpp procBlock.hpp
	$(CC) $(CFLAGS) slices.cpp

viscousFlux.o : viscousFlux.cpp vector3d.hpp tensor.hpp eos.hpp primVars.hpp viscousFlux.hpp input.hpp turbulence.hpp macros.hpp
	$(CC) $(CFLAGS) viscousFlux.cpp

output.o : output.cpp output.hpp procBlock.hpp tensor.hpp vector3d.hpp plot3d.hpp eos.hpp primVars.hpp inviscidFlux.hpp input.hpp turbulence.hpp genArray.hpp
	$(CC) $(CFLAGS) output.cpp

parallel.o : parallel.cpp parallel.hpp primVars.hpp procBlock.hpp vector3d.hpp plot3d.hpp boundaryConditions.hpp resid.hpp
	$(CC) $(CFLAGS) parallel.cpp

matrix.o : matrix.cpp matrix.hpp macros.hpp genArray.hpp
	$(CC) $(CFLAGS) matrix.cpp

genArray.o : genArray.cpp genArray.hpp macros.hpp
	$(CC) $(CFLAGS) genArray.cpp

turbulence.o : turbulence.cpp turbulence.hpp
	$(CC) $(CFLAGS) turbulence.cpp

source.o : source.cpp source.hpp macros.hpp turbulence.hpp primVars.hpp
	$(CC) $(CFLAGS) source.cpp

resid.o : resid.cpp resid.hpp
	$(CC) $(CFLAGS) resid.cpp

kdtree.o : kdtree.cpp kdtree.hpp vector3d.hpp
	$(CC) $(CFLAGS) kdtree.cpp

fluxJacobian.o : fluxJacobian.cpp fluxJacobian.hpp turbulence.hpp vector3d.hpp primVars.hpp eos.hpp input.hpp genArray.hpp matrix.hpp inviscidFlux.hpp uncoupledScalar.hpp
	$(CC) $(CFLAGS) fluxJacobian.cpp

uncoupledScalar.o : uncoupledScalar.cpp uncoupledScalar.hpp genArray.hpp
	$(CC) $(CFLAGS) uncoupledScalar.cpp

clean:
	rm *.o *~ $(CODENAME)
