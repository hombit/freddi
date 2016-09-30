CC = g++
CPP = g++
CPPFLAGS = -std=c++11

LDLIBS = -lboost_program_options

OBJ = nonlinear_diffusion.o opacity_related.o orbit.o spectrum.o


all: freddi
fred: $(OBJ) freddi.o

clean:
	rm -f *.o
