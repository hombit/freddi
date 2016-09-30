CC = g++
CPP = g++
CPPFLAGS = -std=c++11
prefix=/usr/local

LDLIBS = -lboost_program_options

OBJ = nonlinear_diffusion.o opacity_related.o orbit.o spectrum.o


all: freddi
freddi: $(OBJ) freddi.o

install: all
    install -m 0755 freddi $(prefix)/bin

clean:
	rm -f *.o
