CC = g++
CPP = g++
CPPFLAGS = -std=c++11

LDLIBS = -lgsl -lgslcblas -lboost_program_options

OBJ = nonlinear_diffusion.o spectrum.o

all: fred
fred: $(OBJ) fred.o
