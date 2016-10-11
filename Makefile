CC = g++
CPP = g++
CPPFLAGS = -std=c++11
prefix=/usr/local

LDLIBS = -lboost_program_options

OBJ = nonlinear_diffusion.o opacity_related.o orbit.o spectrum.o


all: freddi
freddi: $(OBJ) freddi.o

readme: all
	./freddi --help > ./.freddi_help_message
	sed -e '/\.\/freddi --help/,/~~~/ {//!d;}' -e '/\.\/freddi --help/r .freddi_help_message' Readme.md > .freddi_Readme.md
	mv .freddi_Readme.md Readme.md
	rm -f ./.freddi_help_message

install: all freddi.ini
	install -m 0755 freddi $(prefix)/bin
	install -m 0755 freddi.ini $(prefix)/etc

clean:
	rm -f *.o
