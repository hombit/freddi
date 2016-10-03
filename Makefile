CC = g++
CPP = g++
CPPFLAGS = -std=c++11
prefix=/usr/local

LDLIBS = -lboost_program_options

OBJ = nonlinear_diffusion.o opacity_related.o orbit.o spectrum.o


all: freddi
freddi: $(OBJ) freddi.o

readme: all
	./freddi --help > ./help_message
	sed -e '/freddi --help/r help_message' Readme_template.md > Readme.md
	rm -f ./help_message

install: all
	install -m 0755 freddi $(prefix)/bin

clean:
	rm -f *.o
