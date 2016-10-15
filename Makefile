CC = g++
CPP = g++
prefix=/usr/local
CPPFLAGS = -std=c++11 -D INSTALLPATHPREFIX='"$(prefix)"'

LDLIBS = -lboost_program_options

OBJ = nonlinear_diffusion.o opacity_related.o orbit.o spectrum.o


all: freddi
freddi: $(OBJ) freddi.o
#ini: freddi.ini
#	sed -i '' -e '3 s_[a-zA-Z/]*/etc/_$(prefix)/etc/_' freddi.ini

readme: all
	./freddi --help > ./.freddi_help_message
	sed -e '/\.\/freddi --help/,/~~~/ {//!d;}' -e '/\.\/freddi --help/r .freddi_help_message' Readme.md > .freddi_Readme.md
	mv .freddi_Readme.md Readme.md
	rm -f ./.freddi_help_message

install: all
	install -m 0755 freddi $(prefix)/bin
	install -m 0755 freddi.ini $(prefix)/etc

clean:
	rm -f *.o
