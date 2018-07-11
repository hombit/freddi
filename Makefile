CC = g++
CPP = g++
prefix=/usr/local
override CPPFLAGS += -std=c++11 -O2 -D INSTALLPATHPREFIX='"$(prefix)"'

override LDLIBS += -lboost_program_options

OBJ = arguments.o nonlinear_diffusion.o opacity_related.o orbit.o spectrum.o


all: freddi ini
freddi: $(OBJ) freddi.o
ini: freddi.ini
	sed -e '3 s_[a-zA-Z/]*/etc/_$(prefix)/etc/_' freddi.ini > .freddi.ini
	mv .freddi.ini freddi.ini

readme: all
	./freddi --help > ./.freddi_help_message
	sed -e '/\.\/freddi --help/,/~~~/ {//!d;}' -e '/\.\/freddi --help/r .freddi_help_message' Readme.md > .freddi_Readme.md
	mv .freddi_Readme.md Readme.md
	rm -f ./.freddi_help_message

install: all
	install -m 0755 freddi $(prefix)/bin
	install -m 0755 freddi.ini $(prefix)/etc

test: all
	./freddi

clean:
	rm -f *.o
