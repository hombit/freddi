#include "application.hpp"
#include "freddi_evolution.hpp"


int main(int ac, char *av[]) {
	run_application<FreddiOptions, FreddiEvolution>(ac, av);
	return 0;
}
