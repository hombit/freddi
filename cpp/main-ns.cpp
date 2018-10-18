#include "application.hpp"
#include "freddi.hpp"


int main(int ac, char *av[]) {
	run_application<FreddiNeutronStarOptions, FreddiNeutronStarEvolution>(ac, av);
	return 0;
}
