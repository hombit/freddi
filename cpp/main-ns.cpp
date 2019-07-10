#include "application.hpp"
#include "freddi_evolution.hpp"
#include "output.hpp"


int main(int ac, char *av[]) {
	run_application<FreddiNeutronStarOptions, FreddiNeutronStarFileOutput, FreddiNeutronStarEvolution>(ac, av);
	return 0;
}
