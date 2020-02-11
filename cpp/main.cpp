#include <application.hpp>
#include <freddi_evolution.hpp>
#include <options.hpp>
#include <output.hpp>


int main(int ac, char *av[]) {
	if (! run_application<FreddiFileOutput, FreddiOptions, FreddiEvolution>(ac, av)){
		return 1;
	}
	return 0;
}
