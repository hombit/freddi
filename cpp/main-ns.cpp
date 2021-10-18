#include <application.hpp>
#include <ns/ns_evolution.hpp>
#include <ns/ns_output.hpp>
#include <ns/ns_options.hpp>


int main(int ac, char *av[]) {
	if (! run_application< FreddiNeutronStarFileOutput, FreddiNeutronStarOptions, FreddiNeutronStarEvolution>(ac, av)){
		return 1;
	}
	return 0;
}
