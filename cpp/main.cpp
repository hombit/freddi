#include <iostream>

#include "options.hpp"
#include "freddi.hpp"
#include "output.hpp"

using namespace std;


int main(int ac, char *av[]) {
	auto vm = parseOptions<FreddiOptions>(ac, av);
	if (vm.count("help") > 0) {
		cout << FreddiOptions::description() << endl;
		return 0;
	}
	FreddiOptions args(vm);
	FreddiEvolution freddi(args);
	FreddiFileOutput output(freddi, vm);
	for (int i_t = 0; i_t <= static_cast<int>(freddi.args->calc->time / freddi.args->calc->tau); i_t++) {
		output.dump();
		freddi.step();
	}
}
