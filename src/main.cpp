#include <iostream>

#include "arguments.hpp"
#include "freddi.hpp"
#include "output.hpp"

using namespace std;


int main(int ac, char *av[]) {
	auto vm = parseArguments(ac, av);
	if (vm.count("help") > 0) {
		cout << FreddiArguments::description() << endl;
		return 0;
	}
	FreddiArguments args(vm);
	Freddi freddi(args);
	FreddiFileOutput output(freddi, vm);
	for (int i_t = 0; i_t <= static_cast<int>(freddi.args->calc->time / freddi.args->calc->tau); i_t++) {
		output.dump();
		freddi.step();
	}
}
