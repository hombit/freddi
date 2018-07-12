#include <iostream>

#include "arguments.hpp"
#include "freddi.hpp"

using namespace std;


int main(int ac, char *av[]) {
	auto vm = parseArguments(ac, av);
	if (vm.count("help") > 0) {
		cout << FreddiArguments::description() << endl;
		return 0;
	}
	FreddiArguments args(vm);
	freddi(args);
}
