#ifndef FREDDI_APPLICATION_H
#define FREDDI_APPLICATION_H

#include <iostream>

#include <boost/program_options.hpp>

#include "options.hpp"
#include "output.hpp"

namespace po = boost::program_options;


template <typename Output, typename Options, typename Evolution>
bool run_application(int ac, char *av[]) {
	po::variables_map vm;
	if (! parseOptions<Options>(vm, ac, av)){
		return false;
	}
	Options opts(vm);
	std::shared_ptr<Evolution> freddi{new Evolution(opts)};
	Output output(freddi, vm);
	for (int i_t = 0; i_t <= static_cast<int>(freddi->args().calc->time / freddi->args().calc->tau); i_t++) {
		if (i_t % freddi->args().general->temp_sparsity_output == 0) {
			output.dump();
		}
		freddi->step();
	}
	return true;
}

#endif //FREDDI_APPLICATION_H
