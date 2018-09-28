#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <memory>

#include "common.hpp"
#include "params.hpp"
#include "simulation.hpp"
#include "file.hpp"

// arg: init.prof out_dir
int main(int argc, char **argv){
	namespace fs = std::filesystem;
	simulation_t sim;

	if(argc != 3) return -1;

#ifdef OPENMP
	std::cout << "[OpenMP] num_procs: " << omp_get_num_procs() << std::endl;
#endif

	// output directoryの準備
#ifndef BENCH
	out_dir = argv[2];
	if(!check_outdir(out_dir)) return -1;
#endif

	load_data(argv[1], sim);

	set_param(sim);

	sim_loop(sim);

	return 0;
}


