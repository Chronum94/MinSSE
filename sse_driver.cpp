// #include "readinput.hpp"
#include "Lattice.hpp"
#include "Simulation.hpp"
#include "SimulationInput.hpp"
#include <iostream>

#include "readinput.hpp"

#include "SmallPRNG/prng.h"

int main(int argc, char *argv[]) {
  if (argc > 1) {
    std::cout << argc << std::endl;
    auto sim_input = initialize_from_console(argc, argv);
    std::cout << sim_input << std::endl;
    Lattice<int, float> sim_lattice(sim_input);

    prng<4, uint32_t, xorshift128> myprng;

    sim_lattice.initialize(myprng);

    Simulation mysim(sim_input, sim_lattice, myprng);
    mysim.run();
  }
  else {
    return 0;
    // auto sim_input = read_default_file_input();
  }
}