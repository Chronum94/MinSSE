// #include "readinput.hpp"
#include "Lattice.hpp"
#include "Simulation.hpp"
#include "SimulationInput.hpp"
#include "readinput.hpp"

#include <iostream>
#include <fstream>

#include "SmallPRNG/prng.h"


#ifdef XORSHIFT128
#define PRNGTYPE smallprng::xor128
#elif XORSHIFT32
#define PRNGTYPE smallprng::xor32
#elif SQUARES
#define PRNGTYPE smallprng::improved_squares
#elif AES4
#define PRNGTYPE smallprng::aes_4
#elif LCG
#define PRNGTYPE smallprng::knuth_lcg
#else
#define PRNGTYPE smallprng::improved_squares
#endif



int main(int argc, char *argv[]) {
  if (argc > 1) {
    // std::cout << argc << std::endl;
    auto sim_input = initialize_from_console(argc, argv);

    std::ofstream outfile("out.log");
    outfile << sim_input << std::endl;
    Lattice<int, float> sim_lattice(sim_input);

    // prng<4, uint32_t, xorshift128> myprng;
    // prng<1, uint32_t, xorshift32> myprng;
    // prng<4, uint32_t, squares> myprng;
    PRNGTYPE myprng;

    sim_lattice.initialize(myprng);
    // std::cout << sim_lattice.n_active_sites << std::endl;

    Simulation mysim(sim_input, sim_lattice, myprng);
    mysim.run(outfile);
    outfile.close();
  }
  else {
    return 0;
    // auto sim_input = read_default_file_input();
  }
}