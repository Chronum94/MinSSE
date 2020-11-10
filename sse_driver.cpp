// #include "readinput.hpp"
#include "Lattice.hpp"
#include "Simulation.hpp"
#include "SimulationInput.hpp"
#include "readinput.hpp"

#include <iostream>
#include <fstream>

#include "SmallPRNG/prng.h"


#if XORSHIFT128
using PRNGTYPE = smallprng::xor128;
constexpr auto PRNG_DESC = "XORSHIFT128";
#elif XORSHIFT32
using PRNGTYPE = smallprng::xor32;
constexpr auto PRNG_DESC = "XORSHIFT32";
#elif SQUARES
using PRNGTYPE = smallprng::improved_squares;
constexpr auto PRNG_DESC = "IMPROVED_SQUARES";
#elif AES4
using PRNGTYPE = smallprng::aes_4;
constexpr auto  PRNG_DESC = "AES(4)";
#elif LCG
using PRNGTYPE = smallprng::knuth_lcg;
constexpr auto  PRNG_DESC = "KNUTH_LCG";
#else
using PRNGTYPE = smallprng::improved_squares;
constexpr auto  PRNG_DESC = "IMPROVED_SQUARES";
#endif



int main(int argc, char *argv[]) {
  if (argc > 1) {
    // std::cout << argc << std::endl;
    auto sim_input = initialize_from_console(argc, argv);

    std::ofstream outfile("out.log");
    outfile << sim_input << std::endl;
    outfile << "PRNG: " << PRNG_DESC << "\n";
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