// #include "readinput.hpp"
#include "Lattice.hpp"
#include "Simulation.hpp"
#include "SimulationInput.hpp"
#include <iostream>

#include "SmallPRNG/prng.h"

int main(int argc, char *argv[]) {
  SimulationInput<int, float> sim_input{4, 4, 2, 0.0};
  sim_input.nbins = 3;
  sim_input.msteps = 20;
  sim_input.isteps = 20;
  sim_input.seed = 14746845679U;

  std::cout << sim_input << std::endl;
  // Lattice<uint32_t, float> sim_lattice{sim_input};
  Lattice<int, float> sim_lattice(sim_input);
  Lattice<int, float> trial_lattice;
  std::cout << "Initialized:" << sim_lattice.initialized << std::endl;
  std::cout << "Initialized:" << trial_lattice.initialized << std::endl;

  prng<4, uint32_t, xorshift128> myprng;
  // myprng.set_state(34);
  sim_lattice.initialize(myprng);
//   for (auto e : sim_lattice.spin_mask) {
//     std::cout << "Mask: " << e << "\n";
//   }

  std::cout << "Bonds:\n";
  for (auto bond : sim_lattice.bondsites) {
    std::cout << bond.s1 << " " << bond.s2 << "\t";
    std::cout << static_cast<int>(sim_lattice.spins[bond.s1]) << " "
              << static_cast<int>(sim_lattice.spins[bond.s2]) << "\n";
  }
  // prng_state<4> s;
  // auto attempt = xorshift128(s);

  Simulation mysim(sim_input, sim_lattice, myprng);
  mysim.run();
}