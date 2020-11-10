#include "SimulationInput.hpp"

#include <iostream>
#include <string>
#include <sstream>

auto initialize_from_console(int arg_count, char *arg_array[]) {
  uint16_t lx = std::stoi(arg_array[1]);
  uint16_t ly = std::stoi(arg_array[2]);
  uint16_t lz = std::stoi(arg_array[3]);
  float beta = std::stof(arg_array[4]);
  float dilution = std::stof(arg_array[5]);
  unsigned int nbins = std::stof(arg_array[6]);
  unsigned int msteps = std::stof(arg_array[7]);
  unsigned int isteps = std::stof(arg_array[8]);

  SimulationInput<int, float> sim_input{lx, ly, lz, beta, dilution};
  sim_input.nbins = nbins;
  sim_input.msteps = msteps;
  sim_input.isteps = isteps;

  if (arg_count > 9) {
    sim_input.resultsfile = arg_array[9];
  }

  return sim_input;
}