#pragma once
#include <array>
#include <utility>
#include <vector>

#include "Bond.hpp"
#include "SimulationInput.hpp"
#include "SmallPRNG/prng.h"

template <class IntegerType, class FloatType> struct Lattice {
  Lattice() { initialized = false; }

  Lattice(SimulationInput<IntegerType, FloatType> in) {

    initialized = false;
    sim_input = in;
  }

  template <class PrngType> void initialize(PrngType &prng) {

    auto nx = sim_input.nx;
    auto ny = sim_input.ny;

    spin_mask.resize(nx * ny, false);

    std::vector<IntegerType> masked_spin_indices(nx * ny, -1);
    active_sites = 0;

    IntegerType bond_counter = 0;

    for (auto i = 0; i < nx * ny; i++) {
      auto is_site_occupied = prng.randf() > sim_input.dilution;
      spin_mask[i] = is_site_occupied;
      if (is_site_occupied) {
        masked_spin_indices[i] = active_sites;
        std::cout << i << " " << active_sites << std::endl;
        active_sites += 1;
      }
    }

    IntegerType j = 0;

    while (j < ny) {
      IntegerType i = 0;

      while (i < nx) {

        auto current_spin = (i % nx) + j * nx;
        auto next_spin_x = (i + 1) % nx + j * nx;
        auto next_spin_y = (current_spin + nx) % sim_input.max_sites;

        // If the current spin doesn't exist, move on.
        if (!spin_mask[current_spin]) {
          i += 1;
          continue;
        }

        // If the next spin horizontally exists, add it to the bond,
        // increase bond counter.
        if (spin_mask[next_spin_x]) {
          bondsites.push_back(
              Bond<IntegerType>{masked_spin_indices[current_spin],
                                masked_spin_indices[next_spin_x]});
          bond_counter += 1;
        }

        if (spin_mask[next_spin_y]) {
          bondsites.push_back(
              Bond<IntegerType>{masked_spin_indices[current_spin],
                                masked_spin_indices[next_spin_y]});
          bond_counter += 1;
        }

        i += 1;
      }

      j += 1;
    }

    nbonds = bond_counter;

    for (auto i = 0; i < active_sites; i++) {
      if (prng.randf() < 0.5) {
        spins.push_back(1);
      } else {
        spins.push_back(-1);
      }
    }

    initialized = true;
  }

  bool initialized;
  IntegerType active_sites;
  IntegerType nbonds;
  SimulationInput<IntegerType, FloatType> sim_input;

  std::vector<bool> spin_mask;
  std::vector<int8_t> spins;
  std::vector<Bond<IntegerType>> bondsites;
};