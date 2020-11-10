#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>
#include <vector>


#include "Bond.hpp"
#include "SimulationInput.hpp"

template <class IntegerType, class FloatType> struct Lattice {
  Lattice() { initialized = false; }

  Lattice(SimulationInput<IntegerType, FloatType> in) {

    initialized = false;
    sim_input = in;
  }

  template <class PrngType> void initialize(PrngType &prng) {

    auto nx = sim_input.nx;
    auto ny = sim_input.ny;
    auto nz = sim_input.nz;

    spin_mask.resize(nx * ny * nz, false);

    std::vector<IntegerType> masked_spin_indices(nx * ny * nz, -1);
    n_active_sites = 0;

    make_lattice(prng);

    IntegerType active_site = 0;
    for (auto i = 0; i < sim_input.max_sites; i++) {
      if (spin_mask[i]) {
        masked_spin_indices[i] = active_site;
        active_site += 1;
      }
    }

    std::ofstream spinmask_file("spinmask.out");
    for (auto e: spin_mask) {
      spinmask_file << static_cast<uint32_t>(e) << "\n";
    }
    spinmask_file.close();

    IntegerType bond_counter = 0;

    IntegerType k = 0;
    while(k < nz) {
      IntegerType j = 0;

      while (j < ny) {
        IntegerType i = 0;

        while (i < nx) {

          auto current_spin = (i % nx) + j * nx + k * nx * ny;
          auto next_spin_x = (i + 1) % nx + j * nx;
          auto next_spin_y = (current_spin + nx) % (nx * ny);
          auto next_spin_z = (current_spin + nx * ny) % (nx * ny * nz);

          // If the current spin doesn't exist, move on.
          if (!spin_mask[current_spin]) {
            i += 1;
            continue;
          }

          // If the next spin horizontally exists, add it to the bond,
          // increase bond counter.
          if (spin_mask[next_spin_x]) {
            bondsites.push_back(
                Bond{static_cast<uint16_t>(masked_spin_indices[current_spin]),
                    static_cast<uint16_t>(masked_spin_indices[next_spin_x])});
            bond_counter += 1;
          }

          if (spin_mask[next_spin_y]) {
            bondsites.push_back(
                Bond{static_cast<uint16_t>(masked_spin_indices[current_spin]),
                    static_cast<uint16_t>(masked_spin_indices[next_spin_y])});
            bond_counter += 1;
          }

          // Making a z=1 lattice a 2D lattice.
          if (nx > 1) {
            if (spin_mask[next_spin_z]) {
              bondsites.push_back(
                  Bond{static_cast<uint16_t>(masked_spin_indices[current_spin]),
                      static_cast<uint16_t>(masked_spin_indices[next_spin_z])});
              bond_counter += 1;
            }
          }

          i += 1;
        }

        j += 1;
      }

    k += 1;
    }

    nbonds = bond_counter;

    std::ofstream bonds_file("bonds.out");
    for (auto e: bondsites) {
      bonds_file << e.s1 << "," << e.s2 << "\n";
    }
    bonds_file.close();

    for (auto i = 0; i < n_active_sites; i++) {
      if (prng.randf() < 0.5) {
        spins.push_back(1);
      } else {
        spins.push_back(-1);
      }
    }

    initialized = true;
  }

  template <class PrngType> void make_lattice(PrngType &prng) {

    uint16_t intended_n_activesites = static_cast<uint16_t>(sim_input.n_active_sites);

    std::vector<FloatType> site_floatmask(sim_input.max_sites);
    for (auto i = 0; i < sim_input.max_sites; i++) {
      site_floatmask[i] = prng.randf();
    }
    float lower_bound_frac = 0.0;
    float upper_bound_frac = 1.0;
    float mid_point;
    

    IntegerType iterations = 0;

    do {
      mid_point = (upper_bound_frac + lower_bound_frac) / 2.0;
      n_active_sites = 0;

      std::for_each(site_floatmask.begin(), site_floatmask.end(),[&](float r){n_active_sites += r < mid_point? 1: 0;});
      if (n_active_sites > intended_n_activesites) {
        upper_bound_frac = mid_point;
      }
      else if (n_active_sites < intended_n_activesites) {
        lower_bound_frac = mid_point;
      } else {
        break;
      }

      iterations += 1;
    } while (n_active_sites != intended_n_activesites);

    std::transform(site_floatmask.begin(), site_floatmask.end(), spin_mask.begin(),
    [&](float r){return r < mid_point? 1: 0;});

  }

  bool initialized;
  uint16_t n_active_sites;
  uint16_t nbonds;
  SimulationInput<IntegerType, FloatType> sim_input;

  std::vector<bool> spin_mask;
  std::vector<int8_t> spins;
  std::vector<Bond> bondsites;
};