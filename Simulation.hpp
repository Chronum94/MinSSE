#pragma once

#include <algorithm>
#include <numeric>
#include <cassert>

#include "Bond.hpp"
#include "Lattice.hpp"
#include "Op.hpp"
#include "SimulationInput.hpp"

#ifdef _MSC_VER
#define _INLINE __forceinline
#elif __GNUC__
#define _INLINE __attribute__((always_inline)) inline
#else
#define _INLINE inline
#endif

template <class IntegerType, class FloatType, class LatticeType, class PrngType>
struct Simulation {
  Simulation(SimulationInput<IntegerType, FloatType> in, LatticeType lat,
             PrngType &_prng) {
    sim_input = in;
    lattice = lat;
    prng = _prng;

    auto nbeta = lattice.n_active_sites * sim_input.beta;
    expansion_cutoff =  nbeta > 4? nbeta : 4;
    current_opcount = 0;
    aprob = 0.5 * lattice.nbonds * sim_input.beta;
    dprob = 1.0 / aprob;

    ususc = 0.0;
    energy1 = 0.0;
    specific_heat = 0.0;

    opstring.resize(expansion_cutoff, Op{IDENTITY, 0});
    first_spinop.resize(lattice.n_active_sites, -1);
    last_spinop.resize(lattice.n_active_sites, -1);
    vertexlist.resize(4 * expansion_cutoff, -1);
  }
  SimulationInput<IntegerType, FloatType> sim_input;
  LatticeType lattice;

  PrngType prng;
  uint32_t expansion_cutoff;
  uint32_t current_opcount;
  FloatType aprob;
  FloatType dprob;

  FloatType ususc;
  FloatType energy1;
  FloatType specific_heat;

  std::vector<Op> opstring;
  std::vector<IntegerType> first_spinop;
  std::vector<IntegerType> last_spinop;
  std::vector<IntegerType> vertexlist;

  std::vector<uint32_t> equilibriation_opcounts;

  void diagonal_update() {

    for (uint32_t op_index = 0; op_index < expansion_cutoff; ++op_index) {
      if (opstring[op_index].optype == IDENTITY) {
        auto bond_index = static_cast<IntegerType>(prng.randf() * lattice.nbonds);
        if (lattice.spins[lattice.bondsites[bond_index].s1] !=
            lattice.spins[lattice.bondsites[bond_index].s2]) {

          if (prng.randf() * (expansion_cutoff - current_opcount) <= aprob) {
            opstring[op_index].optype = DIAGONAL;
            opstring[op_index].opbond = bond_index;
            current_opcount += 1;
          }
        }
      }

      else if (opstring[op_index].optype == DIAGONAL) {
        if (aprob * prng.randf()<= (expansion_cutoff - current_opcount + 1)) {
          opstring[op_index].optype = IDENTITY;
          current_opcount -= 1;
        }
      }

      else {
        auto bond_index = opstring[op_index].opbond;
        flip_spins(lattice.bondsites[bond_index]);
      }
    } // endloop
  }

  void link_vertices() {
    std::fill(first_spinop.begin(), first_spinop.end(), -1);
    std::fill(last_spinop.begin(), last_spinop.end(), -1);

    for (int v0 = 0; v0 < 4 * expansion_cutoff; v0 += 4) {
      auto op_index = v0 >> 2;

      if (opstring[op_index].optype != IDENTITY) {
        auto bond_index = opstring[op_index].opbond;

        auto s1 = lattice.bondsites[bond_index].s1;
        auto s2 = lattice.bondsites[bond_index].s2;

        auto v1 = last_spinop[s1];
        auto v2 = last_spinop[s2];

        if (v1 != -1) {
          vertexlist[v1] = v0;
          vertexlist[v0] = v1;
        } else {
          first_spinop[s1] = v0;
        }

        if (v2 != -1) {
          vertexlist[v2] = v0 + 1;
          vertexlist[v0 + 1] = v2;
        } else {
          first_spinop[s2] = v0 + 1;
        }

        last_spinop[s1] = v0 + 2;
        last_spinop[s2] = v0 + 3;
      }

      else {
        std::fill(vertexlist.begin() + v0, vertexlist.begin() + v0 + 4, -1);
      }
    }

    for (auto s1 = 0; s1 < lattice.n_active_sites; s1++) {
      auto v1 = first_spinop[s1];
      if (v1 != -1) {
        auto v2 = last_spinop[s1];
        vertexlist[v2] = v1;
        vertexlist[v1] = v2;
      }
    }
  }

  void loop_update() {
    for (auto v0 = 0; v0 < 4 * expansion_cutoff; v0 += 4) {
      if (vertexlist[v0] < 0) {
        continue;
      }

      auto v1 = v0;
      if (prng.randf() < 0.5) {
        flip_loop(v0, v1);
      } else {
        visit_loop(v0, v1);
      }
    }

    for (auto i = 0; i < lattice.n_active_sites; i++) {
      auto first_op = first_spinop[i];
      if (first_op != -1) {
        if (vertexlist[first_op] == -2) {
          lattice.spins[i] *= -1;
        }
      } else {
        if (prng.randf() < 0.5) {
          lattice.spins[i] *= -1;
        }
      }
    }
  }

  void adjust_expansion_cutoff_linear(uint32_t i) {

    auto new_cutoff = static_cast<uint32_t>(4.0/3.0 * current_opcount);
    if (new_cutoff < expansion_cutoff) {
      return;
    }

    else {
      expansion_cutoff = new_cutoff;
      // std::cout << i << "\t" << expansion_cutoff << "\n";
      opstring.resize(expansion_cutoff, Op{IDENTITY, 0});

      vertexlist.resize(4 * expansion_cutoff, -1);
    }
  }

  void measure() {
    FloatType temp = static_cast<FloatType>(std::accumulate(lattice.spins.begin(), lattice.spins.end(), 0));
    ususc += temp * temp / 4.0;

    energy1 += static_cast<FloatType>(current_opcount);
  }

  void write_results(std::ofstream &outfile, unsigned int bin_number) {
    ususc /= sim_input.msteps;
    ususc*= sim_input.beta / lattice.n_active_sites;

    energy1 /= sim_input.msteps;
    energy1 = - energy1 / (sim_input.beta * lattice.n_active_sites) + (0.25 * lattice.nbonds) / lattice.n_active_sites;


    outfile.width(10);
    outfile.precision(7);
    outfile << std::fixed << bin_number << "\t" << ususc << "\t" << energy1 << "\n";

    if (bin_number % 10 == 0) {
      outfile.flush();
    }

    ususc = 0.0;
    energy1 = 0.0;
  }
  _INLINE void flip_loop(IntegerType v0, IntegerType v1) {
    do {
      assert(v0 >= 0);
      assert(v1 >= 0);
      opstring[v1 >> 2].optype =
          static_cast<optype_t>(opstring[v1 >> 2].optype ^ OFF_DIAGONAL);

      vertexlist[v1] = -2;
      IntegerType v2 = v1 ^ 1;
      v1 = vertexlist[v2];
      vertexlist[v2] = -2;
    } while (v1 != v0);
  }

  _INLINE void visit_loop(IntegerType v0, IntegerType v1) {
    do {
      assert(v0 >= 0);
      assert(v1 >= 0);
      vertexlist[v1] = -1;
      IntegerType v2 = v1 ^ 1;
      v1 = vertexlist[v2];
      vertexlist[v2] = -1;
    } while (v1 != v0);
  }

  _INLINE void flip_spins(Bond &b) {
    lattice.spins[b.s1] *= -1;
    lattice.spins[b.s2] *= -1;
  }

  void run(std::ofstream &outputfile) {

    outputfile << "Initial cutoff: " << expansion_cutoff << "\n";
    for (auto j = 0; j < sim_input.isteps; j++) {
      diagonal_update();
      link_vertices();
      loop_update();
      adjust_expansion_cutoff_linear(j);
    }

    outputfile << "Final cutoff: " << expansion_cutoff << "\n";

    std::ofstream results(sim_input.resultsfile, std::ofstream::out);
    results.width(10);
    results.precision(7);
    
    results << "Bin\t" << "Uni. Susc.\t" << "Tot. en.\n";

    for (auto i = 0; i < sim_input.nbins; i++) {
      for (auto j = 0; j < sim_input.msteps; j++) {
        diagonal_update();
        link_vertices();
        loop_update();
        measure();
      }
      write_results(results, i);
    }

    results.close();
  }
};