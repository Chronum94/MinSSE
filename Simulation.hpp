#pragma once

#include "Lattice.hpp"
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
             PrngType _prng) {
    sim_input = in;
    lattice = lat;
    prng = _prng;

    expansion_cutoff = sim_input.beta_n_integer;
    current_opcount = 0;
    aprob = 0.5 * sim_input.beta_n;
    dprob = 1.0 / aprob;

    optype.resize(expansion_cutoff, 0);
    opbond.resize(expansion_cutoff, 0);

    first_spinop.resize(lattice.active_sites, -1);
    last_spinop.resize(lattice.active_sites, -1);

    vertexlist.resize(4 * expansion_cutoff, -1);
  }
  SimulationInput<IntegerType, FloatType> sim_input;
  LatticeType lattice;

  PrngType prng;
  IntegerType expansion_cutoff;
  IntegerType current_opcount;
  FloatType aprob;
  FloatType dprob;

  // 0 for identity, 1 for diagonal, 2 for off-diagonal.
  std::vector<IntegerType> optype;
  std::vector<IntegerType> opbond;

  std::vector<int8_t> first_spinop;
  std::vector<int8_t> last_spinop;

  std::vector<IntegerType> vertexlist;

  void diagonal_update() {

    for (uint32_t op_index = 0; op_index < expansion_cutoff; ++op_index) {
      auto bond = static_cast<IntegerType>(prng.randf() * lattice.nbonds);
      if (optype[op_index] == 0) {
        // auto bondsites = lattice.bondsites;
        if (lattice.spins[lattice.bondsites[bond].s1] != lattice.spins[lattice.bondsites[bond].s2]) {
          
          if (prng.randf() * (expansion_cutoff - current_opcount) <=
          aprob) {
            optype[op_index] = 1; // diagonal operation.
            opbond[op_index] = bond;
            current_opcount += 1;
            std::cout << "A, opcount: " << current_opcount << "\n";
          }
        }
      }

      else if (optype[op_index] == 1)
      {
        if (prng.randf() <= (expansion_cutoff - current_opcount) * dprob) {
          optype[op_index] = 0;
          // opbond[i] = 0; This can be left alone, since optype[i] == 0 is
          // an
          // identity op.
          current_opcount -= 1;
          std::cout << "D, opcount: " << current_opcount << "\n";
        }
      }

      // If the optype is off-diagonal.
      else
      {
        flip_spins(lattice.bondsites[bond]);
      }
    } // endloop
  }

  void link_vertices() {

    for (uint32_t v0 = 0; v0 < 4 * expansion_cutoff; v0 += 4) {
      auto op_index = v0 >> 2;
      auto bond = opbond[op_index];

      if (optype[op_index] != 0) {
        auto s1 = lattice.bondsites[bond].s1;
        auto s2 = lattice.bondsites[bond].s2;

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
        vertexlist[v0] = -1;
        vertexlist[v0 + 1] = -1;
        vertexlist[v0 + 2] = -1;
        vertexlist[v0 + 3] = -1;
      }
    }

    for (auto s1 = 0; s1 < lattice.active_sites; s1 += 1) {
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

    for (auto i = 0; i < lattice.active_sites; i += 1) {
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

  _INLINE void flip_loop(IntegerType v0, IntegerType v1) {
    do {
      optype[v1 >> 2] ^= 1;
      vertexlist[v1] = -2;
      IntegerType v2 = v1 ^ 1;
      v1 = vertexlist[v2];
      vertexlist[v2] = -2;
    } while (v1 != v0);
  }

  _INLINE void visit_loop(IntegerType v0, IntegerType v1) {
    do {
      vertexlist[v1] = -1;
      IntegerType v2 = v1 ^ 1;
      v1 = vertexlist[v2];
      vertexlist[v2] = -1;
    } while (v1 != v0);
  }

  _INLINE void flip_spins(Bond<IntegerType> &b) {
    lattice.spins[b.s1] *= -1;
    lattice.spins[b.s2] *= -1;
  }

  void run() {
    std::cout << "Running..." << std::endl;
    for (auto i = 0; i < sim_input.nbins; i++) {
      std::cout << i << std::endl;
      diagonal_update();
      // link_vertices();
      // loop_update();
    }
  }
};