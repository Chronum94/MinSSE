#pragma once

#include <algorithm>
#include <cassert>
#include <numeric>

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
    expansion_cutoff = nbeta > 4 ? nbeta : 4;
    current_opcount = 0;
    aprob = 0.5 * lattice.nbonds * sim_input.beta;
    dprob = 1.0 / aprob;

    ususc = 0.0;
    energy1 = 0.0;

    opstring.resize(expansion_cutoff, Op{IDENTITY, 0});
    first_spinop.resize(lattice.n_active_sites, -1);
    last_spinop.resize(lattice.n_active_sites, -1);
    vertexlist.resize(4 * expansion_cutoff, -1);

    spins_copy.resize(lattice.n_active_sites, 0);
    snapshot_acc.resize(lattice.n_active_sites, 0);
    local_susceptibilities.resize(lattice.n_active_sites, 0.0);
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
  std::vector<FloatType> local_susceptibilities;
  std::vector<int8_t> spins_copy;
  std::vector<IntegerType> snapshot_acc;
  IntegerType intrasweep_measurement_interval;
  IntegerType local_susc_measurement_count = 0;

  uint64_t coin_pool = 0;
  uint32_t coin_bits_left = 0;

  // Batched coin flips: one PRNG draw serves 64 fair coin decisions.
  _INLINE bool coin_flip() {
    if (coin_bits_left == 0) {
      coin_pool = prng.rand_64();
      coin_bits_left = 64;
    }
    bool heads = coin_pool & 1;
    coin_pool >>= 1;
    coin_bits_left -= 1;
    return heads;
  }

  std::vector<Op> opstring;
  std::vector<IntegerType> first_spinop;
  std::vector<IntegerType> last_spinop;
  std::vector<IntegerType> vertexlist;

  // Fused diagonal update + vertex linking: the op type at index i is final
  // once the diagonal update passes it, so we can build the linked vertex
  // list in the same sweep instead of re-walking the opstring.
  void diagonal_update_and_link() {
    std::fill(first_spinop.begin(), first_spinop.end(), -1);
    std::fill(last_spinop.begin(), last_spinop.end(), -1);

    // Hoist members into locals/raw pointers: stores through the vectors
    // would otherwise force reloads of prng state and counters (aliasing).
    auto local_prng = prng;
    auto opcount = current_opcount;
    const auto cutoff = expansion_cutoff;
    const auto local_aprob = aprob;
    const uint64_t aprob_x2_32 =
        static_cast<uint64_t>(local_aprob * 4294967296.0);
    const uint64_t nbonds = lattice.nbonds;
    Op *const ops = opstring.data();
    const Bond *const bonds = lattice.bondsites.data();
    int8_t *const spins = lattice.spins.data();
    IntegerType *const vlist = vertexlist.data();
    IntegerType *const firstop = first_spinop.data();
    IntegerType *const lastop = last_spinop.data();

    for (uint32_t op_index = 0; op_index < cutoff; ++op_index) {
      IntegerType bond_index;

      if (ops[op_index].optype == IDENTITY) {
        // Lemire-style bounded pick: avoids float multiply + float->int cast.
        bond_index = static_cast<IntegerType>(
            (static_cast<uint64_t>(local_prng.rand_32()) * nbonds) >> 32);
        auto bond = bonds[bond_index];
        // All-integer acceptance test: u/2^32 * (M-n) <= aprob, scaled by 2^32.
        if (spins[bond.s1] != spins[bond.s2] &&
            static_cast<uint64_t>(local_prng.rand_32()) * (cutoff - opcount) <=
                aprob_x2_32) {
          ops[op_index].optype = DIAGONAL;
          ops[op_index].opbond = bond_index;
          opcount += 1;
        } else {
          vlist[op_index << 2] = -1;
          continue;
        }
      }

      else if (ops[op_index].optype == DIAGONAL) {
        if (local_aprob * local_prng.randf() <= (cutoff - opcount + 1)) {
          ops[op_index].optype = IDENTITY;
          opcount -= 1;
          vlist[op_index << 2] = -1;
          continue;
        }
        bond_index = ops[op_index].opbond;
      }

      else {
        bond_index = ops[op_index].opbond;
        spins[bonds[bond_index].s1] *= -1;
        spins[bonds[bond_index].s2] *= -1;
      }

      // Link this (non-identity) op's four vertex legs.
      IntegerType v0 = op_index << 2;
      auto s1 = bonds[bond_index].s1;
      auto s2 = bonds[bond_index].s2;

      auto v1 = lastop[s1];
      auto v2 = lastop[s2];

      if (v1 != -1) {
        vlist[v1] = v0;
        vlist[v0] = v1;
      } else {
        firstop[s1] = v0;
      }

      if (v2 != -1) {
        vlist[v2] = v0 + 1;
        vlist[v0 + 1] = v2;
      } else {
        firstop[s2] = v0 + 1;
      }

      lastop[s1] = v0 + 2;
      lastop[s2] = v0 + 3;
    } // endloop

    prng = local_prng;
    current_opcount = opcount;

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
    IntegerType *const vlist = vertexlist.data();
    Op *const ops = opstring.data();

    for (IntegerType v0 = 0; v0 < 4 * (IntegerType)expansion_cutoff; v0 += 4) {
      if (vlist[v0] < 0) {
        continue;
      }

      IntegerType v1 = v0;
      if (coin_flip()) {
        do {
          ops[v1 >> 2].optype =
              static_cast<optype_t>(ops[v1 >> 2].optype ^ OFF_DIAGONAL);
          vlist[v1] = -2;
          IntegerType v2 = v1 ^ 1;
          v1 = vlist[v2];
          vlist[v2] = -2;
        } while (v1 != v0);
      } else {
        do {
          vlist[v1] = -1;
          IntegerType v2 = v1 ^ 1;
          v1 = vlist[v2];
          vlist[v2] = -1;
        } while (v1 != v0);
      }
    }

    for (auto i = 0; i < lattice.n_active_sites; i++) {
      auto first_op = first_spinop[i];
      if (first_op != -1) {
        if (vertexlist[first_op] == -2) {
          lattice.spins[i] *= -1;
        }
      } else {
        if (coin_flip()) {
          lattice.spins[i] *= -1;
        }
      }
    }
  }

  void adjust_expansion_cutoff_linear(uint32_t i) {

    auto new_cutoff = static_cast<uint32_t>(4.0 / 3.0 * current_opcount);
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
    
    FloatType sum_of_spins = static_cast<FloatType>(
        std::accumulate(lattice.spins.begin(), lattice.spins.end(), 0));
    // std::cout << "\t" << temp << "\n";
    
    ususc += sum_of_spins * sum_of_spins / 4.0;
    energy1 += static_cast<FloatType>(current_opcount);


    // We would like to not flip the original spins... I think.
    std::copy(lattice.spins.begin(), lattice.spins.end(), spins_copy.begin());
    // Accumulate spin snapshots in integers; multiply by sum_of_spins once
    // per sweep instead of once per snapshot (all values are exact ints).
    std::fill(snapshot_acc.begin(), snapshot_acc.end(), 0);
    const Op *const ops = opstring.data();
    const Bond *const bonds = lattice.bondsites.data();
    int8_t *const scopy = spins_copy.data();
    IntegerType *const acc = snapshot_acc.data();
    const IntegerType nsites = lattice.n_active_sites;
    const uint32_t cutoff = expansion_cutoff;

    IntegerType countdown = 0; // fires at op_index 0, matching op_index % interval == 0
    for (uint32_t op_index = 0; op_index < cutoff; op_index++) {
      const auto op = ops[op_index];
      if (op.optype == OFF_DIAGONAL) {
        auto bond = bonds[op.opbond];
        scopy[bond.s1] *= -1;
        scopy[bond.s2] *= -1;
      }

      if (countdown == 0) {
        countdown = intrasweep_measurement_interval;
        local_susc_measurement_count += 1;
        for (IntegerType i = 0; i < nsites; i++) {
          acc[i] += scopy[i];
        }
      }
      countdown -= 1;
    }

    for (auto i = 0; i < lattice.n_active_sites; i++) {
      local_susceptibilities[i] += snapshot_acc[i] * sum_of_spins;
    }
  }

  void finalize_and_write_results(std::ofstream &outfile, std::ofstream &locsuscfile, unsigned int bin_number) {
    ususc /= sim_input.msteps;
    ususc *= sim_input.beta / lattice.n_active_sites;

    energy1 /= sim_input.msteps;
    energy1 = -energy1 / (sim_input.beta * lattice.n_active_sites) +
              (0.25 * lattice.nbonds) / lattice.n_active_sites;

    

    outfile.width(9);
    outfile.precision(7);
    outfile << std::fixed << bin_number << " " << ususc << " " << energy1
            << "\n";

    locsuscfile.width(8);
    locsuscfile.precision(6);



    for (auto e: local_susceptibilities) {
      e /= local_susc_measurement_count;
      auto f = e * sim_input.beta / (4.0);
      locsuscfile << f << " ";
    }
    locsuscfile << "\n";
    if (bin_number % 10 == 0) {
      outfile.flush();
      locsuscfile.flush();
    }

    ususc = 0.0;
    energy1 = 0.0;
    local_susc_measurement_count = 0;
    std::fill(local_susceptibilities.begin(), local_susceptibilities.end(), 0.0);
  }
  void run(std::ofstream &outputfile) {

    outputfile << "Initial cutoff: " << expansion_cutoff << "\n";
    for (auto j = 0; j < sim_input.isteps; j++) {
      diagonal_update_and_link();
      loop_update();
      adjust_expansion_cutoff_linear(j);
    }

    FloatType interval = lattice.n_active_sites / 4.0;
    intrasweep_measurement_interval =  interval > 1 ? static_cast<IntegerType>(interval) : 1;
    outputfile << "Final cutoff: " << expansion_cutoff << "\n";
    outputfile << "Locsusc m.interval: " << intrasweep_measurement_interval << "\n";

    std::ofstream results_outstream("results.dat", std::ofstream::out);
    std::ofstream local_susc_outstream("locsusc.dat", std::ofstream::out);
    results_outstream.width(10);
    results_outstream.precision(7);

    results_outstream << "Bin\t"
            << "Uni. Susc.\t"
            << "Tot. en.\n";

    for (auto i = 0; i < sim_input.nbins; i++) {
      for (auto j = 0; j < sim_input.msteps; j++) {
        diagonal_update_and_link();
        loop_update();
        measure();
      }
      finalize_and_write_results(results_outstream, local_susc_outstream, i);
    }

    results_outstream.close();
  }
};