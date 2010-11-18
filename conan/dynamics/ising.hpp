/*
 * Conan - COmplex Network ANalisys
 * Copyright (C) 2008-2009  Ricardo Honorato Zimmer [rikardo.horo@gmail.com]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ISING_HPP
#define ISING_HPP
#include <conan/config.hpp>
#include <conan/utils.hpp>

#define K_B 0.08617 // constante de boltzmann en meV

#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;


namespace conan {

  /**
   * @brief A suitable particle for Ising model.
   */
  struct spin
  { // a suitable particle for Ising model, but easily can be generalized
    typedef std::vector<int> neighbors_vector;

    std::vector<double> J;        // the interaction constant (can be different for different neighbors)
    int position[3];              // where the particle is
    int value;                    // the value of spin usually +1 or -1
    neighbors_vector neighbors;   // a list of neighbors
    size_t num_neighbors;           // number of neighbors (with who interact)
  };

  /**
   * @brief Class designed for doing rewiring of a graph based on Ising models.
   */
  template <class Graph>
  struct ising
  {
    typedef typename ublas::matrix<decimal, ublas::column_major> cmatrix;

    ising(
        Graph& g
        )
    {
      _T = 0;
      adj_matrix = get_adj_matrix<Graph, cmatrix>(g);
      size = boost::num_vertices(g);
      num_particles = size_t((size * (size - 1)) / 2);
      system = std::vector<spin>(num_particles);
      size_t k = 0;
      for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < i; ++j)
        { // iterate over all adjacency matrix elements
          system[k].value = int(2 * adj_matrix(i, j) - 1);
          //system[k].num_neighbors = 2 * size - 4;
          system[k].neighbors = spin::neighbors_vector(system[k].num_neighbors);
          system[k].position[0] = i;
          system[k].position[1] = j;
          system[k++].position[2] = 0;
        }

      for (size_t i = 0; i < num_particles; ++i)
      {
        for (size_t j = 0; j < num_particles; ++j)
        {
          bool is_neighbor = false;
          if (i != j &&
             (system[i].position[0] == system[j].position[0] || // if column coincides
              system[i].position[1] == system[j].position[1] || // or row coincides
              system[i].position[0] == system[j].position[1] || // or the lower part of adj_matrix contains the information
              system[i].position[1] == system[j].position[0]    // idem...
             ))
            is_neighbor = true;

          if (is_neighbor)
          {
            system[i].neighbors.push_back(j);
            system[i].J.push_back(1);
          }
        }

        system[i].num_neighbors = system[i].neighbors.size();
      }

      return;
    }

    void equilibrate_system(
        decimal T,
        size_t num_steps = 50000
        )
    {
      _T = T;
      beta = 1.0 / (K_B * T);
      num_edges_now = 0;
      energy_now = 0;
      mag_now = 0;
      metropolis(num_steps);
      return;
    }

    size_t num_edges(
        size_t index = 0
        )
    {
      size_t _num_edges_size = _num_edges.size();
      if (_num_edges_size == 0)
        throw std::runtime_error("Montecarlo simulation has not been run");
      else if (_num_edges_size == 1)
        return _num_edges[0];
      else if (_num_edges_size < index)
        throw std::runtime_error("Index out of range");
      else
        return _num_edges[index];
    }

    decimal magnetization(
        size_t index = 0
        )
    {
      size_t _magnetization_size = _magnetization.size();
      if (_magnetization_size == 0)
        throw std::runtime_error("Montecarlo simulation has not been run");
      else if (_magnetization_size == 1)
        return _magnetization[0];
      else if (_magnetization_size < index)
        throw std::runtime_error("Index out of range");
      else
        return _magnetization[index];
    }

    decimal energy(
        size_t index = 0
        )
    {
      size_t _energy_size = _energy.size();
      if (_energy_size == 0)
        throw std::runtime_error("Montecarlo simulation has not been run");
      else if (_energy_size == 1)
        return _energy[0];
      else if (_energy_size < index)
        throw std::runtime_error("Index out of range");
      else
        return _energy[index];
    }

    Graph run_montecarlo(
        decimal T,
        size_t num_steps = 50000
        )
    {
      if (T != _T)
        std::cerr << "ising::run_montecarlo(T = " << T << ", num_steps = " << num_steps << "): Warning: "
                  << "Temperature for equilibrating the system was different from the montecarlo temperature" << std::endl;

      // initialize
      _num_edges.clear();
      _magnetization.clear();
      _energy.clear();

      _T = T;
      beta = 1.0 / (K_B * T);
      num_edges_now = 0;
      energy_now = 0;
      mag_now = 0;
      for (size_t i = 0; i < size; ++i)
      {
        mag_now += system[i].value;
        energy_now += hamiltonian(system[i].value, i);
        if (system[i].value == 1)
          ++num_edges_now;
      }
      energy_now /= 2;

      // run simulation
      metropolis(num_steps);

      _num_edges.push_back(num_edges_now);
      _magnetization.push_back(mag_now);
      _energy.push_back(energy_now);

      reconstruct_adj_matrix();
      return make_graph_from_adj_matrix<Graph, cmatrix>(adj_matrix);
    }

    std::vector<Graph> simulate(
        decimal T_init,
        decimal T_final,
        decimal deltaT,
        size_t num_steps = 50000,
        size_t num_steps_equilibrium = 50000
        )
    {
      std::vector<Graph> out_vec;
      _num_edges.clear();
      _magnetization.clear();
      _energy.clear();

      for (decimal T = T_init; T < T_final; T += deltaT)
      {
        // equilibrate system
        equilibrate_system(T, num_steps_equilibrium);

        // initialize
        _T = T;
        beta = 1.0 / (K_B * T);
        num_edges_now = 0;
        energy_now = 0;
        mag_now = 0;
        for (size_t i = 0; i < size; ++i)
        {
          mag_now += system[i].value;
          energy_now += hamiltonian(system[i].value, i);
          if (system[i].value == 1)
            ++num_edges_now;
        }
        energy_now /= 2;

        // run simulation
        metropolis(num_steps);

        reconstruct_adj_matrix();
        out_vec.push_back(make_graph_from_adj_matrix<Graph, cmatrix>(adj_matrix));

        _num_edges.push_back(num_edges_now);
        _magnetization.push_back(mag_now);
        _energy.push_back(energy_now);
      }
      return out_vec;
    }

  protected:
    void metropolis(size_t num_steps)
    {
      // create random number generators for metropolis-montecarlo
      typename boost::uniform_real<decimal> percent(0, 1);
      typename boost::uniform_int<size_t> int_dist(0, num_particles - 1);
#ifdef DETERMINISTIC_RNG
      typename boost::mt19937 rng(time(NULL) - num_steps - num_particles - int(beta * 100));
      typename boost::variate_generator< boost::mt19937&, boost::uniform_real<decimal> >
        probability(rng, percent);
      typename boost::variate_generator< boost::mt19937&, boost::uniform_int<size_t> >
        random_target(rng, int_dist);
#else
      typename boost::random_device rng;
      typename boost::variate_generator< boost::random_device&, boost::uniform_real<decimal> >
        probability(rng, percent);
      typename boost::variate_generator< boost::random_device&, boost::uniform_int<size_t> >
        random_target(rng, int_dist);
#endif

      for (size_t step = 0; step < num_steps; ++step)
      {
        for (size_t i = 0; i < num_particles; ++i)
        {
          size_t target = random_target(); // choose a random target
          int new_value = -system[target].value;
          decimal E_old = hamiltonian(system[target].value, target),
                  E_new = hamiltonian(new_value, target),
                  flip = 1.0; // probability of spin flipping, set -by default- to change

          if (E_old < E_new)
            flip = exp(beta * (E_old - E_new));

          if ( (flip - probability()) >= 0 )
          {
            system[target].value = -system[target].value; // update the spin value
            mag_now += -system[target].value + new_value;
            energy_now += E_new - E_old;
            if (system[target].value == 1)
              ++num_edges_now;
            else
              --num_edges_now;
          }

        } // for i < num_particles
      } // for step < neq
      return;
    }

    void reconstruct_adj_matrix()
    {
      for (size_t i = 0; i < num_particles; ++i)
      {
        size_t k = system[i].position[0],
               l = system[i].position[1];
        adj_matrix(l, k) = adj_matrix(k, l) = (system[i].value + 1) / 2;
      }
      return;
    }

    decimal hamiltonian(int value, size_t target)
    {
      decimal E = 0;
      for (size_t i = 0; i < system[target].num_neighbors; ++i)
      { // E = sum(J * S1 * S2) over all neighbors
        int j = system[target].neighbors[i]; // index of i-th neighbor
        E += -system[target].J[i] * value * system[j].value;
      }
      return E;
    }

    decimal _T; // last temperature
    size_t num_edges_now;
    decimal energy_now;
    decimal mag_now;

    std::vector<size_t> _num_edges;
    std::vector<decimal> _magnetization;
    std::vector<decimal> _energy;

  public:
    size_t num_particles; // Ntotal
    decimal beta;	      // 1 / (Boltzmann constant * _T)
    size_t size;      	  // size of the adjacency matrix
    cmatrix adj_matrix;
    std::vector<spin> system;
  };

} // conan

#endif // ISING_HPP
