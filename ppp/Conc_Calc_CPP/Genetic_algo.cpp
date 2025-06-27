#ifndef PROFILING
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h> 
namespace nb = nanobind;
#endif

#include <omp.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "Calc_CPP.h"
#include "Params.h"

using namespace std;

int nt_global;

//Genetic Individual
struct Individual {
    int tau;
    float rho;
    float fitness;
};

//Global RNG setup
mt19937 rng(random_device{}()); //fixed seed for better reproducability
uniform_int_distribution<int> tau_dist(5, 200);
uniform_real_distribution<float> rho_dist(0.1f, 5.0f);
uniform_int_distribution<int> tau_mut_dist(-10, 10);
uniform_real_distribution<float> rho_mut_dist(-0.2f, 0.2f);

//fitness evaluation function
float evaluate(int tau, float rho) {
    params p = params(nt_global, tau, rho);
    auto [Ct, eps_field] = calc(p); //return ct and eps_field in tuple
    //calc total generated co2
    int width_CO2 = (p.nt-1) * 1 * p.ny * p.nx;  // species 1 = CO₂
    auto begin = Ct.begin() + width_CO2; //skipping past conc(CO)
    auto end = begin + (p.nx * p.ny);

    float total_CO2 = accumulate(begin, end, 0.0f);

    //calc total absorbed CO2
    float total_ABS_CO2 = accumulate(eps_field.begin(), eps_field.end(), 0.0f);
    
    if (total_CO2 < 1e-6f) return 0.0f;  // avoid division by near-zero
    return total_ABS_CO2 / total_CO2;
}

Individual mutate(const Individual& parent) {
    Individual child = parent;
    child.tau = clamp(child.tau + tau_mut_dist(rng), 5, 200); //clamp to make sure they stay in range
    child.rho = clamp(child.rho + rho_mut_dist(rng), 0.0f, 10.0f);
    return child;
}

tuple<int , int> run_ga(int nt_given) {//main
    nt_global = nt_given;
    const int population_size = 10;
    const int generations = 50;

    vector<Individual> population(population_size);

     // Step 1: Initial random population
    for (int i = 0; i < population_size; i++) {
        Individual ind{tau_dist(rng), rho_dist(rng), 0.0f};
        ind.fitness = evaluate(ind.tau, ind.rho);
        population[i] = ind; 
    }
    
    Individual best_all_time = population[0];

    for (int gen = 0; gen < generations; gen++) {
        //Step 2: Select top 3
        sort(population.begin(), population.end(),[](const Individual& a, const Individual& b) {
                      return a.fitness > b.fitness;
                  });

        vector<Individual> top3(population.begin(), population.begin() + 3);
        if (top3[0].fitness > best_all_time.fitness) {
            best_all_time = top3[0];
        }

        //Step 3: Create new generation by mutating top 3
        std::vector<Individual> next_gen(population_size);
        for (int i = 0; i < population_size; i++) {
            const Individual& parent = top3[i % 3]; //i % 3 gives 0, 1, 2, 0, 1, 2...
            Individual child = mutate(parent);
            child.fitness = evaluate(child.tau, child.rho);
            next_gen[i] = child;
        }

        population = std::move(next_gen);

        cout << "Gen " << gen + 1 << " Best fitness = " << top3[0].fitness << "  (tau = " << top3[0].tau << ", rho = " << top3[0].rho << ")\n";
    }

    std::cout << "\n Final best after " << generations << " generations:\n";
    std::cout << "Fitness = " << best_all_time.fitness << "  |  tau = " << best_all_time.tau << "  |  rho = " << best_all_time.rho << "\n";


    return {best_all_time.tau, best_all_time.rho};
}


#ifndef PROFILING
//make nb module
NB_MODULE(Genetic_algo, m) {
    m.def("Genetic_algo", &run_ga ,"genetic algo");
};
#endif