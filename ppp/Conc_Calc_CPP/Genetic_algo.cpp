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
#include <fstream>

#include <sys/resource.h> //parallelisation and optimization too good pc crashed at over temp.

#include "Calc_CPP.h"
#include "Params.h"

using namespace std;

int nt_global;

//mem limiter
void limit_memory(std::size_t max_bytes) {
    struct rlimit rl;
    rl.rlim_cur = max_bytes;
    rl.rlim_max = max_bytes;
    setrlimit(RLIMIT_AS, &rl);
}

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
float evaluate(int tau, float rho, ofstream& log_file) {
    omp_set_num_threads(2);
    params p = params(nt_global, tau, rho);
    auto [Ct, eps_field] = calc(p); //return ct and eps_field in tuple
    //stability checks
    // NaNs or Infs in Ct or eps_field or overflow values
    bool is_unstable = !std::all_of(Ct.begin(), Ct.end(), [](float x) { return std::isfinite(x); }) || !std::all_of(eps_field.begin(), eps_field.end(), [](float x) { return std::isfinite(x); }) || (*std::max_element(Ct.begin(), Ct.end()) > 1e6f);
    float total_CO2 = 0.0f;
    float total_ABS_CO2 = 0.0f;
    //abort when unstable
    if (is_unstable) {
        #pragma omp critical
        {
        log_file << "[UNSTABLE] tau = " << tau << ", rho = " << rho << " fitness = 0 (due to numerical instability)" << "\n";
        }
        return 0;
    }
    //performace eval
    //calc total generated co2
    if(!is_unstable){
        for (int t = 0; t < p.nt; t++){
            size_t ct_base   = ((t * p.ns + 1) * p.ny) * p.nx; //index ofset per t
            size_t eps_base  = t * p.ny * p.nx;

            auto begin_ct  = Ct.begin() + ct_base;
            auto end_ct    = begin_ct + (p.nx * p.ny);
            auto begin_eps = eps_field.begin() + eps_base;

            total_CO2 += std::accumulate(begin_ct, end_ct, 0.0f) * p.dt;
            total_ABS_CO2 += std::inner_product(begin_ct, end_ct, begin_eps, 0.0f) * p.dt; //calculate actual abs. doesn't need an end point for eps => bound to ct
            if (t % 1000 == 0){
                #pragma omp critical
                {
                log_file << "[t=" << t << "] total_CO2=" << total_CO2 << ", total_ABS_CO2=" << total_ABS_CO2 << ", fitness=" << (total_ABS_CO2 / total_CO2) << "\n";
                }
            }
        }
    }

    
    if (total_CO2 < 1e-6f) return 0.0f;  // avoid division by near-zero
    return total_ABS_CO2 / total_CO2;
}

Individual mutate(const Individual& parent) {
    thread_local mt19937 local_rng(random_device{}());
    thread_local uniform_int_distribution<int> tau_dist(-10, 10);
    thread_local uniform_real_distribution<float> rho_dist(-0.2f, 0.2f);

    Individual child = parent;
    child.tau = clamp(child.tau + tau_dist(local_rng), 5, 200);
    child.rho = clamp(child.rho + rho_dist(local_rng), 0.0f, 10.0f);
    return child;
};


tuple<int , int> run_ga(int nt_given) {//main
    limit_memory(15L * 1024 * 1024 * 1024); // 15 GB cap
    nt_global = nt_given;
    ofstream log_file("ga_log.txt");
    const int population_size = 10;
    const int generations = 25;

    vector<Individual> population(population_size);

     // Step 1: Initial random population
    #pragma omp parallel for shared(population)
    for (int i = 0; i < population_size; i++) {
        Individual ind{tau_dist(rng), rho_dist(rng), 0.0f};
        ind.fitness = evaluate(ind.tau, ind.rho, log_file);
        #pragma omp critical
        {
        population[i] = ind; 
        }
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
        #pragma omp parallel for shared(next_gen, log_file) firstprivate(gen, top3) // 
        for (int i = 0; i < population_size; i++) {
            const Individual& parent = top3[i % 3]; //i % 3 gives 0, 1, 2, 0, 1, 2...
            Individual child = mutate(parent);
            child.fitness = evaluate(child.tau, child.rho, log_file);
            #pragma omp critical 
            {
            next_gen[i] = child;
            log_file << "gen " << gen << ", child " << i << ", tau = " << child.tau << ", rho = " << child.rho << ", fitness = " << child.fitness << "\n";
            }
        }

        population = std::move(next_gen);
        

        log_file << "Gen " << gen + 1 << " Best fitness = " << top3[0].fitness << "  (tau = " << top3[0].tau << ", rho = " << top3[0].rho << ")\n";
    }

    log_file << "\n Final best after " << generations << " generations:\n";
    log_file << "Fitness = " << best_all_time.fitness << "  |  tau = " << best_all_time.tau << "  |  rho = " << best_all_time.rho << "\n";
    log_file.close();


    return {best_all_time.tau, best_all_time.rho};
}


#ifndef PROFILING
//make nb module
NB_MODULE(Genetic_algo, m) {
    m.def("Genetic_algo", &run_ga ,"genetic algo");
};
#endif