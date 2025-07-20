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
#include <mpi.h>
#include <sys/resource.h>

#include "Calc_CPP.h"
#include "Params.h"

using namespace std;

int nt_global;

static struct {
    int tau_min = 1;
    int tau_max = 100;

    float rho_min = -2.0f; //20
    float rho_max = 2.0f;  //20

    float rho_lower_clamp = 0.1f;
    float rho_upper_clamp = 2.0f;

    int tau_clamp_min = 1;
    int tau_clamp_max = 100;
} cfg;

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

//fitness evaluation function
float evaluate(params& p, int tau, float rho, ofstream& log_file) { //params only intialized once
    //omp_set_num_threads(2);
    p.setTau(tau);
    p.setRho(rho);
    auto [Ct, eps_field] = calc(p); //return ct and eps_field in tuple
    //stability checks
    // NaNs or Infs in Ct or eps_field or overflow values
    bool is_unstable = !std::all_of(Ct.begin(), Ct.end(), [](float x) { return std::isfinite(x); }) || !std::all_of(eps_field.begin(), eps_field.end(), [](float x) { return std::isfinite(x); }) || (*std::max_element(Ct.begin(), Ct.end()) > 1e6f);
    float total_CO2 = 0.0f;
    float total_ABS_CO2 = 0.0f;
    //abort when unstable and return fitness 0
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
            #ifdef PROFILING
            if (t % 1000 == 0){
                #pragma omp critical
                {
                log_file << "[t=" << t << "] total_CO2=" << total_CO2 << ", total_ABS_CO2=" << total_ABS_CO2 << ", fitness=" << (total_ABS_CO2 / total_CO2) << "\n";
                }
            }
            #endif
        }
    }
    // avoid division by near-zero
    if (total_CO2 < 1e-6f) return 0.0f;
    //discourage very small rho's
    if (rho <= 0.05f) {
        #pragma omp critical
        {
        log_file << "[DISCOURAGED] rho too small => forcing fitness = 0\n";
        }
    return 0.0f;
    }

    return total_ABS_CO2 / total_CO2;
}

Individual mutate(const Individual& parent) {
    thread_local mt19937 local_rng(random_device{}());
    thread_local uniform_int_distribution<int> tau_dist(cfg.tau_min, cfg.tau_max);
    thread_local uniform_real_distribution<float> rho_dist(cfg.rho_min, cfg.rho_max);

    Individual child = parent;
    child.tau = clamp(child.tau + tau_dist(local_rng), cfg.tau_clamp_min, cfg.tau_clamp_max);
    child.rho = clamp(child.rho + rho_dist(local_rng), cfg.rho_lower_clamp, cfg.rho_upper_clamp);
    return child;
};

Individual crossover(const Individual& a, const Individual& b) {
    thread_local mt19937 local_rng(random_device{}());
    uniform_real_distribution<float> dist(0.0f, 1.0f);
    Individual child;
    //crossover
    float ratio = dist(local_rng);
    child.tau = static_cast<int>(ratio * a.tau + (1.0f - ratio) * b.tau);
    child.rho = ratio * a.rho + (1.0f - ratio) * b.rho;
    return child;
}   

//top 4 mpi exchanger function
void exchange_top_individuals(vector<Individual>& population, int rank) {
    //extract top 4 locally
    vector<Individual> top4_local(population.begin(), population.begin() + 4);

    //buffer for sending: tau, rho, fitness per individual
    float send_buf[12]; //array init
    for (int i = 0; i < 4; i++) {
        send_buf[i * 3 + 0] = static_cast<float>(top4_local[i].tau); //tau is an int
        send_buf[i * 3 + 1] = top4_local[i].rho;
        send_buf[i * 3 + 2] = top4_local[i].fitness;
    }
    float recv_buf[12];

    int peer = rank ^ 1; //bitwise XOR operator to switch between ranks to talk to other rank => switches 1 and 0

    //exchange
    MPI_Sendrecv(send_buf, 12, MPI_FLOAT, peer, 0, recv_buf, 12, MPI_FLOAT, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //merge and select best 4
    array<Individual, 8> merged; //stack allocation
    for (int i = 0; i < 4; i++) {
        merged[i] = top4_local[i];
    }
    for (int i = 0; i < 4; i++) {
        merged[i + 4] = {static_cast<int>(recv_buf[i * 3 + 0]), recv_buf[i * 3 + 1], recv_buf[i * 3 + 2]}; //assigning i 4=>7 wiht the 4 recieved individuals (assing the tau,rho, fitness to the individuals inside the call)
    };

    //sort and update pop
    sort(merged.begin(), merged.end(), [](const Individual& a, const Individual& b) {return a.fitness > b.fitness;});
    for (int i = 0; i < 4; i++) {
        population[i] = merged[i];
    }
}

tuple<int , float> run_ga(params& p) {//main
    limit_memory(15L * 1024 * 1024 * 1024); // 15 GB cap
    //init mpi
    int initialized = 0;
    MPI_Initialized(&initialized); //checks if MPI has been init, set inititialzed to 1 if so
    if (!initialized) {
        MPI_Init(nullptr, nullptr); //prevent from double init due to nanobind python call
    }   
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //prevent accidental runs with wrong rank (2 processes should be used)
    if (size != 2) {
        if (rank == 0) cerr << "use mpirun -np 2 \n"; //error message
        MPI_Finalize();
        return {-1, -1.0f}; //neg values break the simulation visual
    }

    nt_global = p.nt;
    ofstream log_file("ga_log.txt");
    const int population_size = 10;
    const int generations = 50;

    vector<Individual> population(population_size);

     // Step 1: Initial random population
    #pragma omp parallel for shared(population)
    for (int i = 0; i < population_size; i++) {
        thread_local mt19937 local_rng(random_device{}());
        thread_local uniform_int_distribution<int> tau_dist(cfg.tau_min, cfg.tau_max);
        thread_local uniform_real_distribution<float> rho_dist(0.1, cfg.rho_max); //no inital neg rho allowed
        thread_local params p_local = p; // thread-local copy of params
        Individual ind{tau_dist(local_rng), rho_dist(local_rng), 0.0f};
        ind.fitness = evaluate(p_local, ind.tau, ind.rho, log_file);
        #pragma omp critical
        {
        population[i] = ind; 
        }
    }
    
    Individual best_all_time = population[0];
    vector<Individual> top4(4);
    for (int gen = 0; gen < generations; gen++) {
        //Step 2: get best 4
        sort(population.begin(), population.end(),[](const Individual& a, const Individual& b) {return a.fitness > b.fitness;}); //sorting before if blocks
        if ((gen % 4 == 0) && gen != 0) {
            exchange_top_individuals(population, rank); //outputs population 0->3 as the merged best
        };
        for (int ii = 0; ii < 4; ii++) {
            top4[ii] = population[ii];
        };
        if (top4[0].fitness > best_all_time.fitness) {
            best_all_time = top4[0];
        };

        //Step 3: new gen
        vector<Individual> next_gen(population_size);
        #pragma omp parallel for shared(next_gen, log_file) firstprivate(gen, top4) // 
        for (int i = 0; i < population_size; i++) {
            thread_local params p_local = p; // thread-local copy of params
            const Individual& parent = top4[i % 3]; //i % 3 gives 0, 1, 2, 0, 1, 2...
            const Individual& parent2 = top4[(i + 1) % 3];
            Individual child_cross = crossover(parent, parent2);
            Individual child = mutate(child_cross);
            child.fitness = evaluate(p_local, child.tau, child.rho, log_file);
            #pragma omp critical    
            {
            next_gen[i] = child;
            }
        }
        population = std::move(next_gen);    
    };//end generations loop

    //selecting LOCAL best after final step (since loop ends with generation of last new gen)
    sort(population.begin(), population.end(),[](const Individual& a, const Individual& b) {return a.fitness > b.fitness;}); //sorting before if blocks
    exchange_top_individuals(population, rank); //outputs population 0->3 as the merged best
    for (int ii = 0; ii < 4; ii++) {
        top4[ii] = population[ii];
    };
    if (top4[0].fitness > best_all_time.fitness) {
    best_all_time = top4[0];
    };
    
    //synchronizing final best
    //prepare local best_all_time buffer arrays for exchange
    float local_best_data[3] = {static_cast<float>(best_all_time.tau),best_all_time.rho,best_all_time.fitness};
    float peer_best_data[3];

    //exchange
    int peer = rank ^ 1;
    MPI_Sendrecv(local_best_data, 3, MPI_FLOAT, peer, 0,peer_best_data, 3, MPI_FLOAT, peer, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //peer best
    Individual peer_best = {static_cast<int>(peer_best_data[0]),peer_best_data[1],peer_best_data[2]};
    
    //final best
    if (peer_best.fitness > best_all_time.fitness) {
        best_all_time = peer_best;
    }

    log_file << "\n Final best after " << generations << " generations:\n";
    log_file << "Fitness = " << best_all_time.fitness << "  |  tau = " << best_all_time.tau << "  |  rho = " << best_all_time.rho << "\n";
    log_file.close();

    MPI_Finalize();

    //only return one real value to python wrapper
    if (rank == 0) {
        return {best_all_time.tau, best_all_time.rho};
    } else {    
        return {0, 0};
    }
}

#ifndef PROFILING
//make nb module
NB_MODULE(Genetic_algo, m) {
    m.def("Genetic_algo", &run_ga ,"genetic algo");
};
#endif