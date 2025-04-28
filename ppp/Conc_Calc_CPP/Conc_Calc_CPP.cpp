/*
 *  C++ source file for module ppp.Conc_Calc_CPP
 */


// See http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html for examples on how to use pybind11.
// The example below is modified after http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#More-on-working-with-numpy-arrays
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <vector>
#include <omp.h>
#include <iostream>
namespace nb = nanobind;


#include <iostream>
#include <omp.h>

#include <iostream>
#include <omp.h>
#include <nanobind/nanobind.h>

namespace nb = nanobind;

// A simple function that does parallel addition using OpenMP
int parallel_add(int a, int b) {
    int result = 0;
    
    // Parallel for loop
    #pragma omp parallel
    {
        // Parallel threads will sum the values and store in 'result'
        #pragma omp single
        {
            result = a + b;
            int thread_id = omp_get_thread_num();
            int total_threads = omp_get_num_threads();
            std::cout << "Thread " << thread_id << " of " << total_threads << " is working!" << std::endl;
        }
    }
    
    return result;
}

NB_MODULE(Conc_Calc_CPP, m) {
    m.def("parallel_add", &parallel_add, "A function that adds two numbers in parallel using OpenMP");
}
