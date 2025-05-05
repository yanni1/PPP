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
#include "Params.h"

namespace nb = nanobind;
using namespace std;

vector<float> CC_CPP(int nt_def) {
    //init params
    params p = params(nt_def);
    cout << "nt set to: " << nt_def << '\n' << endl;

    //init C matrix (all values set to 0)
    vector<vector<vector<float> > > C(3, vector<vector<float> > (p.ny, vector<float>(p.nx, 0.0))); // C will contain C_CO C_CO2 and C_O2 => C[species,[y,[x,conc_value]]]
    //C[0] = C_CO
    //C[1] = C_CO2
    //C[2] = C_O2
    // boundary conditions

    // CO concetration at y = 0 => C0_reservoir
    fill(C[0][0].begin(), C[0][0].end(), p.CO_reservoir);
    
    //start conc update
    for (int j = 1; j < p.ny+1; j++){
        for (int i = 1; i < p.nx+1; i++){
            // circular O2 vent
            int x_dist2 = (i-p.a) * (i-p.a);
            int y_dist2 = (j-p.b) * (j-p.b);
            int distSq = x_dist2 + y_dist2;
            if (distSq <= p.r2){
                C[2][j][i] = p.O2_reservoir; // write o2 vent in C_O2
            };
            //deep copy of C for new time step
            vector<vector<vector<float> > > Cn(C);

            //diffusion
            float D2x_CO = float((C[0][j][i+1] + C[0][j][i-1] - 2*C[0][j][i]) / (p.dx * p.dx));
            float D2y_CO = float((C[0][j+1][i] + C[0][j-1][i] - 2*C[0][j][i]) / (p.dy * p.dy));
            float Diff_CO = p.D * (D2x_CO + D2y_CO);

            float D2x_CO2 = float((C[1][j][i+1] + C[1][j][i-1] - 2*C[1][j][i]) / (p.dx * p.dx));
            float D2y_CO2 = float((C[1][j+1][i] + C[1][j-1][i] - 2*C[1][j][i]) / (p.dy * p.dy));
            float Diff_CO2 = p.D * (D2x_CO2 + D2y_CO2);

            float D2x_O2 = float((C[2][j][i+1] + C[2][j][i-1] - 2*C[2][j][i]) / (p.dx * p.dx));
            float D2y_O2 = float((C[2][j+1][i] + C[2][j-1][i] - 2*C[2][j][i]) / (p.dy * p.dy));
            float Diff_O2 = p.D * (D2x_O2 + D2y_O2);

            //advection
            float dy_CO = float((C[0][j][i] - C[0][j-1][i]) / p.dy);
            float advec_CO = p.v[i] * dy_CO;

            float dy_CO2 = float((C[1][j][i] - C[1][j-1][i]) / p.dy);
            float advec_CO2 = p.v[i] * dy_CO2;

            float dy_O2 = float((C[2][j][i] - C[2][j-1][i]) / p.dy);
            float advec_O2 = p.v[i] * dy_CO2;

            //reaction


        };
    };
    return p.v;


}


//make nb module
NB_MODULE(Conc_Calc_CPP, m) {
    m.def("CC_CPP", &CC_CPP, "updates conc each timestep");
};
