/*
 *  C++ source file for module ppp.Conc_Calc_CPP
 */


// See http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html for examples on how to use pybind11.
// The example below is modified after http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#More-on-working-with-numpy-arrays

#ifndef PROFILING
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <nanobind/stl/vector.h>
namespace nb = nanobind;
#endif

#include <vector>
#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Reac_Rate_CPP.h"

using namespace std;

void CC_CPP(const params& p, vector<vector<vector<float> > >& C, vector<vector<vector<float> > >& Cn) { //make the passed vectors pointers (buffers)

    
    //deep copy of C for new time step
    //vector<vector<vector<float> > > Cn(C);
    // instead of this Cn is modified directly so no large vectors needs to be returned => work directly into the two preallocated buffers
    //start conc update
    for (int j = 1; j < p.ny-1; j++){
        for (int i = 1; i < p.nx-1; i++){
            // circular O2 vent
            int x_dist2 = (i-p.a) * (i-p.a);
            int y_dist2 = (j-p.b) * (j-p.b);
            int distSq = x_dist2 + y_dist2;
            if (distSq <= p.r2){
                C[2][j][i] = p.O2_reservoir; // write o2 vent in C_O2
            };

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
            float advec_O2 = p.v[i] * dy_O2;

            //reaction
            float k = reaction_rate_cpp(C, j, i, p);
            float reac_CO = k * C[0][j][i];
            float reac_CO2 = -1 * float(reac_CO); //right side of reaction
            float reac_O2 = float(reac_CO);
            
            //update conc
            Cn[0][j][i] = C[0][j][i] + (p.dt * (Diff_CO - advec_CO - reac_CO));
            Cn[1][j][i] = C[1][j][i] + (p.dt * (Diff_CO2 - advec_CO2 - reac_CO2));
            Cn[2][j][i] = C[2][j][i] + (p.dt * (Diff_O2 - advec_O2 - reac_O2));

        };
    };
    //edges
    for (int j = 0; j < p.ny; j++) {
        Cn[0][j][0] = Cn[0][j][1]; //setting vertical edge (left) conc to the same as column next to it
        Cn[1][j][0] = Cn[1][j][1];
        Cn[2][j][0] = Cn[2][j][1];
        Cn[0][j][p.nx - 1] = Cn[0][j][p.nx - 2]; //setting vertical edge (right) conc to the same as column next to it
        Cn[1][j][p.nx - 1] = Cn[1][j][p.nx - 2];
        Cn[2][j][p.nx - 1] = Cn[2][j][p.nx - 2];
        //re set CO boundary condition
        fill(C[0][0].begin(), C[0][0].end(), p.CO_reservoir);
    }
    return;
}

#ifndef PROFILING
//make nb module
NB_MODULE(Conc_Calc_CPP, m) {
    m.def("CC_CPP", &CC_CPP, "updates conc each timestep");
};
#endif
