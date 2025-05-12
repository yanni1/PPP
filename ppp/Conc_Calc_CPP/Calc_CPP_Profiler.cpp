#include <gperftools/profiler.h>

#include <vector>
//#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"


int main(int, char**)
{
    ProfilerStart("Calc_CPP_profile.log");
    int nt_given = 1;
    params p = params(nt_given);
    cout << "nt set to: " << nt_given << '\n' << endl;
    vector<vector<vector<vector<float>>>> Ct(nt_given, vector<vector<vector<float>>>(3, vector<vector<float>>(p.ny, vector<float>(p.nx, 0.0)))); // Ct will contain C_CO C_CO2 and C_O2 => C[time, species,[y,[x,conc_value]]]
    vector<vector<vector<float> > > C(3, vector<vector<float> > (p.ny, vector<float>(p.nx, 0.0))); // C will contain C_CO C_CO2 and C_O2 => C[species,[y,[x,conc_value]]] 
    //C[0] = C_CO
    //C[1] = C_CO2
    //C[2] = C_O2
    // CO concetration at y = 0 => C0_reservoir
    fill(C[0][0].begin(), C[0][0].end(), p.CO_reservoir);

    for (int n = 0; n < nt_given; n++) {
        cout << "timestep:" << n << '\n' << endl;
        C = CC_CPP(nt_given, C);  // evolve concentration
        Ct[n] = C;                  // store current state
    };

    ProfilerStop();

    return 0;
}