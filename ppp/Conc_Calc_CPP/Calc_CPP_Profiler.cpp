#include <vector>
//#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"


int main(int, char**)
{
    int nt_given = 1000;
    params p = params(nt_given);
    cout << "nt set to: " << nt_given << '\n' << endl;
    vector<vector<vector<vector<float>>>> Ct(nt_given, vector<vector<vector<float>>>(3, vector<vector<float>>(p.ny, vector<float>(p.nx, 0.0)))); // Ct will contain C_CO C_CO2 and C_O2 => C[time, species,[y,[x,conc_value]]]
    vector<vector<vector<float> > > C(3, vector<vector<float> > (p.ny, vector<float>(p.nx, 0.0))); // C will contain C_CO C_CO2 and C_O2 => C[species,[y,[x,conc_value]]] 
    vector<vector<vector<float> > > Cn(3, vector<vector<float> > (p.ny, vector<float>(p.nx, 0.0))); // output buffer => now only made once :)
    fill(C[0][0].begin(), C[0][0].end(), p.CO_reservoir);

    for (int n = 0; n < nt_given; n++) {
        cout << "timestep:" << n << '\n' << endl;
        CC_CPP(p, C, Cn);  // evolve concentration directly in preallocated vectors, only passed as pointers => using 2 buffers essentially
        Ct[n] = C;                  // store current timestep "snapshot" into Ct
        swap(C, Cn); // reuse memory
    }                

    return 0;
}