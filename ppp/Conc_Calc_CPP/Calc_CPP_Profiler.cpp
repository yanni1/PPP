#include <vector>
//#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"


int main(int, char**)
{
    int nt_given = 1;
    params p = params(nt_given);
    cout << "nt set to: " << nt_given << '\n' << endl;
    vector<float> Ct(nt_given * p.ns * p.ny * p.nx, 0.0f); //flattened vector values set to 0f (index = (s * ny + y) * nx + x)
    vector<float> C(p.ns * p.ny * p.nx, 0.0f); //flattened input buffer
    vector<float> Cn(p.ns * p.ny * p.nx, 0.0f); //flattened output buffer
    fill(C.begin() + p.idx(0, 0, 0), C.begin() + p.idx(0, 0, 0) + p.nx, p.CO_reservoir);

    for (int n = 0; n < nt_given; n++) {
        cout << "timestep:" << n << '\n' << endl;
        CC_CPP(p, C, Cn); 
        for (int s = 0; s < p.ns; ++s) {
            for (int y = 0; y < p.ny; ++y) {
                for (int x = 0; x < p.nx; ++x) {
                    int ct_idx = ((n * p.ns + s) * p.ny + y) * p.nx + x;
                    int c_idx = (s * p.ny + y) * p.nx + x;
                    Ct[ct_idx] = C[c_idx];
                }
            }
        }
        swap(C, Cn); // reuse buffers for next step without copying => C@n -> C@n+1 (input) and Cn now holds old state (C@n) => does not matter, gets replaced by computed data anyway               
    }    

    return 0;
}