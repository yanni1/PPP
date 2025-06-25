#include <vector>
#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"


int main(int, char**)
{
    int nt_given = 10000;
    params p = params(nt_given);
    cout << "nt set to: " << nt_given << '\n' << endl;
    vector<float> Ct(nt_given * p.ns * p.ny * p.nx, 0.0f); //flattened vector values set to 0f (index = (s * ny + y) * nx + x)
    vector<float> C(p.ns * p.ny * p.nx, 0.0f); //flattened input buffer
    vector<float> Cn(p.ns * p.ny * p.nx, 0.0f); //flattened output buffer
    fill(C.begin() + p.idx(0, 0, 0), C.begin() + p.idx(0, 0, 0) + p.nx, p.CO_reservoir);

     for (int n = 0; n < nt_given; n++) {
        cout << "timestep:" << n << '\n' << endl;
        CC_CPP(p, C, Cn);
        #pragma omp parallel for simd
        for (int idx = 0; idx < p.ns * p.ny * p.nx; idx++) {
            Ct[n * p.ns * p.ny * p.nx + idx] = C[idx];
        }
        swap(C, Cn);
    };                
    return 0;
}