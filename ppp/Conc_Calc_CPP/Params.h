// Params.h

#ifndef Conc_Calc_CPP_PARAMS_H
#define Conc_Calc_CPP_PARAMS_H

#include <iostream>
#include <vector>
using namespace std;

class params {
public:
    float dx; //m
    float dy; //m
    float dt; //s
    float Lx; //m
    float Ly; //m
    float D; //m^2/s
    float k0; //s^-1
    float vmax; //m/s
    float CO_reservoir; //mol/m^3
    float O2_reservoir; //mol/m^3
    int nx; //number of gridpoints in x
    int ny; //number of gridpoints in y
    int nt; //number of time steps
    vector<float> x;
    vector<float> y;
    int b;
    int a;
    float r;
    float r2;
    vector<float> v;
    //carbon capture params
    float alpha;
    float mu_CO;
    float mu_CO2;
    float sigma_CO;
    float sigma_CO2;
    int x0_cc;
    int y0_cc;
    int semiMaj;
    int semiMin;
    //precomputed const fot flat indexing => inline fct
    int ns;
    int nxny = 0;
    //GA params
    int tau;
    float rho;

    // Constructor (defined in cpp file)
    params(int nt);

    //setters
    void setTau(int new_tau);
    void setRho(float new_rho);
    
    // member function declarations (moved to .cpp) due to compiler issues
    int nxf();
    int nyf();
    vector<float> xf();
    vector<float> yf();
    int bf();
    int af();
    float rf();
    float r2f();
    vector<float> vf();
    int x0_ccf();
    int y0_ccf();
    int semiMajf();
    int semiMinf();
    //inline idx calculator
    inline int idx(int s, int j, int i) const {
        return s * nxny + j * nx + i;
    }
};

#endif //Conc_Calc_CPP_PARAMS_H
