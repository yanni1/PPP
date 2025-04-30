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

    // Constructor (defined in cpp file)
    params(int nt);
    
    // member function declarations (moved to .cpp) due to compiler issues
    int nxf();
    int nyf();
    vector<float> xf();
    vector<float> yf();
    int bf();
    int af();
    float rf();
    float r2f();
};

#endif //Conc_Calc_CPP_PARAMS_H
