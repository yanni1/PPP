// Params.cpp
#include <iostream>
#include "Params.h"
using namespace std;

params::params(int nt, int tau, float rho) {
    params::dx = 0.01;
    params::dy = 0.01;
    params::dt = 0.0005;
    params::Lx = 0.5;
    params::Ly = 2.0;
    params::D = 0.002;
    params::k0 = 1000.0;
    params::vmax = 0.5;
    params::CO_reservoir = 1.2;
    params::O2_reservoir = 1.0;
    params::nx = params::nxf();
    params::ny = params::nyf();
    params::nxny = nx * ny;
    params::ns = 3;
    params::nt = nt;
    params::x = params::xf();
    params::y = params::yf();
    params::b = params::bf();
    params::a = params::af();
    params::r = params::rf();
    params::r2 = params::r2f();
    params::v = params::vf();
    params::alpha = 0.9;
    params::mu_CO = 0.4; 
    params::mu_CO2 = 0.6;
    params::sigma_CO = 0.5; //0.5
    params::sigma_CO2 = 0.5;
    params::x0_cc = params::x0_ccf();
    params::y0_cc = params::y0_ccf();
    params::semiMaj = params::semiMajf();
    params::semiMin = params::semiMinf();
    params::tau = tau;
    params::rho = rho;


};
//actual mem functions
int params::nxf() { return(int(Lx/dx)); };
int params::nyf() { return(int(Ly/dy)); };
int params::x0_ccf() { return(int(nx/2)); };
int params::y0_ccf() { return(int((Ly - 0.1)/dy)); };
int params::semiMajf() { return(int(nx/2)); }; //25
int params::semiMinf() { return(int(0.1/dy)); };  //10

vector<float> params::xf() { 
    vector<float> x(nx+1); 
    for (int i = 0; i < nx+1; i++){
        float ddx = float(i) * float(dx);
        x[i] = ddx;
    };
    return x;
};

vector<float> params::yf() {
    vector<float> y(ny+1); 
    for (int i = 0; i < ny+1; i++){
        float ddy = float(i) * float(dy);
        y[i] = ddy;
    };
    return y;
};

vector<float> params::vf() {
    vector<float> v(nx + 1); // new local vector
    for (int i = 0; i <= nx; i++) {
        float xx = (x[i] - (Lx * 0.5)) / Lx;
        float xx2 = xx * xx;
        v[i] = (1 - 0.5f * xx2);
    }
    return v;
}
/*
vector<float> params::vf() {
    for (int i = 0; i < nx+1; i++) {
        float xx = (x[i] - (Lx * 0.5)) / Lx;
        float xx2 = xx * xx;
        float xfin = (1 - 0.5*xx2);
        v.push_back(xfin);
    };
    return v;
}
    remove v.push_back and init the vectro with right size beforehand and then allocate values
    */

// o2 vent circle parameters
int params::bf() { return (ny*0.1); }; //center of circle y
int params::af() { return (nx*0.25); }; //center of circle x
float params::rf() { return (float(nx)*float(0.125)); }; //radius of circle
float params::r2f() { return (r * r); };

