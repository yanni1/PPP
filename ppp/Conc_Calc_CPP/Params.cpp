// Params.cpp
#include <iostream>
#include "Params.h"
using namespace std;

params::params(int nt) {
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
    params::nt = nt;
    params::x = params::xf();
    params::y = params::yf();
    params::b = params::bf();
    params::a = params::af();
    params::r = params::rf();
    params::r2 = params::r2f();
    params::v = params::vf();
};
//actual mem functions
int params::nxf() { return(int(Lx/dx)); };
int params::nyf() { return(int(Ly/dy)); };

vector<float> params::xf() { 
    for (int i = 0; i < nx+1; i++){
        float ddx = float(i) * float(dx);
        x.push_back(ddx);
    };
    return x;
};

vector<float> params::yf() { 
    for (int i = 0; i < ny+1; i++){
        float ddy = float(i) * float(dy);
        y.push_back(ddy);
    };
    return y;
};

vector<float> params::vf() {
    for (int i = 0; i < nx+1; i++) {
        float xx = (x[i] - (Lx * 0.5)) / Lx;
        v.push_back(xx);
    };
    return v;
}

// o2 vent circle parameters
int params::bf() { return (ny*0.1); }; //center of circle y
int params::af() { return (nx*0.25); }; //center of circle x
float params::rf() { return (float(nx)*float(0.125)); }; //radius of circle
float params::r2f() { return (r * r); };

