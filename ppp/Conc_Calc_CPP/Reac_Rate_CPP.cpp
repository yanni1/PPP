#include "Params.h"
#include "Reac_Rate_CPP.h"
using namespace std;

float reaction_rate_cpp(vector<vector<vector<float> > > C, int j, int i){
    //init params
    params p = params(100);
    float k_b = (p.k0 * 0.25) * C[1][j][i];
    float k_f = (p.k0) * C[0][j][i] * C[2][j][i];
    float k_consumption = k_f - k_b;
    return k_consumption;
}