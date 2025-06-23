#include "Params.h"
#include "Reac_Rate_CPP.h"
using namespace std;

float reaction_rate_cpp(const vector<float> &C, int j, int i, const params& p){ 
    float k_b = (p.k0 * 0.25) * C[p.idx(1,j,i)];
    float k_f = (p.k0) * C[p.idx(0,j,i)] * C[p.idx(2,j,i)];
    float k_consumption = k_f - k_b;
    return k_consumption;
}