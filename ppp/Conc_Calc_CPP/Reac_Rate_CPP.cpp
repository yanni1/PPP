#include "Params.h"
#include "Reac_Rate_CPP.h"
using namespace std;

float reaction_rate_cpp(const vector<vector<vector<float> > > &C, int j, int i, const params& p){ 
    /*You're passing C by value, which creates a full deep copy of the entire 3D vector C every time this function is called.

If C is size [3][ny][nx] and you're calling this inside loops (e.g., over i, j), then you're duplicating potentially tens or hundreds of MB per call—which explains your massive memory usage (1.28 GiB).*/

    float k_b = (p.k0 * 0.25) * C[1][j][i];
    float k_f = (p.k0) * C[0][j][i] * C[2][j][i];
    float k_consumption = k_f - k_b;
    return k_consumption;
}