#pragma once //prevents duplicate calls
#include <vector>
using namespace std;

float reaction_rate_cpp(const vector<vector<vector<float> > > &C, int j, int i, const params& p);
