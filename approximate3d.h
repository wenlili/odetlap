#include "solver.h"

// Approximate data from known points
bool approximate(cusp::array1d<float,HOST> &h_x, const int nslcs, const int nrows, const int ncols, const cusp::array1d<float4,HOST> &known, const float R);
