#include <cusp/array1d.h>
#include <cusp/coo_matrix.h>
#include <cusp/detail/timer.h>
#include <cusp/ell_matrix.h>
#include <cusp/krylov/cg.h>
#include <cusp/monitor.h>
#include <cusp/multiply.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/transpose.h>

typedef cusp::host_memory HOST;
typedef cusp::device_memory DEVICE;

// Solve Ax = b
bool solve(cusp::coo_matrix<int,float,DEVICE> &A, cusp::array1d<float,DEVICE> &x, cusp::array1d<float,DEVICE> &b);
