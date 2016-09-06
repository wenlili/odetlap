#include "solver.h"

// Solve Ax = b
bool solve(cusp::coo_matrix<int,float,DEVICE> &A, cusp::array1d<float,DEVICE> &x, cusp::array1d<float,DEVICE> &b)
{
    // Compute AtA and Atb
    cusp::coo_matrix<int,float,DEVICE> At;
    cusp::transpose(A, At);
    cusp::ell_matrix<int,float,DEVICE> An;
    cusp::multiply(At, A, An);
    cusp::array1d<float,DEVICE> bn(x.size());
    cusp::multiply(At, b, bn);

    // Solve AtAx = Atb
    cusp::monitor<float> monitor(bn, 10000, 1e-6);
    cusp::precond::aggregation::smoothed_aggregation<int,float,DEVICE> M(An);
    cusp::krylov::cg(An, x, bn, monitor, M);
    return monitor.converged();
}
