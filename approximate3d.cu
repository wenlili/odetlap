#include "approximator.h"

// Append a nonzero element to A
static void append(cusp::coo_matrix<int,float,HOST> &A, const int index, const int row, const int col, const float val)
{
    A.row_indices[index] = row;
    A.column_indices[index] = col;
    A.values[index] = val;
}

// Approximate data from known points
bool approximate(cusp::array1d<float,HOST> &h_x, const int nslcs, const int nrows, const int ncols, const cusp::array1d<float4,HOST> &known, const float R)
{
    const int size = nslcs*nrows*ncols; // number of points
    const int m = size + known.size();  // number of equations
    const int n = size;                 // number of variables
    const int nnz = 7*size - 2*(nslcs*nrows + nslcs*ncols + nrows*ncols) + known.size(); // number of nonzero elements

    // Build A
    cusp::coo_matrix<int,float,HOST> h_A(m, n, nnz);
    int index = 0;
    for (int s = 0; s < nslcs; s++)
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++) {
                int i = s*nrows*ncols + r*ncols + c;
                float sum = 0;
                if (s > 0)       { append(h_A, index++, i, i-nrows*ncols, -R); sum += R; }
                if (s < nslcs-1) { append(h_A, index++, i, i+nrows*ncols, -R); sum += R; }
                if (r > 0)       { append(h_A, index++, i, i-ncols, -R);       sum += R; }
                if (r < nrows-1) { append(h_A, index++, i, i+ncols, -R);       sum += R; }
                if (c > 0)       { append(h_A, index++, i, i-1, -R);           sum += R; }
                if (c < ncols-1) { append(h_A, index++, i, i+1, -R);           sum += R; }
                append(h_A, index++, i, i, sum);
            }
    for (int i = 0; i < known.size(); i++)
        append(h_A, index++, n+i, int(known[i].x)*nrows*ncols+int(known[i].y)*ncols+int(known[i].z), 1);
    cusp::coo_matrix<int,float,DEVICE> A(h_A);

    // Build b
    cusp::array1d<float,HOST> h_b(m, 0);
    for (int i = 0; i < known.size(); i++)
        h_b[n+i] = known[i].w;
    cusp::array1d<float,DEVICE> b(h_b);

    // Solve Ax = b
    cusp::array1d<float,DEVICE> x(n, 0);
    bool result = solve(A, x, b);
    h_x = x;
    return result;
}
