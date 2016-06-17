#ifndef SPARSE_HPP
#define SPARSE_HPP

#include "Eigen/Sparse"
#include <vector>
using namespace Eigen;
using namespace std;

// helper class for incremental construction of sparse matrices
class SpMatrix {
public:
    int m, n;
    vector< Triplet<double> > triplets;
    SpMatrix(int m, int n): m(m), n(n) {}
    void add(int i, int j, double value);
    void addBlock(int i, int j, int p, int q, MatrixXd block);
    VectorXd solve(const VectorXd &b);
};

#endif
