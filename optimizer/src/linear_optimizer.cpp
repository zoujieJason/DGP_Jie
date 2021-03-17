//
// Created by jie zou on 2021/3/16.
//
#include "../include/linear_optimizer.h"
int opt::LeastSquaresSolveATA(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd &b, Eigen::MatrixXd &results) {
  results.resize(A.cols(), b.cols());
  results.setZero();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A.transpose() * A);
  if(solver.info()!= Eigen::Success) {
    std::cerr << "pmp_jie::least_squares_solve_ATA cerr: matrix decomposition failed!" << std::endl;
    return 0;
  }
  results = solver.solve(A.transpose() * b);
  return 1;
}
