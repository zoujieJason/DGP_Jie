//
// Created by jie zou on 2020/7/20.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_LINEAR_SOLVER_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_LINEAR_SOLVER_H_

#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace opt {
class LinearSolver {
 public:
  LinearSolver();
  int SetMatrix(const std::vector<Eigen::Triplet<double>> &triplets,
                 const int &rows,
                 const int &cols,
                 const double &coeff);
  int SetMatrix(const std::vector<Eigen::Triplet<double>> &triplets,
                 const int &rows,
                 const int &cols,
                 const double &coeff,
                 const Eigen::MatrixXd &constant_mat);
  int Solve(const int &rows, const int &cols, Eigen::MatrixXd &result_mat);

 private:
  Eigen::SparseMatrix<double> A_;
  Eigen::SparseMatrix<double> B_;
};
int LeastSquaresSolve(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd &b, Eigen::MatrixXd &results);
int LinearBlockSolve(const Eigen::SparseMatrix<double> &invL,
                     const Eigen::SparseMatrix<double> &Jt,
                     const Eigen::SparseMatrix<double> &star,
                     const Eigen::VectorXd &b1,
                     const Eigen::VectorXd &b2,
                     Eigen::VectorXd &x,
                     Eigen::VectorXd &lambda);
}
#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_LINEAR_SOLVER_H_
