//
// Created by jie zou on 2020/7/20.
//
#include "../include/linear_solver.h"
int opt::LeastSquaresSolve(const Eigen::SparseMatrix<double> &A,
                           const Eigen::MatrixXd &b,
                           Eigen::MatrixXd &results) {
  results.setZero(A.cols(), b.cols());
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A.transpose() * A);
  if(solver.info()!= Eigen::Success) {
    std::cerr << "LeastSquaresSolve cerr: matrix decomposition failed!" << std::endl;
    return 0;
  }
  results = solver.solve(A.transpose() * b);
  return 0;
}
int opt::LinearBlockSolve(const Eigen::SparseMatrix<double> &invL,
                              const Eigen::SparseMatrix<double> &Jt,
                              const Eigen::SparseMatrix<double> &star,
                              const Eigen::VectorXd &b1,
                              const Eigen::VectorXd &b2,
                              Eigen::VectorXd &x,
                              Eigen::VectorXd &lambda) {
//  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> result2_solver;
//  Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> result2_solver;
//  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> result2_solver;
//  Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> result2_solver;
//  const auto eivals = Eigen::MatrixXd(Jt.transpose() * invL * Jt).eigenvalues();
//  std::cout << "The eigenvalues of the 3x3 matrix of ones are: " << std::endl << eivals.real() << std::endl;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> result2_solver;
  result2_solver.compute(Jt.transpose() * invL * Jt - star);
  if(result2_solver.info()!= Eigen::Success) {
    std::cerr << "LinearBlockSolve cerr: matrix decomposition failed!" << std::endl;
    return 0;
  }
  lambda = result2_solver.solve(Jt.transpose() * invL * b1 - b2);
  x = invL * (b1 - Jt * lambda);
  return 0;
}
opt::LinearSolver::LinearSolver()=default;
int opt::LinearSolver::SetMatrix(const std::vector<Eigen::Triplet<double>> &triplets,
                                 const size_t &rows,
                                 const size_t &cols,
                                 const double &coeff) {
  if(A_.rows()==0||A_.cols()==0) {
    A_.resize(cols, cols);
    A_.setZero();
  }
  Eigen::SparseMatrix<double> temp_mat(rows, cols);
  temp_mat.setFromTriplets(triplets.begin(), triplets.end());
  A_ = A_ + coeff * temp_mat.transpose() * temp_mat;
  return 0;
}
int opt::LinearSolver::SetMatrix(const std::vector<Eigen::Triplet<double>> &triplets,
                                 const size_t &rows,
                                 const size_t &cols,
                                 const double &coeff,
                                 const Eigen::MatrixXd &constant_mat) {
  if(A_.rows()==0||A_.cols()==0) {
    A_.resize(cols, cols);
    A_.setZero();
  }
  if(B_.rows()==0||B_.cols()==0) {
    B_.resize(cols, constant_mat.cols());
    B_.setZero();
  }
  Eigen::SparseMatrix<double> temp_mat_A(rows, cols);
  temp_mat_A.setFromTriplets(triplets.begin(), triplets.end());
  A_ = A_ + coeff * temp_mat_A.transpose() * temp_mat_A;
  B_ = B_ + coeff * temp_mat_A.transpose() * constant_mat;
  return 0;
}
int opt::LinearSolver::Solve(const size_t &rows, const size_t &cols, Eigen::MatrixXd &result_mat) {
  assert(A_.cols() == rows && B_.cols() == cols);
  assert(A_.size() != 0 && B_.size() != 0);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_);
  if(solver.info()!= Eigen::Success) {
    std::cerr << "LINEARSOLVER FUNCTION Solve: matrix decomposition failed!" <<std::endl;
    return 0;
  }
  result_mat.resize(rows, cols);
  result_mat = solver.solve(B_);
  return 0;
}
