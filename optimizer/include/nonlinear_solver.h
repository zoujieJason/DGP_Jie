//
// Created by jie zou on 2020/12/14.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_NONLINEAR_SOLVER_H_
#define DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_NONLINEAR_SOLVER_H_

#include <iostream>
#include <vector>

#include "function.h"
#include "fast_ATA.h"

namespace opt {

class Base_Optimizer {
 public:
  virtual int solve(double *x_ptr, Function &func, int &num_iter) = 0;
 protected:
  virtual int compute_grad(const double *x_ptr, Function &func, double *f_ptr) = 0;
  virtual int solve_step(double *h_ptr, const double &mu) = 0;
  Eigen::SparseMatrix<double> JT_,JTJ_;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> slv_;
  Eigen::VectorXd g_, h_;
};

class LM_Gauss_Newton: public opt::Base_Optimizer {
 public:
  explicit LM_Gauss_Newton(const Function &func);
  int solve(double *x_ptr, Function &func, int &num_iter) override;
 protected:
  int compute_grad(const double *x_ptr, Function &func, double *f_ptr) override;
  int get_diag_ptr();
  int solve_step(double *h_ptr, const double &mu) override;

  std::vector<double *> diag_ptr_;
};

class More_LM_Gauss_Newton: public opt::LM_Gauss_Newton {
 public:
  explicit More_LM_Gauss_Newton(const Function &func);
  int solve(double *x_ptr, Function &func, int &num_iter) override;
 protected:
  int solve_step(double *h_ptr, const double &mu) override;
  void set_JTJ_diag(const double *diag_ptr);
  void compute_diag(double *diag);

  Eigen::VectorXd diag_;
};

class Dog_Leg: public opt::Base_Optimizer {
 public:
  explicit Dog_Leg(const Function &func);
  int solve(double *x_ptr, Function &func, int &num_iter) override;
 protected:
  int compute_grad(const double *x_ptr, Function &func, double *f_ptr) override;
  int solve_step(double *h_ptr, const double &mu) override;
  double dog_leg_step(const double &alpha,
                      const double *h_gn_ptr,
                      double *h_ptr,
                      const double &delta,
                      const double &Fx);
};

}

#endif //DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_NONLINEAR_SOLVER_H_
