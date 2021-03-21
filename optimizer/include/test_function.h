//
// Created by jie zou on 2020/12/18.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_TEST_FUNCTION_H_
#define DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_TEST_FUNCTION_H_

#include "function.h"

class TestFunc: public opt::Function {
 public:
  TestFunc(const double *t, const double *y, const size_t &n);

  size_t dim_of_x() const override;
  size_t dim_of_f() const override;

  int jac(const double *x_ptr, Eigen::SparseMatrix<double> &JT) override;
  int val(const double *x_ptr, double *f_ptr) override;
  double energy(const double *x_ptr) override;

 private:
  Eigen::VectorXd t_;
  Eigen::VectorXd y_;

  size_t dim_of_x_{0};
  size_t dim_of_f_{0};

};

#endif //DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_TEST_FUNCTION_H_
