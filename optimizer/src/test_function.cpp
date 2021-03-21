//
// Created by jie zou on 2020/12/18.
//
#include "../include/test_function.h"

TestFunc::TestFunc(const double *t, const double *y, const size_t &n) {
  Eigen::Map<const Eigen::VectorXd> temp_t(t,n);
  t_ = temp_t;
  Eigen::Map<const Eigen::VectorXd> temp_y(y,n);
  y_ = temp_y;
  dim_of_x_ = 4;
  dim_of_f_ = n;
}
size_t TestFunc::dim_of_x() const {
  return dim_of_x_;
}
size_t TestFunc::dim_of_f() const {
  return dim_of_f_;
}
int TestFunc::jac(const double *x_ptr, Eigen::SparseMatrix<double> &JT) {
  //fi = (yi - {a*ti^3+b*ti^2+c^ti+d});
  //x  = [a,b,c,d];
  //JT_i = (-ti^3, -ti^2, -ti, 1.0);
  std::vector<Eigen::Triplet<double>> triplets;
  for(size_t i = 0; i < dim_of_f_; ++i) {
    triplets.emplace_back(0, i, -t_(i)*t_(i)*t_(i));
    triplets.emplace_back(1, i, -t_(i)*t_(i));
    triplets.emplace_back(2, i, -t_(i));
    triplets.emplace_back(3, i, -1.0);
  }
  JT.resize(dim_of_x_, dim_of_f_);
  JT.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::VectorXd w(dim_of_f_);
  w.setConstant(std::sqrt(2));
  JT = JT*w.asDiagonal();
  return 0;
}
int TestFunc::val(const double *x_ptr, double *f_ptr) {
  const Eigen::Map<const Eigen::VectorXd> curr_x(x_ptr, dim_of_x_);
  Eigen::Map<Eigen::VectorXd> f(f_ptr, dim_of_f_);
  for(size_t i = 0; i < dim_of_f_; ++i) {
    f(i) = y_(i)-(curr_x(0)*t_(i)*t_(i)*t_(i)+curr_x(1)*t_(i)*t_(i)+curr_x(2)*t_(i)+curr_x(3));
  }
  return 0;
}
double TestFunc::energy(const double *x_ptr) {
  double Fx(0.0);
  const Eigen::Map<const Eigen::VectorXd> curr_x(x_ptr, dim_of_x_);
  for(size_t i = 0; i < dim_of_f_; ++i) {
    const double temp = y_(i)-(curr_x(0)*t_(i)*t_(i)*t_(i)+curr_x(1)*t_(i)*t_(i)+curr_x(2)*t_(i)+curr_x(3));
    Fx += temp*temp;
  }
  return Fx;
}
