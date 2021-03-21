//
// Created by jie zou on 2020/12/14.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_FUNCTION_H_
#define DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_FUNCTION_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace opt {

class Function {
 public:

  virtual size_t dim_of_x() const = 0;
  virtual size_t dim_of_f() const = 0;

  virtual int jac(const double *x_ptr, Eigen::SparseMatrix<double> &JT) = 0;
  virtual int val(const double *x_ptr, double *f_ptr) = 0;
  virtual double energy(const double *x_ptr) = 0;

};

//Local interior point method energy/function
class LocalFunction: public Function {
 public:
  virtual void set_vid(const size_t &vid) = 0;
};

class GradMethod {
 public:
  virtual size_t dim_of_x() const = 0;
  virtual void grad(const double *x_ptr, double *g) = 0;
  virtual double energy(const double *x_ptr) = 0;

};
//Interior point method local energy/function
class IPMLocalFunction: public GradMethod  {
 public:
  virtual void set_vid(const size_t &vid) = 0;
  virtual bool is_interior_point(const double *x_ptr) = 0;

};

}


#endif //DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_FUNCTION_H_
