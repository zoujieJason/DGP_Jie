//
// Created by jie zou on 2020/12/18.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_ANGLE_BASED_OPTIMIZATION_XU_H_
#define DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_ANGLE_BASED_OPTIMIZATION_XU_H_

#include "polygon_mesh.h"
#include "function.h"
#include "tools.h"

#include "nonlinear_solver.h"

namespace pmp_jie {

class Xu06: public opt::LocalFunction {
 public:
  Xu06(pmp_jie::Base *base_mesh, const size_t &vid);

  size_t dim_of_x() const override;
  size_t dim_of_f() const override;

  int jac(const double *x_ptr, Eigen::SparseMatrix<double> &JT) override;
  double energy(const double *x_ptr) override;
  int val(const double *x_ptr, double *f_ptr) override;
  void set_vid(const size_t &vid) override;

 private:
  Eigen::Vector3d get_bisector_line(const size_t &vid, const size_t &adj_vid) const;
  Eigen::Vector3d angle_bisector_line(const size_t &vid, const size_t &adj_vid) const;


  pmp_jie::CPMesh *mesh_ptr_;
  size_t vid_{0};
  size_t dim_of_x_{0};
  size_t dim_of_f_{0};

};

double LM_optimize_vert_Xu06(pmp_jie::Base *mesh_ptr, const size_t &vid);

}

#endif //DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_ANGLE_BASED_OPTIMIZATION_XU_H_
