//
// Created by jie zou on 2020/7/27.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_ANGLE_RECONSTRUCT_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_ANGLE_RECONSTRUCT_H_

#include <memory>
#include <cinolib/meshes/meshes.h>

#include "linear_solver.h"

namespace pmp_jie {
class AngleRecons {
 public:
  AngleRecons();

  int LoadMesh(const cinolib::Trimesh<> &trimesh);
  int To2DSpace(const Eigen::VectorXd &alpha, Eigen::MatrixXd &verts) const;

 private:
  size_t SetTriplets(const Eigen::VectorXd &alpha, std::vector<Eigen::Triplet<double>> &triplets) const;
  Eigen::Matrix2d GetRorateMatrix(const Eigen::VectorXd &alpha, const size_t &pid) const;
  Eigen::Vector3d GetPolyAngles(const Eigen::VectorXd &alpha, const size_t &pid) const;

  //data
  cinolib::Trimesh<> trimesh_;
  size_t num_bod_verts_{0};

};
}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_ANGLE_RECONSTRUCT_H_
