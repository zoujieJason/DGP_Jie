//
// Created by jie zou on 2021/3/16.
//

#ifndef DGP_JIE_PMP_JIE_INCLUDE_BARYCENTRIC_MAPPING_H_
#define DGP_JIE_PMP_JIE_INCLUDE_BARYCENTRIC_MAPPING_H_

#include <cinolib/meshes/meshes.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "tools.h"
#include "linear_optimizer.h"

namespace pmp_jie {
//Mapping surface to circle.
enum WEIGHT_TYPE {
  UNIFORM,
  HARMONIC,
  COTLAPLACE
};
Eigen::MatrixXd BaryMapping(const cinolib::Trimesh<> &trimesh, const WEIGHT_TYPE &w_type, const Eigen::MatrixXd &bnd_uv);

Eigen::MatrixXd CircleBoundary(const size_t &n, const double &scale= 1.0);
}


#endif //DGP_JIE_PMP_JIE_INCLUDE_BARYCENTRIC_MAPPING_H_
