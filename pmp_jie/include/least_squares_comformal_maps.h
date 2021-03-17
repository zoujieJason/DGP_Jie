//
// Created by jie zou on 2021/3/16.
//

#ifndef DGP_JIE_PMP_JIE_INCLUDE_LEAST_SQUARES_COMFORMAL_MAPS_H_
#define DGP_JIE_PMP_JIE_INCLUDE_LEAST_SQUARES_COMFORMAL_MAPS_H_

#include <cinolib/meshes/meshes.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "tools.h"

namespace pmp_jie {

Eigen::MatrixXd LSCM(const cinolib::Trimesh<> &trimesh);

}

#endif //DGP_JIE_PMP_JIE_INCLUDE_LEAST_SQUARES_COMFORMAL_MAPS_H_
