//
// Created by jie zou on 2020/8/11.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_DEFORMATION_TRANSFER_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_DEFORMATION_TRANSFER_H_

#include "correspond.h"

namespace pmp_jie {

class DeformationTransfer: public pmp_jie::Correspond {
 public:
  DeformationTransfer();

  int Transfer(const std::string &file);

  int WriteDeformedResult(const std::string &filename_with_path="") const;

 protected:

  int ReadDeformedTarget(const std::string &file);

  int SetAffineMatrix(Eigen::MatrixXd &affine_matrix) const;
  Eigen::Matrix3d GetPolyAffineMatrix(const size_t &pid) const;

  cinolib::Trimesh<> deformed_target_;
  Eigen::MatrixXd result_;

};

}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_DEFORMATION_TRANSFER_H_
