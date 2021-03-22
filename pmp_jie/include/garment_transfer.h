//
// Created by jie zou on 2021/3/21.
//

#ifndef DGP_JIE_PMP_JIE_INCLUDE_GARMENT_TRANSFER_H_
#define DGP_JIE_PMP_JIE_INCLUDE_GARMENT_TRANSFER_H_

#include "deformation_transfer.h"

namespace pmp_jie {

class GarmentTransfer: public pmp_jie::DeformationTransfer {
 public:
  GarmentTransfer();

  int PositionTransfer(const std::string &file, const std::vector<size_t> &fixed_vidx);
  int Fit(const std::string &filename);
  int WriteGarment(const std::string &filename_with_path) const;

 private:
  size_t SetConstraintTerm(const std::map<size_t, size_t> &fixed_points, std::vector<Eigen::Triplet<double>> &triplets) const;
  int EigenVerts(Eigen::MatrixXd &verts) const;
  Eigen::Vector3d GetVert(const Eigen::MatrixXd &verts, const size_t &vid) const;
  size_t SetInterpenetrationTermAndMatrix(const Eigen::MatrixXd &cloth_verts,
                                          std::vector<size_t> &interpenetrations,
                                          std::vector<Eigen::Triplet<double>> &triplets,
                                          Eigen::MatrixXd &mat);
  int SetFirstTermTriplets(const size_t &row,
                           const size_t &vid,
                           const Eigen::Vector3d &coeffs,
                           std::vector<Eigen::Triplet<double>> &triplets);
  size_t SetSmoothWarpingTerm(std::vector<Eigen::Triplet<double>> &triplets);
  int SetSmoothWarpingConstantsMatrix(const Eigen::MatrixXd &cloth_verts,
                                      Eigen::MatrixXd &mat);
  Eigen::Vector3d GetLocalVertexDeformation(const Eigen::MatrixXd &mesh_verts,
                                            const int &vidx);
  int SetTriplets(const int &vidx,
                        const int &adj_vidx,
                        const Eigen::Vector3d &coeffs,
                        std::vector<Eigen::Triplet<double>> &triplets);
  size_t SetDampingTerm(std::vector<Eigen::Triplet<double>> &triplets);
  int SetDampingConstantsMatrix(const Eigen::MatrixXd &cloth_verts,
                                Eigen::MatrixXd &mat);

  cinolib::Trimesh<> cloth_;
};

}
#endif //DGP_JIE_PMP_JIE_INCLUDE_GARMENT_TRANSFER_H_
