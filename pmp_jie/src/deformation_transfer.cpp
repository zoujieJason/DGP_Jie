//
// Created by jie zou on 2020/8/11.
//
#include "deformation_transfer.h"
pmp_jie::DeformationTransfer::DeformationTransfer()=default;
int pmp_jie::DeformationTransfer::Transfer(const std::string &file) {
  ReadDeformedTarget(file);

  const auto cols = num_polys_ + num_verts_;
  std::vector<Eigen::Triplet<double>> idterm_triplets;
  const auto idterm_rows = SetIdentityTerm(idterm_triplets);
  Eigen::MatrixXd affine_matrix;
  SetAffineMatrix(affine_matrix);

  opt::LinearSolver solver;
  solver.SetMatrix(idterm_triplets,idterm_rows,cols,1.0,affine_matrix);
  Eigen::MatrixXd  result;
  solver.Solve(cols, 3, result);
  WriteDeformedTrimesh(result);

  return 0;
}
int pmp_jie::DeformationTransfer::ReadDeformedTarget(const std::string &file) {
  deformed_target_.clear();
  deformed_target_.load(file.c_str());
  return 0;
}
int pmp_jie::DeformationTransfer::SetAffineMatrix(Eigen::MatrixXd &affine_matrix) const {
  affine_matrix.resize(3*p2p_.size(),3);
  int pid = 0;
  for(const auto &pair: p2p_) {
    const auto target_pidx        = pair.second;
    const auto poly_affine_matrix = GetPolyAffineMatrix(target_pidx);
    for(int i = 0; i < 3; ++i) {
      affine_matrix.row(3*pid + i) = poly_affine_matrix.col(i);
    }
    ++pid;
  }
  return 0;
}
Eigen::Matrix3d pmp_jie::DeformationTransfer::GetPolyAffineMatrix(const size_t &pid) const {
  return pmp_jie::GetPolyGradientMatrix(target_, pid)*pmp_jie::GetPolyGradientMatrix(deformed_target_, pid, false);
}
int pmp_jie::DeformationTransfer::WriteDeformedTrimesh(const Eigen::MatrixXd &result) const {
  const auto filename = GetLoadFilename();
  const std::string ofile = pmp_jie::GetFilePath(filename)+"deformed_"+pmp_jie::GetFileName(filename, true);
  std::ofstream file(ofile, std::ofstream::out);
  if(!file.is_open()) {
    std::cout << "OPEN \"" <<  ofile << "\" FAILED." << std::endl;
    return __LINE__;
  }

  for(size_t vid = 0; vid < trimesh_.num_verts(); ++vid) {
    file << "v " << result(vid,0) << " " << result(vid,1) << " " << result(vid,2) << "\n";
  }
  for(size_t pid = 0; pid < trimesh_.num_polys(); ++pid) {
    file << "f ";
    for(const auto &vid: trimesh_.adj_p2v(pid)) {
      file << vid + 1 << " ";
    }
    file << "\n";
  }
  file.close();
  return 0;
}
