//
// Created by jie zou on 2020/7/27.
//
#include "../include/angle_reconstruct.h"
pmp_jie::AngleRecons::AngleRecons()=default;
int pmp_jie::AngleRecons::LoadMesh(const cinolib::Trimesh<> &trimesh) {
  trimesh_ = trimesh;
  num_bod_verts_ = trimesh_.get_boundary_vertices().size();
  return 0;
}
size_t pmp_jie::AngleRecons::SetTriplets(const Eigen::VectorXd &alpha, std::vector<Eigen::Triplet<double>> &triplets) const {
  const auto I = Eigen::Matrix2d::Identity();
  for(size_t pid = 0; pid < trimesh_.num_polys(); ++pid) {
    const auto R = GetRorateMatrix(alpha, pid);
    Eigen::MatrixXd coeff_matrix;
    coeff_matrix.resize(2, 6);
    coeff_matrix.leftCols(2)      = I-R;
    coeff_matrix.middleCols(2, 2) = -I;
    coeff_matrix.rightCols(2)     = R;
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      const auto curr_coeffs_matrix = coeff_matrix.middleCols(2 * trimesh_.poly_vert_offset(pid, adj_vid),2);
      //x, y
      triplets.emplace_back(2 * pid + 0, 2 * adj_vid + 0, curr_coeffs_matrix(0, 0));
      triplets.emplace_back(2 * pid + 0, 2 * adj_vid + 1, curr_coeffs_matrix(0, 1));
      //x ,y
      triplets.emplace_back(2 * pid + 1, 2 * adj_vid + 0, curr_coeffs_matrix(1, 0));
      triplets.emplace_back(2 * pid + 1, 2 * adj_vid + 1, curr_coeffs_matrix(1, 1));
    }
  }
  return 2 * trimesh_.num_polys();
}
Eigen::Matrix2d pmp_jie::AngleRecons::GetRorateMatrix(const Eigen::VectorXd &alpha, const size_t &pid) const {
  const auto   angles    = GetPolyAngles(alpha, pid);
  const double sin_alpha = std::sin(angles(0));
  const double sin_beta  = std::sin(angles(1));
  const double sin_gamma = std::sin(angles(2));
  const double cos_alpha = std::cos(angles(0));

  const double length    = sin_gamma / sin_beta;
  Eigen::Matrix2d R;
  R << cos_alpha, sin_alpha, -1.0 * sin_alpha, cos_alpha;
  R = R * length;
  return R;
}
Eigen::Vector3d pmp_jie::AngleRecons::GetPolyAngles(const Eigen::VectorXd &alpha, const size_t &pid) const {
  return {alpha(3 * pid + 0), alpha(3 * pid + 1) ,alpha(3 * pid + 2)};
}
int pmp_jie::AngleRecons::To2DSpace(const Eigen::VectorXd &alpha, Eigen::MatrixXd &verts) const {
  std::vector<Eigen::Triplet<double>> triplets;
  SetTriplets(alpha, triplets);
  const Eigen::Map<const Eigen::Vector3d> v0(trimesh_.vert(0).ptr());
  const Eigen::Map<const Eigen::Vector3d> v1(trimesh_.vert(1).ptr());
  Eigen::Vector4d fixed_points;
  fixed_points.setZero();
  fixed_points(0) = v0(0);
  fixed_points(1) = v0(1);
  fixed_points(2) = v1(0);
  fixed_points(3) = v1(1);

  const auto rows = 2*trimesh_.num_polys();
  const auto cols = 2*trimesh_.num_verts();
  Eigen::MatrixXd res;
  Eigen::SparseMatrix<double> M(rows,cols);
  M.setFromTriplets(triplets.begin(), triplets.end());
  const auto A = M.block(0,4,rows,cols-4);
  const auto B = M.block(0,0,rows,4);
  opt::LeastSquaresSolve(A, -B*fixed_points, res);


  verts.setZero(trimesh_.num_verts(), 2);
  verts.row(0) = Eigen::Vector2d(v0(0), v0(1));
  verts.row(1) = Eigen::Vector2d(v1(0), v1(1));
  for(size_t i = 2; i < trimesh_.num_verts(); ++i) {
    verts.row(i) = Eigen::Vector2d(res(2*(i-2)+0), res(2*(i-2)+1));
  }
  return 0;
}