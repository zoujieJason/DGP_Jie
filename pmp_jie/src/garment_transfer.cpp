//
// Created by jie zou on 2021/3/21.
//
#include "../include/garment_transfer.h"
pmp_jie::GarmentTransfer::GarmentTransfer()=default;
int pmp_jie::GarmentTransfer::PositionTransfer(const std::string &file, const std::vector<size_t> &fixed_vidx) {
  ReadDeformedTarget(file);
  const auto cols = num_polys_ + num_verts_;
  std::vector<Eigen::Triplet<double>> idterm_triplets;
  const auto idterm_rows = SetIdentityTerm(idterm_triplets);
  Eigen::MatrixXd affine_matrix;
  SetAffineMatrix(affine_matrix);

  opt::LinearSolver solver;
  solver.SetMatrix(idterm_triplets,idterm_rows,cols,1.0,affine_matrix);
  solver.Solve(cols, 3, result_);

  std::map<size_t ,size_t> fixed_points;
  Eigen::MatrixXd fixed_points_mat;
  fixed_points_mat.resize(fixed_vidx.size(), 3);
  fixed_points_mat.setZero();
  size_t i = 0;
  for(const auto &vidx: fixed_vidx) {
    const double distance       = e2d_.find(vidx)->second;
    const double correspondence = v2v_.find(vidx)->second;
    const Eigen::Map<const Eigen::Vector3d> normal(deformed_target_.vert_data(correspondence).normal.ptr());
    const Eigen::Map<const Eigen::Vector3d> position(deformed_target_.vert(correspondence).ptr());
    fixed_points_mat.row(i) = position + distance * normal;
    fixed_points.insert(std::pair<int,int>(vidx, correspondence));
    ++i;
  }

//  set_fixed_points_matrix(deformed_src_mesh,fixed_points,fixed_points_mat);
  std::vector<Eigen::Triplet<double>> constraint_term_triplets;
  const auto constraint_term_rows  = SetConstraintTerm(fixed_points, constraint_term_triplets);
  solver.SetMatrix(constraint_term_triplets,constraint_term_rows,cols,1.0,fixed_points_mat);

  //result's rows and cols
    solver.Solve(cols, 3, result_);
  return 0;
}
size_t pmp_jie::GarmentTransfer::SetConstraintTerm(const std::map<size_t, size_t> &fixed_points,
                                                   std::vector<Eigen::Triplet<double>> &triplets) const{
  size_t rows = 0;
  for(const auto &pair: fixed_points) {
    triplets.emplace_back(rows, pair.first, 1.0);
    ++rows;
  }
  assert(rows == fixed_points.size());
  return rows;
}
int pmp_jie::GarmentTransfer::Fit(const std::string &filename) {
  std::vector<size_t> interpenetrations;
  assert(v2v_.size() == num_verts_);
  cloth_.load(filename.c_str());
  Eigen::MatrixXd cloth_verts;
  EigenVerts(cloth_verts);

  const size_t cols = 3 * num_verts_;

  const double constants_pairs[3][2] = {{4.0, 0.8},{2.0, 0.6},{1.0, 0.4}};
  for(const auto &pair: constants_pairs) {
    opt::LinearSolver solver;
    std::vector<Eigen::Triplet<double>> interpenetration_term_triplets;
    Eigen::MatrixXd interpenetration_constants_mat;
    const auto interpenetration_term_rows = SetInterpenetrationTermAndMatrix(cloth_verts,
                                                                             interpenetrations,
                                                                             interpenetration_term_triplets,
                                                                             interpenetration_constants_mat);
    solver.SetMatrix(interpenetration_term_triplets,interpenetration_term_rows,cols,1.0,interpenetration_constants_mat);

    std::vector<Eigen::Triplet<double>> smooth_warping_term_triplets;
    const auto smooth_warping_term_rows   = SetSmoothWarpingTerm(smooth_warping_term_triplets);
    Eigen::MatrixXd smooth_warping_constants_mat;
    SetSmoothWarpingConstantsMatrix(cloth_verts, smooth_warping_constants_mat);
    solver.SetMatrix(smooth_warping_term_triplets,smooth_warping_term_rows,cols,pair[0],smooth_warping_constants_mat);

    std::vector<Eigen::Triplet<double>> damping_term_triplets;
    const auto damping_term_rows          = SetDampingTerm(damping_term_triplets);
    Eigen::MatrixXd damping_constants_mat;
    SetDampingConstantsMatrix(cloth_verts, damping_constants_mat);
    solver.SetMatrix(damping_term_triplets,damping_term_rows,cols,pair[1],damping_constants_mat);

    solver.Solve(3*num_verts_, 1, cloth_verts);
  }

  result_.setZero(3 * num_verts_, 1);
  result_ = cloth_verts;

  return 0;
}
int pmp_jie::GarmentTransfer::EigenVerts(Eigen::MatrixXd &verts) const {
  verts.setZero(3 * num_verts_, 1);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    const Eigen::Map<const Eigen::Vector3d> position(cloth_.vert(vid).ptr());
    verts(3 * vid + 0) = position(0);
    verts(3 * vid + 1) = position(1);
    verts(3 * vid + 2) = position(2);
  }
  return 0;
}
size_t pmp_jie::GarmentTransfer::SetInterpenetrationTermAndMatrix(const Eigen::MatrixXd &cloth_verts,
                                                                  std::vector<size_t> &interpenetrations,
                                                                  std::vector<Eigen::Triplet<double>> &triplets,
                                                                  Eigen::MatrixXd &mat) {
  interpenetrations.clear();
  std::vector<double> coeffs;
  size_t interpenetration_term_rows = 0;
  int cloth_vid = 0;
  auto it = e2d_.begin();
  for(const auto &pair: v2v_) {
    assert(cloth_vid == pair.first);
    assert((*it).first == cloth_vid);
    const auto body_vid = pair.second;
    const Eigen::Map<const Eigen::Vector3d> closest_point_normal(deformed_target_.vert_data(body_vid).normal.ptr());
    const Eigen::Map<const Eigen::Vector3d> closest_point(deformed_target_.vert(body_vid).ptr());
    const auto reverse_normal = -1.0 * closest_point_normal;
    const double epsilon = (*it).second;
    const double constant = reverse_normal.dot(closest_point) - epsilon;
    const Eigen::Vector3d cloth_vert = GetVert(cloth_verts, cloth_vid);
    const double penetration = reverse_normal.dot(cloth_vert - closest_point);
    if(penetration > 0) {
      SetFirstTermTriplets(interpenetration_term_rows, cloth_vid, reverse_normal, triplets);
      ++interpenetration_term_rows;
      coeffs.push_back(constant);
      interpenetrations.push_back(cloth_vid);
    }
    ++cloth_vid;
    ++it;
  }
  assert(interpenetration_term_rows == coeffs.size());
  int i = 0;
  mat.setZero(interpenetration_term_rows, 1);
  for(const auto &val: coeffs) {
    mat(i) = val;
    ++i;
  }
  return interpenetration_term_rows;
}
Eigen::Vector3d pmp_jie::GarmentTransfer::GetVert(const Eigen::MatrixXd &verts, const size_t &vid) const {
  assert(verts.cols() == 1);
  assert(verts.rows() % 3 == 0);
  return {verts(3 * vid + 0),verts(3 * vid + 1),verts(3 * vid + 2)};
}
int pmp_jie::GarmentTransfer::SetFirstTermTriplets(const size_t &row,
                                                   const size_t &vid,
                                                   const Eigen::Vector3d &coeffs,
                                                   std::vector<Eigen::Triplet<double>> &triplets) {
  //set x,y,z
  for(int i = 0; i < 3; ++i) {
    triplets.emplace_back(row, 3 * vid + i, coeffs(i));
  }
  return 0;
}
size_t pmp_jie::GarmentTransfer::SetSmoothWarpingTerm(std::vector<Eigen::Triplet<double>> &triplets) {
  for(int vid = 0; vid < num_verts_; ++vid) {
    const Eigen::Vector3d vid_coeffs(1.0, 1.0, 1.0);
    SetTriplets(vid, vid, vid_coeffs, triplets);
    const auto num_adj = cloth_.adj_v2v(vid).size();
    for(const auto &adj_vid: cloth_.adj_v2v(vid)) {
      Eigen::Vector3d adj_vid_coeffs; adj_vid_coeffs.setConstant(-1.0/(double)num_adj);
      SetTriplets(vid, adj_vid, adj_vid_coeffs, triplets);
    }
  }

  return 3 * num_verts_;
}
int pmp_jie::GarmentTransfer::SetSmoothWarpingConstantsMatrix(const Eigen::MatrixXd &cloth_verts,
                                                              Eigen::MatrixXd &mat) {
  mat.resize(3 * num_verts_, 1);
  for(int vidx = 0; vidx < num_verts_; ++vidx) {
    const auto deformation = GetLocalVertexDeformation(cloth_verts, vidx);
    for (int i = 0; i < 3; ++i) {
      mat(3 * vidx + i) = deformation(i);
    }
  }
  return 0;
}
Eigen::Vector3d pmp_jie::GarmentTransfer::GetLocalVertexDeformation(const Eigen::MatrixXd &mesh_verts,
                                                                    const int &vidx) {
  const auto vidx_position = GetVert(mesh_verts, vidx);
  Eigen::Vector3d sum_adj_vidx_position(0,0,0);
  int num_adj_verts = 0;
  for(const auto &adj_vid: cloth_.adj_v2v(vidx)) {
    const auto adj_vidx_position = GetVert(mesh_verts, adj_vid);
    sum_adj_vidx_position += adj_vidx_position;
    ++num_adj_verts;
  }
  return vidx_position - sum_adj_vidx_position / num_adj_verts;
}
int pmp_jie::GarmentTransfer::SetTriplets(const int &vidx,
                                          const int &adj_vidx,
                                          const Eigen::Vector3d &coeffs,
                                          std::vector<Eigen::Triplet<double>> &triplets) {
  //set x,y,z
  for(int i = 0; i < 3; ++i) {
    triplets.emplace_back(3 * vidx + i, 3 * adj_vidx + i, coeffs(i));
  }
  return 0;
}
size_t pmp_jie::GarmentTransfer::SetDampingTerm(std::vector<Eigen::Triplet<double>> &triplets) {
  for(int vid = 0; vid < num_verts_; ++vid) {
    const Eigen::Vector3d coeffs(1.0,1.0,1.0);
    SetTriplets(vid, vid, coeffs, triplets);
  }
  return 3*num_verts_;
}
int pmp_jie::GarmentTransfer::SetDampingConstantsMatrix(const Eigen::MatrixXd &cloth_verts, Eigen::MatrixXd &mat) {
  mat.setZero(3 * num_verts_, 1);
  for(int vid = 0; vid < num_verts_; ++vid) {
    const auto position = GetVert(cloth_verts, vid);
    mat(3 * vid + 0) = position(0);
    mat(3 * vid + 1) = position(1);
    mat(3 * vid + 2) = position(2);
  }
  return 0;
}
int pmp_jie::GarmentTransfer::WriteGarment(const std::string &filename_with_path) const {
  const auto filename = GetLoadFilename();
  std::string ofile = filename_with_path;
  if(filename_with_path.empty()) {
    ofile = pmp_jie::GetFilePath(filename)+"deformed_"+pmp_jie::GetFileName(filename, true);
  }
  std::ofstream file(ofile, std::ofstream::out);
  if(!file.is_open()) {
    std::cout << "OPEN \"" <<  ofile << "\" FAILED." << std::endl;
    return __LINE__;
  }

  for(size_t vid = 0; vid < cloth_.num_verts(); ++vid) {
    file << "v " << result_(3*vid+0) << " " << result_(3*vid+1) << " " << result_(3*vid+2) << "\n";
  }
  for(size_t pid = 0; pid < cloth_.num_polys(); ++pid) {
    file << "f ";
    for(const auto &vid: cloth_.adj_p2v(pid)) {
      file << vid + 1 << " ";
    }
    file << "\n";
  }
  file.close();
  return 0;
}
