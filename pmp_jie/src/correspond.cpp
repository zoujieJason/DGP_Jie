//
// Created by jie zou on 2020/8/7.
//
#include "correspond.h"
pmp_jie::Correspond::Correspond()=default;
int pmp_jie::Correspond::LoadTarget(const std::string &target) {
  target_.load(target.c_str());
  BuildTree();
  return 0;
}
int pmp_jie::Correspond::UsrAllocateVars() {
  return 0;
}
int pmp_jie::Correspond::BuildTree() {
  for (int pidx = 0; pidx < target_.num_polys(); ++pidx) {
    Eigen::Matrix3d V;
    V.setZero();
    int i = 0;
    for(const auto &vidx: target_.adj_p2v(pidx)) {
      const Eigen::Map<const Eigen::Vector3d> position(target_.vert(vidx).ptr());
      V.col(i) = position;
      ++i;
    }
    Point_3 v0 = Eigen2CGAL(V.col(0));
    Point_3 v1 = Eigen2CGAL(V.col(1));
    Point_3 v2 = Eigen2CGAL(V.col(2));
    tris_.emplace_back(v0, v1, v2);
  }
  tree_.rebuild(tris_.begin(),tris_.end());
  return 0;
}
Point_3 pmp_jie::Correspond::Eigen2CGAL(const Eigen::Vector3d &vec) const {
  return {vec.x(), vec.y(), vec.z()};
}
Eigen::Vector3d pmp_jie::Correspond::CGAL2Eigen(const Point_3 &vec) const {
  return {vec.x(), vec.y(), vec.z()};
}
Point_3 pmp_jie::Correspond::GetClosestPoint(const Point_3 &query) const {
  Point_3 closest_point = tree_.closest_point(query);
  return closest_point;
}
size_t pmp_jie::Correspond::GetClosestPolyId(const Point_3 &query) const{
  Point_and_primitive_id pp = tree_.closest_point_and_primitive(query);
  CIterator id = pp.second;
  std::size_t index = std::distance(tris_.begin(),id);
  assert( tris_[index] == *id );
  return index;
}
int pmp_jie::Correspond::LoadFixedPoints(const std::string &mapfile) {
  fixed_points_.clear();
  std::ifstream file(mapfile);
  if(!file.is_open()) {
    std::cout << "OPEN \"" <<  mapfile << "\" FAILED." << std::endl;
    return __LINE__;
  }
  std::string buffer;
  while (std::getline(file, buffer)) {
    std::istringstream str(buffer);
    std::string symbol;
    str >> symbol;
    if(symbol == "p") {
      std::string x, y;
      str >> x >> y;
      fixed_points_.insert(std::pair<int,int>(std::stoi(x), std::stoi(y)));
    }
  }
  file.close();
  return 0;
}
int pmp_jie::Correspond::Deform(const std::string &mapfile) {
  LoadFixedPoints(mapfile);
  if(fixed_points_.empty()) {
    std::cerr << "[STAG]: CORRESPONDENCE - No Fixed Points Map. " << std::endl;
    return 0;
  }
  //init
  const auto cols = num_verts_ + num_polys_;
  opt::LinearSolver init_solver;
  std::vector<Eigen::Triplet<double>> smoothness_term_triplets;
  const auto smoothness_term_rows = SetSmoothnessTerm(smoothness_term_triplets);
  init_solver.SetMatrix(smoothness_term_triplets,smoothness_term_rows,cols,1.0);

  std::vector<Eigen::Triplet<double>> identity_term_triplets;
  const auto identiry_term_rows = SetIdentityTerm(identity_term_triplets);
  Eigen::MatrixXd identity_matrix;
  SetIdentityMatrix(identiry_term_rows,identity_matrix);
  init_solver.SetMatrix(identity_term_triplets,identiry_term_rows,cols,0.001,identity_matrix);

  std::vector<Eigen::Triplet<double>> constraint_term_triplets;
  const auto constraint_term_rows = SetConstraintTerm(constraint_term_triplets);
  Eigen::MatrixXd fixed_verts;
  SetConstraintMatrix(fixed_verts);
  init_solver.SetMatrix(constraint_term_triplets,constraint_term_rows,cols,1000000.0,fixed_verts);

  Eigen::MatrixXd init_result;
  init_solver.Solve(cols,3,init_result);
  UpdateClosestVaildPoints(init_result);

  //solve problem iteratively
  const double constants[4] = {1.0, 100.0, 1000.0, 5000.0};
  for(const auto &constant: constants) {
    opt::LinearSolver solver;
    solver.SetMatrix(smoothness_term_triplets,smoothness_term_rows,cols,1.0);
    solver.SetMatrix(identity_term_triplets,identiry_term_rows,cols,0.001,identity_matrix);

    std::vector<Eigen::Triplet<double>> valid_points_term_triplets;
    const auto valid_points_term_rows = SetVaildPointsTerm(valid_points_term_triplets);
    solver.SetMatrix(valid_points_term_triplets,valid_points_term_rows,cols,constant, vaild_points_);

    Eigen::MatrixXd result;
    solver.Solve(cols, 3, result);
    UpdateClosestVaildPoints(result);
  }

  FindCorrespondence();
  return 0;
}
size_t pmp_jie::Correspond::SetSmoothnessTerm(std::vector<Eigen::Triplet<double>> &triplets) const {
  int curr_row = 0;
  for(int pid = 0; pid < num_polys_; ++pid) {
    for(const auto &adj_pid: trimesh_.adj_p2p(pid)) {
      SetPolyTriplets(pid,    curr_row, 1.0,triplets);
      SetPolyTriplets(adj_pid,curr_row,-1.0,triplets);
      curr_row += 3;
    }
  }
  return curr_row;
}
size_t pmp_jie::Correspond::SetIdentityTerm(std::vector<Eigen::Triplet<double>> &triplets) const {
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    SetPolyTriplets(pid,3*pid,1.0,triplets);
  }
  return 3 * num_polys_;
}
size_t pmp_jie::Correspond::SetConstraintTerm(std::vector<Eigen::Triplet<double>> &triplets) const {
  int rows = 0;
  for(const auto &pair: fixed_points_) {
    triplets.emplace_back(rows, pair.first, 1.0);
    ++rows;
  }
  assert(rows == fixed_points_.size());
  return rows;
}
size_t pmp_jie::Correspond::SetPolyTriplets(const size_t &pid,
                                            const size_t &curr_row,
                                            const double &weight,
                                            std::vector<Eigen::Triplet<double>> &triplets) const {
  const auto inv_gradient_mat = pmp_jie::GetPolyGradientMatrix(trimesh_, pid);
  const auto tri = trimesh_.adj_p2v(pid);
  //set parameter 3 times for each face
  for(int time = 0; time < 3; ++time) {
    const Eigen::Vector3d coeffs = inv_gradient_mat.col(time);
    triplets.emplace_back(curr_row + time, tri[0],         -1.0*weight*coeffs.sum());
    triplets.emplace_back(curr_row + time, tri[1],         weight*coeffs(0));
    triplets.emplace_back(curr_row + time, tri[2],         weight*coeffs(1));
    triplets.emplace_back(curr_row + time, num_verts_+pid, weight*coeffs(2));
  }
  return 0;
}
int pmp_jie::Correspond::SetIdentityMatrix(const size_t &rows, Eigen::MatrixXd &matrix) const {
  assert(rows % 3 == 0);
  matrix.resize(rows,3);
  matrix.setOnes();
  Eigen::Matrix3d I;
  I.setIdentity();
  for(int i = 0; i < rows / 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      matrix.row(3*i+j) = I.row(j);
    }
  }
  return 0;
}
int pmp_jie::Correspond::SetConstraintMatrix(Eigen::MatrixXd &matrix) const {
  matrix.resize(fixed_points_.size(), 3);
  matrix.setZero();
  size_t vid(0);
  for(const auto &pair: fixed_points_) {
    const auto target_vidx = pair.second;
    const Eigen::Map<const Eigen::Vector3d> position(target_.vert(target_vidx).ptr());
    matrix.row(vid) = position;
    ++vid;
  }
  return 0;
}
int pmp_jie::Correspond::UpdateClosestVaildPoints(const Eigen::MatrixXd &curr_result) {
  vaild_points_.resize(num_verts_, 3);
  vaild_points_.setZero();

  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!fixed_points_.count(vid)) {
      vaild_points_.row(vid) = CGAL2Eigen(GetClosestPoint(Eigen2CGAL(curr_result.row(vid))));
    }
    else {
      const auto target_vidx = fixed_points_.find(vid)->second;
      const Eigen::Map<const Eigen::Vector3d> position(target_.vert(target_vidx).ptr());
      vaild_points_.row(vid) = position;
    }
  }
  return 0;
}
size_t pmp_jie::Correspond::SetVaildPointsTerm(std::vector<Eigen::Triplet<double>> &triplets) const {
  for(int i = 0; i < num_verts_; ++i) {
    triplets.emplace_back(i,i,1.0);
  }
  return num_verts_;
}
int pmp_jie::Correspond::WriteCorrespondMesh(const std::string &filename_with_path) {
  const auto filename = GetLoadFilename();
  std::string ofile = filename_with_path;
  if(filename_with_path.empty()) {
    ofile = pmp_jie::GetFilePath(filename)+"corres_"+pmp_jie::GetFileName(filename, true);
  }
  std::ofstream file(ofile, std::ofstream::out);
  if(!file.is_open()) {
    std::cout << "OPEN \"" <<  ofile << "\" FAILED." << std::endl;
    return __LINE__;
  }

  for(size_t vid = 0; vid < num_verts_; ++vid) {
    file << "v " << vaild_points_(vid,0) << " " << vaild_points_(vid,1) << " " << vaild_points_(vid,2) << "\n";
  }

  for(size_t pid = 0; pid < num_polys_; ++pid) {
    file << "f ";
    for(const auto &vid: trimesh_.adj_p2v(pid)) {
      file << vid + 1 << " ";
    }
    file << "\n";
  }
  file.close();
  return 0;
}
int pmp_jie::Correspond::FindCorrespondence() {
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    Eigen::Vector3d centroid(0.0,0.0,0.0);
    for(const auto &adj_vidx: trimesh_.adj_p2v(pid)) {
      centroid += vaild_points_.row(adj_vidx);
    }
    centroid /= 3.0;
    p2p_.insert(std::pair<size_t, size_t>(pid, GetClosestPolyId(Eigen2CGAL(centroid))));
  }
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    const auto cor_pid = GetClosestPolyId(Eigen2CGAL(vaild_points_.row(vid)));
    const auto cor_vid = GetClosestVertId(cor_pid, vaild_points_.row(vid));
    assert(cor_vid != -1);
    v2v_.insert(std::pair<size_t, size_t>(vid, cor_vid));
    const Eigen::Map<const Eigen::Vector3d> pos(trimesh_.vert(vid).ptr());
    e2d_.insert(std::pair<size_t, double>(vid, (pos-CGAL2Eigen(GetClosestPoint(Eigen2CGAL(pos)))).norm()));
  }
  return 0;
}
int pmp_jie::Correspond::GetClosestVertId(const size_t &cor_pid, const Eigen::Vector3d &pos) const {
  double min_distance = 99999.9;
  int min_distance_adj_vidx = -1;
  for(const auto &adj_vidx: target_.adj_p2v(cor_pid)) {
    const Eigen::Map<const Eigen::Vector3d> adj_vert_position(target_.vert(adj_vidx).ptr());
    if(min_distance > (adj_vert_position - pos).norm()) {
      min_distance = (adj_vert_position - pos).norm();
      min_distance_adj_vidx = adj_vidx;
    }
  }
  return min_distance_adj_vidx;
}
