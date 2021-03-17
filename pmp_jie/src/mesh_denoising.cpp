//
// Created by jie zou on 2020/8/18.
//
#include "../include/mesh_denoising.h"
pmp_jie::BilFilter::BilFilter()=default;
int pmp_jie::BilFilter::Smoothing() {
  NormalFilter();
  Recovery();
  return 0;
}
int pmp_jie::BilFilter::UsrAllocateVars() {
  InitNormals();
  InitVerts();
  volume_ = trimesh_.mesh_volume();
  std::cout << "Mesh Volume: " << volume_ << std::endl;
  return 0;
}
int pmp_jie::BilFilter::InitNormals(){
  normals_.resize(3, num_polys_);
  normals_.setZero();
  size_t pid(0);
  for(const auto &vec: trimesh_.vector_poly_normals()) {
    const Eigen::Map<const Eigen::Vector3d> normal(vec.ptr());
    normals_.col(pid) = normal;
    ++pid;
  }
  return 0;
}
double pmp_jie::BilFilter::SolveCurrNormalIteration() {
  Eigen::MatrixXd normals;
  normals.resize(3, num_polys_);
  normals.setZero();
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    double kp(0.0);
    const auto ci = GetPolyCenter(pid);
    const auto ni = normals_.col(pid);
    for(const auto &adj_pid: trimesh_.adj_p2p(pid)) {
      const auto cj = GetPolyCenter(adj_pid);
      const auto Aj = trimesh_.poly_area(adj_pid);
      const auto nj = normals_.col(adj_pid);
      const double wj = Aj * pmp_jie::GaussValue((ci-cj).norm(),sigma_) * pmp_jie::GaussValue((ni-nj).norm(),sigma_);
      kp += wj;
      normals.col(pid) = normals.col(pid) + wj * nj;
    }
    normals.col(pid) = normals.col(pid) / kp;
  }
  const auto error = (normals - normals_).norm();
  normals_ = normals;
  return error;
}
int pmp_jie::BilFilter::InitVerts() {
  verts_.resize(3,num_verts_);
  verts_.setZero();
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    const Eigen::Map<const Eigen::Vector3d> position(trimesh_.vert(vid).ptr());
    verts_.col(vid) = position;
  }
  return 0;
}
int pmp_jie::BilFilter::NormalFilter() {
  size_t it(0);
  double error(0.0);
  do {
    error = SolveCurrNormalIteration();
    ++it;
  }while(error > epsilon_ && it < 100);
  std::cout << "Smoothing Iteration: " << it  << ". Normal Error: " << error << std::endl ;
  return 0;
}
int pmp_jie::BilFilter::Recovery() {
  Eigen::MatrixXi tris;
  GetTriangles(tris);
  size_t it(0);
  double error(0.0);
  do {
    error = SolveCurrVertIteration(tris);
    ++it;
  }while(error > epsilon_ && it < 100);
  std::cout << "Smoothing Iteration: " << it  << ". Vert Error: " << error << std::endl ;
  return 0;
}
double pmp_jie::BilFilter::SolveCurrVertIteration(const Eigen::MatrixXi &tris) {
  Eigen::MatrixXd verts;
  verts.resize(3,num_verts_);
  verts.setZero();
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    size_t cnt(0);
    const auto xi = verts_.col(vid);
    for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
      const auto nj = normals_.col(adj_pid);
      const auto cj = GetPolyCenter(adj_pid);
      verts.col(vid) = verts.col(vid) + nj.dot(cj-xi) * nj;
      ++cnt;
    }
    verts.col(vid) = verts.col(vid) / (double)cnt + verts_.col(vid);
  }
  const auto error = (verts - verts_).norm();
  //rejust volume
//  const auto curr_volume = dgp::GetVolume(verts, tris);
//  std::cout << "curr_volume: " << curr_volume << std::endl;
//  verts_ = (curr_volume / volume_) * verts;
  verts_ = verts;
  return error;
}
Eigen::Vector3d pmp_jie::BilFilter::GetPolyCenter(const size_t &pid) const {
  assert(pid>=0 && pid<num_polys_);
  Eigen::Vector3d c;
  c.setZero();
  for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
    const auto position = verts_.col(adj_vid);
    c = c + position;
  }
  return c / 3.0;
}
int pmp_jie::BilFilter::Write() const {
  const auto filename = GetLoadFilename();
  const std::string ofile = pmp_jie::GetFilePath(filename)+"de_"+pmp_jie::GetFileName(filename,true);
  std::ofstream file(ofile, std::ofstream::out);
  if(!file.is_open()) {
    std::cout << "OPEN \"" <<  ofile << "\" FAILED." << std::endl;
    return __LINE__;
  }

  for(size_t vid = 0; vid < trimesh_.num_verts(); ++vid) {
    file << "v " << verts_(0,vid) << " " << verts_(1,vid) << " " << verts_(2,vid) << "\n";
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
