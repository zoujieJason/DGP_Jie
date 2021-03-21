//
// Created by jie zou on 2020/12/18.
//
#include "../include/angle_based_optimization_Xu.h"
pmp_jie::Xu06::Xu06(pmp_jie::Base *base_mesh, const size_t &vid) {
  mesh_ptr_ = dynamic_cast<pmp_jie::CPMesh *>(base_mesh);
  vid_ = vid;
  dim_of_x_ = 2;
  dim_of_f_ = mesh_ptr_->mesh_.adj_v2v(vid).size();
}
size_t pmp_jie::Xu06::dim_of_x() const {
  return dim_of_x_;
}
size_t pmp_jie::Xu06::dim_of_f() const {
  return dim_of_f_;
}
int pmp_jie::Xu06::jac(const double *x_ptr, Eigen::SparseMatrix<double> &JT) {
  const auto v0 = mesh_ptr_->vert_data(vid_);
  Eigen::Map<const Eigen::VectorXd> vert(x_ptr, dim_of_x_);
  mesh_ptr_->update_vert(vid_, Eigen::Vector3d(vert.x(),vert.y(),0.0));

  const auto adj_vind = mesh_ptr_->mesh_.adj_v2v(vid_);
  std::vector<Eigen::Triplet<double>> triplets;
  for(size_t i = 0; i < adj_vind.size(); ++i) {
    auto abc = get_bisector_line(vid_, adj_vind[i]);
    const auto dot = abc.dot(Eigen::Vector3d(vert(0), vert(1), 1.0));
    const auto normalize = std::sqrt(abc(0)*abc(0)+abc(1)*abc(1));
    if(dot<0) abc = -abc;
    triplets.emplace_back(0,i,abc(0)/normalize);
    triplets.emplace_back(1,i,abc(1)/normalize);
  }
  JT.resize(2,dim_of_f_);
  JT.setFromTriplets(triplets.begin(), triplets.end());

  mesh_ptr_->update_vert(vid_,v0);
  return 0;
}
double pmp_jie::Xu06::energy(const double *x_ptr) {
  const auto v0 = mesh_ptr_->vert_data(vid_);
  Eigen::Map<const Eigen::VectorXd> vert(x_ptr, dim_of_x_);
  mesh_ptr_->update_vert(vid_, Eigen::Vector3d(vert.x(),vert.y(),0.0));

  double s(0.0);
  for(const auto &adj_vid: mesh_ptr_->mesh_.adj_v2v(vid_)) {
    const auto abc = get_bisector_line(vid_,adj_vid);
    const auto dot = abc.dot(Eigen::Vector3d(vert(0), vert(1), 1.0));
    const auto normalize = abc(0)*abc(0)+abc(1)*abc(1);
    s += (dot*dot/normalize);
  }

  mesh_ptr_->update_vert(vid_,v0);
  return s;
}
int pmp_jie::Xu06::val(const double *x_ptr, double *f_ptr) {
  const auto v0 = mesh_ptr_->vert_data(vid_);
  Eigen::Map<const Eigen::VectorXd> vert(x_ptr, dim_of_x_);
  mesh_ptr_->update_vert(vid_, Eigen::Vector3d(vert.x(),vert.y(),0.0));

  Eigen::Map<Eigen::VectorXd> f(f_ptr,dim_of_f_); f.setZero();
  const auto adj_vind = mesh_ptr_->mesh_.adj_v2v(vid_);
  for(size_t i = 0; i < adj_vind.size(); ++i) {
    const auto abc       = get_bisector_line(vid_,adj_vind[i]);
    const auto dot       = abc.dot(Eigen::Vector3d(vert(0), vert(1), 1.0));
    const auto normalize = std::sqrt(abc(0)*abc(0)+abc(1)*abc(1));
    f(i) = std::abs(dot)/normalize;
  }

  mesh_ptr_->update_vert(vid_,v0);
  return 0;
}
void pmp_jie::Xu06::set_vid(const size_t &vid) {
  vid_ = vid;
  dim_of_f_ = mesh_ptr_->mesh_.adj_v2v(vid).size();
}
Eigen::Vector3d pmp_jie::Xu06::get_bisector_line(const size_t &vid, const size_t &adj_vid) const {
  const auto adj_vert = mesh_ptr_->vert_data(adj_vid);
  const auto bisector = angle_bisector_line(vid, adj_vid);
  return {bisector(1),-bisector(0),-adj_vert(0)*bisector(1)+adj_vert(1)*bisector(0)};
}
Eigen::Vector3d pmp_jie::Xu06::angle_bisector_line(const size_t &vid, const size_t &adj_vid) const {
  const auto &mesh = mesh_ptr_->mesh_;
  const auto &polys = mesh.adj_e2p(mesh.edge_id(vid, adj_vid));
  size_t left_poly(0), right_poly(0);
  if(mesh.poly_verts_are_CCW(polys[0],vid, adj_vid)) {
    left_poly = polys[1];
    right_poly = polys[0];
  }
  else {
    left_poly = polys[0];
    right_poly = polys[1];
  }
  const auto left_curr = mesh.poly_vert_offset(left_poly,adj_vid);
  const auto right_curr = mesh.poly_vert_offset(right_poly,adj_vid);
  const Eigen::Vector3d adj_vert = mesh_ptr_->vert_data(adj_vid);
  const Eigen::Vector3d prev_vert = mesh_ptr_->vert_data(mesh.poly_vert_id(left_poly, (left_curr + 1) % 3));
  const Eigen::Vector3d next_vert = mesh_ptr_->vert_data(mesh.poly_vert_id(right_poly, (right_curr - 1 + 3) % 3));
  const Eigen::Vector3d u = prev_vert-adj_vert;
  const Eigen::Vector3d v = next_vert-adj_vert;
  const auto angle = mesh.poly_angle_at_vert(left_poly,adj_vid)+mesh.poly_angle_at_vert(right_poly,adj_vid);
  const auto _u = pmp_jie::PlanRotate({u(0), u(1)},angle/2.0,true).normalized();
  const auto _v = pmp_jie::PlanRotate({v(0), v(1)},angle/2.0).normalized();
  const auto dot = std::max(-1.0,std::min(1.0,_u.dot(_v)));
  if(std::abs(dot-1.0) > 1e-6) {
    std::cerr << dot << std::endl;
    std::cerr << "FUNCTION AngleBisectorUnitPoint : std::abs(dot-1.0) > epsilon_" << std::endl;
  }
  return {_u(0),_u(1),0.0};
}
double pmp_jie::LM_optimize_vert_Xu06(pmp_jie::Base *mesh_ptr, const size_t &vid) {

  const auto v0 = mesh_ptr->vert_data(vid);
  Eigen::Vector2d v(v0.x(),v0.y());

  pmp_jie::Xu06 func(mesh_ptr,vid);

  opt::LM_Gauss_Newton solver(func);
  int num_iters = 1;
  solver.solve(&v[0],func,num_iters);

  mesh_ptr->update_vert(vid,Eigen::Vector3d(v.x(),v.y(),0.0));
  return func.energy(&v0[0])-func.energy(&v[0]);
}
