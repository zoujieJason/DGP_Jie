//
// Created by jie zou on 2021/3/21.
//
#include "../include/mesh_optimization_Renka16.h"
double pmp_jie::optimize_vert_angle_based_Renka16(pmp_jie::Base *mesh_ptr, const size_t &vid) {
  auto *polymesh_ptr_ = dynamic_cast<pmp_jie::CPMesh *>(mesh_ptr);

  auto &trimesh = polymesh_ptr_->mesh_ref();

  const cinolib::vec3d n(0.0,0.0,1.0);
  const auto v0 = trimesh.vert(vid);
  cinolib::vec3d e1; {
    for(const auto &adj_vid: trimesh.adj_v2v(vid)) {
      e1 = trimesh.vert(adj_vid) - v0;
      if(!e1.is_degenerate()) {
        e1.normalize();
        break; }
    }
  }
  if(e1.is_degenerate()) {
    return 0;
  }
  cinolib::vec3d e2 = n.cross(e1);

  const auto e0 = area_energy(trimesh, vid);

  double a = 0.0, b = 0.0, c = 0.0;
  double f1 = 0.0, f2 = 0.0;
  for(const auto &adj_pid: trimesh.adj_v2p(vid)) {
    int off = trimesh.poly_vert_offset(adj_pid, vid);
    const auto vnext = trimesh.poly_vert_id(adj_pid, (off+1)%3);
    off = trimesh.poly_vert_offset(adj_pid, vnext);
    const auto vprev = trimesh.poly_vert_id(adj_pid, (off+1)%3);
    const auto ek = trimesh.vert(vprev) - trimesh.vert(vnext);
    const double wk = 1.0/ek.length();
    if(ek.length()==0) {
      std::cerr << "ek.length()==0." << std::endl;
      return 0;
    }
    a += wk * ek.dot(e2) * ek.dot(e2);
    b += wk * ek.dot(e1) * ek.dot(e2);
    c += wk * ek.dot(e1) * ek.dot(e1);
    const auto w  = trimesh.vert(vnext)-v0;
    const auto rk = ek.length_squared()*(w) - w.dot(ek)*ek;
    f1 += wk * rk.dot(e1);
    f2 += wk * rk.dot(e2);
  }
  if(b==0) {
    std::cerr << "b==0." << std::endl;
    return 0;
  }
  if(c*a/b-b==0) {
    std::cerr << "c*a/b-b==0" << std::endl;
    return 0;
  }
  const double s1 = (f2+c*f1/b)/(c*a/b-b);
  const double s2 = (a*s1-f1)/b;
  const cinolib::vec3d vnew =v0+s1*e1+s2*e2;
  if((vnew-v0).length() < 1e-8) {
    return 0;
  }
  trimesh.vert(vid) = vnew;
  return area_energy(trimesh,vid) - e0;
}
double pmp_jie::area_energy(const cinolib::Polygonmesh<> &trimesh, const size_t &vid) {
  double sum = 0.0;
  for(const auto &adj_pid: trimesh.adj_v2v(vid)) {
    const auto area = trimesh.poly_area(adj_pid);
    sum += area*area;
  }
  return 0.5*sum;
}
double pmp_jie::optimize_vert_uniform_area_Renka16(pmp_jie::Base *mesh_ptr, const size_t &vid) {
  auto *polymesh_ptr_ = dynamic_cast<pmp_jie::CPMesh *>(mesh_ptr);

  auto &trimesh = polymesh_ptr_->mesh_ref();

  const cinolib::vec3d n(0.0,0.0,1.0);
  const auto v0 = trimesh.vert(vid);
  cinolib::vec3d e1; {
    for(const auto &adj_vid: trimesh.adj_v2v(vid)) {
      e1 = trimesh.vert(adj_vid) - v0;
      if(!e1.is_degenerate()) {
        e1.normalize();
        break; }
    }
  }
  if(e1.is_degenerate()) {
    return 0;
  }
  cinolib::vec3d e2 = n.cross(e1);

  const auto e0 = area_energy(trimesh, vid);

  double a = 0.0, b = 0.0, c = 0.0;
  double f1 = 0.0, f2 = 0.0;
  for(const auto &adj_pid: trimesh.adj_v2p(vid)) {
    int off = trimesh.poly_vert_offset(adj_pid, vid);
    const auto vnext = trimesh.poly_vert_id(adj_pid, (off+1)%3);
    off = trimesh.poly_vert_offset(adj_pid, vnext);
    const auto vprev = trimesh.poly_vert_id(adj_pid, (off+1)%3);
    const auto ek = trimesh.vert(vprev) - trimesh.vert(vnext);
    const double wk = 1.0;
    if(ek.length()==0) {
      std::cerr << "ek.length()==0." << std::endl;
      return 0;
    }
    a += wk * ek.dot(e2) * ek.dot(e2);
    b += wk * ek.dot(e1) * ek.dot(e2);
    c += wk * ek.dot(e1) * ek.dot(e1);
    const auto w  = trimesh.vert(vnext)-v0;
    const auto rk = ek.length_squared()*(w) - w.dot(ek)*ek;
    f1 += wk * rk.dot(e1);
    f2 += wk * rk.dot(e2);
  }
  if(b==0) {
    std::cerr << "b==0." << std::endl;
    return 0;
  }
  if(c*a/b-b==0) {
    std::cerr << "c*a/b-b==0" << std::endl;
    return 0;
  }
  const double s1 = (f2+c*f1/b)/(c*a/b-b);
  const double s2 = (a*s1-f1)/b;
  const cinolib::vec3d vnew =v0+s1*e1+s2*e2;
  if((vnew-v0).length() < 1e-8) {
    return 0;
  }
  trimesh.vert(vid) = vnew;
  return area_energy(trimesh,vid) - e0;
}
