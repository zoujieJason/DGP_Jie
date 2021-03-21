//
// Created by jie zou on 2020/12/18.
//
#include "../include/polygon_mesh.h"
size_t pmp_jie::CPMesh::num_verts() const {
  return mesh_.num_verts();
}
bool pmp_jie::CPMesh::vert_is_boundary(const size_t &vid) const {
  return mesh_.vert_is_boundary(vid);
}
Eigen::Vector3d pmp_jie::CPMesh::vert_data(const size_t &vid) const {
  return Eigen::Map<const Eigen::Vector3d>(mesh_.vert(vid).ptr());
}
void pmp_jie::CPMesh::update_vert(const size_t &vid, const Eigen::Vector3d &vert) {
  mesh_.vert(vid).x() = vert.x();
  mesh_.vert(vid).y() = vert.y();
  mesh_.vert(vid).z() = vert.z();
}
void pmp_jie::CPMesh::read(const std::string &filename) {
  mesh_.load(filename.c_str());
  nv_ = mesh_.verts_per_poly(0);
  //info
  {
    if(nv_==3) {
      std::cout << "Trimesh";
    }
    if(nv_==4) {
      std::cout << "Quadmesh";
    }
    std::cout << "( " << mesh_.mesh_data().filename << " )." << std::endl;
  }
  update_optimal_vars();
}
void pmp_jie::CPMesh::write(const std::string &prefix, const std::string &path) const {
  const std::string out = path + "/" + prefix + cinolib::get_file_name(mesh_.mesh_data().filename);
  mesh_.save(out.c_str());
}
void pmp_jie::CPMesh::update_optimal_vars() {
  //vertex average.
  opt_vars_.setConstant(nv_ * mesh_.num_polys(), 0.0);
  for(size_t vid = 0; vid < mesh_.num_verts(); ++vid) {
    int num_adj = mesh_.adj_v2p(vid).size();
    double sum = vert_angle(vid);
    for(const auto &adj_pid: mesh_.adj_v2p(vid)) {
      if(poly_is_corner(adj_pid)) {
        const auto alpha = poly_angle_at_vert(adj_pid, vid);
        if(alpha>=0) {
          --num_adj;
          sum -= alpha;
          opt_vars_(angle_id(adj_pid, vid)) = alpha;
        }
      }
    }
    //average
    const double opt_angle = sum/(double)num_adj;
    for(const auto &adj_pid: mesh_.adj_v2p(vid)) {
      const auto aid = angle_id(adj_pid, vid);
      if(opt_vars_(aid)==0.0) {
        opt_vars_(aid) = opt_angle;
      }
    }
  }
}
double pmp_jie::CPMesh::optimal_var(const size_t &pid, const size_t &vid) const {
  const auto aid = angle_id(pid, vid);
  return opt_vars_(aid);
}
size_t pmp_jie::CPMesh::angle_id(const size_t &pid, const size_t &vid) const {
  return nv_*pid + mesh_.poly_vert_offset(pid,vid);
}
double pmp_jie::CPMesh::vert_angle(const size_t &vid) const {
  if(!mesh_.vert_is_boundary(vid)) return 2.0*M_PI;
  const auto adj_verts = boundary_adj_v2v(vid);
  assert(adj_verts.size()==2);
  cinolib::vec3d u = mesh_.vert(adj_verts[0])-mesh_.vert(vid);
  cinolib::vec3d v = mesh_.vert(adj_verts[1])-mesh_.vert(vid);
  double angle = u.angle_rad(v);
  const cinolib::vec3d n(0.0,0.0,1.0);
  if(u.cross(v).dot(n)<0) {
    angle = 2.0*M_PI-angle;
  }
  return angle;
}
std::vector<uint> pmp_jie::CPMesh::boundary_adj_v2v(const size_t &vid) const {
  if(!mesh_.vert_is_boundary(vid)) return {};
  const auto list = mesh_.vert_ordered_verts_link(vid);
  const auto lvid = list[0];
  const auto rvid = list[list.size() - 1];

  //verify
  bool left_flag = false;
  bool right_flag = false;
  for(const auto &adj_pid: mesh_.adj_v2p(vid)) {
    if(mesh_.poly_contains_vert(adj_pid, lvid)) {
      if(mesh_.poly_verts_are_CCW(adj_pid, lvid, vid)) {
        left_flag = true;
      }
    }
    if(mesh_.poly_contains_vert(adj_pid, rvid)) {
      if(mesh_.poly_verts_are_CCW(adj_pid, vid, rvid)) {
        right_flag = true;
      }
    }
  }
  if(!(left_flag && right_flag)) {
    std::cerr << "!(left_flag && right_flag)" << std::endl;
  }
  assert(left_flag && right_flag);
  return {lvid, rvid};
}
bool pmp_jie::CPMesh::poly_is_corner(const size_t &pid) const {
  for(const auto &adj_vid: mesh_.adj_p2v(pid)) {
    if(!mesh_.vert_is_boundary(adj_vid)) {
      return false;
    }
  }
  return true;
}
double pmp_jie::CPMesh::poly_angle_at_vert(const size_t &pid, const size_t &vid) const {
  const cinolib::vec3d n(0.0,0.0,1.0);
  const auto tri = poly_tri_verts(pid, vid);
  cinolib::vec3d u = mesh_.vert(tri[1])-mesh_.vert(tri[0]);
  cinolib::vec3d v = mesh_.vert(tri[2])-mesh_.vert(tri[0]);
  double angle = u.angle_rad(v);
  if(u.cross(v).dot(n)<0) {
    angle = 2.0*M_PI-angle;
  }
  return angle;
}
std::vector<uint> pmp_jie::CPMesh::poly_tri_verts(const size_t &pid, const size_t &vid) const {
  std::vector<uint> tri;
  const int offset = mesh_.poly_vert_offset(pid, vid);
  const auto vj    = mesh_.poly_vert_id(pid, (offset+1)%nv_);
  const auto vk    = mesh_.poly_vert_id(pid, (offset-1+nv_)%nv_);
  return {static_cast<unsigned int>(vid), vj, vk};
}
int pmp_jie::CPMesh::get_poly_ccw_vert(const size_t &pid, const size_t &vid) const {
  const int offset = mesh_.poly_vert_offset(pid, vid);
  const int vnext  = mesh_.poly_vert_id(pid, (offset+1)%nv_);
  if(mesh_.poly_verts_are_CCW(pid, vnext, vid)) {
    return vnext;
  }
  else {
    std::cerr << "FUNCTION get_poly_ccw_vert: error. " << __LINE__ << std::endl;
    return -1;
  }
}
int pmp_jie::CPMesh::get_poly_cw_vert(const size_t &pid, const size_t &vid) const {
  const int offset = mesh_.poly_vert_offset(pid, vid);
  const int vprev  = mesh_.poly_vert_id(pid, (offset-1+nv_)%nv_);
  if(mesh_.poly_verts_are_CCW(pid, vid, vprev)) {
    return vprev;
  }
  else {
    std::cerr << "FUNCTION get_poly_cw_vert: error. " << __LINE__ << std::endl;
    return -1;
  }
}
Eigen::Vector3d pmp_jie::CPMesh::arccos_grad_poly_at_vert(const size_t &pid,
                                                               const size_t &vid,
                                                               const size_t &dvid) const {
  const auto tri = poly_tri_verts(pid, vid);
  const Eigen::Vector3d vi = vert_data(tri[0]);
  const Eigen::Vector3d vj = vert_data(tri[1]);
  const Eigen::Vector3d vk = vert_data(tri[2]);
  const double a = (vk-vj).norm();
  const double b = (vk-vi).norm();
  const double c = (vj-vi).norm();
  double h = std::max(-1.0, std::min(1.0,(b*b+c*c-a*a)/(2.0*b*c)));
  //collinear: theta = 0, h = -1.0; theta = pi, h = 1.0.
  double arc  = -1.0/(std::sqrt(1-h*h));
  ///debug.
  if(std::abs(h)==1.0) {
    arc=-1.0/(1e-6);
    return {0.0,0.0,0.0};
  }
  //d arccos/d vk;
  const Eigen::Vector3d dk = arc*((vj-vi)/(b*c)-h*(vk-vi)/(b*b));
  const Eigen::Vector3d dj = arc*((vk-vi)/(b*c)-h*(vj-vi)/(c*c));
  if(dvid == tri[2]) {
    assert(!pmp_jie::Vec3dIsNan(dk));
    assert(!pmp_jie::Vec3dIsInf(dk));
    return dk;
  }
  //d arccos/d vj
  if(dvid == tri[1]) {
    assert(!pmp_jie::Vec3dIsNan(dj));
    assert(!pmp_jie::Vec3dIsInf(dj));
    return dj;
  }
  assert(dvid == tri[0]);
  assert(!pmp_jie::Vec3dIsNan(-dk-dj));
  assert(!pmp_jie::Vec3dIsInf(-dk-dj));
  return -dk-dj;
}
double pmp_jie::CPMesh::vert_1ring_min_length(const size_t &vid) const {
  double min_len(std::numeric_limits<double>::max());
  for(const auto &eid: mesh_.adj_v2e(vid)) {
    min_len = std::min(min_len, mesh_.edge_length(eid));
  }
  return min_len;
}
bool pmp_jie::CPMesh::check_vert_is_local_interior_point(const size_t &vid, const Eigen::Vector3d &vert) {
  const Eigen::Vector3d n(0.0, 0.0, 1.0);
  for (const auto &adj_pid: mesh_.adj_v2p(vid)) {
    const auto tri = poly_tri_verts(adj_pid, vid);
    Eigen::Vector3d u = vert_data(tri[1])-vert; u.normalize();
    Eigen::Vector3d v = vert_data(tri[2])-vert; u.normalize();
    if (u.cross(v).dot(n) <= 0) return false;
  }
  return true;
}
const std::vector<uint> &pmp_jie::CPMesh::adj_v2v(const size_t &vid) const {
  return mesh_.adj_v2v(vid);
}
bool pmp_jie::CPMesh::check_vert_is_local_interior_point(const size_t &vid) const {
  const Eigen::Vector3d n(0.0, 0.0, 1.0);
  for (const auto &adj_pid: mesh_.adj_v2p(vid)) {
    const auto tri = poly_tri_verts(adj_pid, vid);
    Eigen::Vector3d u = vert_data(tri[1])-vert_data(tri[0]); u.normalize();
    Eigen::Vector3d v = vert_data(tri[2])-vert_data(tri[0]); u.normalize();
    if (u.cross(v).dot(n) <= 0) return false;
  }
  return true;
}
