//
// Created by jie zou on 2020/7/23.
//
#include "../include/mesh_io.h"
pmp_jie::CTriMesh::CTriMesh()=default;
int pmp_jie::CTriMesh::Load(const std::string &file) {
  trimesh_.load(file.c_str());
  AllocateVars();
  return 0;
}
std::string pmp_jie::CTriMesh::GetLoadFilename() const {
  return trimesh_.mesh_data().filename;
}
int pmp_jie::CTriMesh::AllocateVars() {
  AllocateBaseVars();
  UsrAllocateVars();
  return 0;
}
int pmp_jie::CTriMesh::AllocateBaseVars() {
  //constants
  num_bod_verts_ = trimesh_.get_boundary_vertices().size();
  num_int_verts_ = trimesh_.num_verts() - num_bod_verts_;
  num_polys_     = trimesh_.num_polys();
  num_verts_     = trimesh_.num_verts();
  num_angles_    = 3 * num_polys_;
  return 0;
}
size_t pmp_jie::CTriMesh::GetAngleId(const size_t &pid, const size_t &vid) const {
  return 3 * pid + trimesh_.poly_vert_offset(pid, vid);
}
double pmp_jie::CTriMesh::SumVertAngle(const size_t &vid) const {
  double sum_angle(0.0);
  for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
    sum_angle += std::max(trimesh_.poly_angle_at_vert(adj_pid, vid), 2.0 * M_PI / 360.0) ;
  }
  return sum_angle;
}
int pmp_jie::CTriMesh::GetPolyCCWVert(const size_t &pid, const size_t &vid) const {
  const auto curr_voffset = trimesh_.poly_vert_offset(pid, vid);
  const auto next_voffset = (curr_voffset + 1) % 3;
  const auto next_vid     = trimesh_.poly_vert_id(pid, next_voffset);
  if(trimesh_.poly_verts_are_CCW(pid,next_vid,vid)) {
    return next_vid;
  }
  else {
    return -1;
  }
}
int pmp_jie::CTriMesh::GetBVertNeighbor(const size_t &vid, size_t &lvid, size_t &rvid) const {
  if(!trimesh_.vert_is_boundary(vid)) {
    return -1;
  }
  const auto list = trimesh_.vert_ordered_verts_link(vid);
  lvid   = list[0];
  rvid  = list[list.size() - 1];

  //verify
  bool left_flag = false;
  bool right_flag = false;
  for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
    if(trimesh_.poly_contains_vert(adj_pid, lvid)) {
      if(trimesh_.poly_verts_are_CCW(adj_pid, lvid, vid)) {
        left_flag = true;
      }
    }
    if(trimesh_.poly_contains_vert(adj_pid, rvid)) {
      if(trimesh_.poly_verts_are_CCW(adj_pid, vid, rvid)) {
        right_flag = true;
      }
    }
  }
  if(!(left_flag && right_flag)) {
    std::cerr << "!(left_flag && right_flag)" << std::endl;
  }
  assert(left_flag && right_flag);
  return 1;
}
Eigen::Vector2d pmp_jie::CTriMesh::GetBVertNeighborLens(const size_t &vid) const {
  Eigen::Vector2d lens;
  lens.setOnes();
  if(trimesh_.vert_is_boundary(vid)) {
    size_t left_vid(0), right_vid(0);
    GetBVertNeighbor(vid, left_vid, right_vid);
    const auto left_eid  = trimesh_.edge_id(vid, left_vid);
    const auto right_eid = trimesh_.edge_id(vid, right_vid);
    assert(left_eid != -1 && right_eid != -1);
    lens(0) *= trimesh_.edge_length(left_eid);
    lens(1) *= trimesh_.edge_length(right_eid);
  }
  return lens;
}
int pmp_jie::CTriMesh::GetTriangles(Eigen::MatrixXi &tris) const {
  tris.resize(3,num_polys_);
  tris.setZero();
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    Eigen::Vector3i tri;
    tri.setZero();
    const auto vec = trimesh_.adj_p2v(pid);
    for(size_t vid = 0; vid < 3; ++vid) {
      tri(vid) = (int)vec[vid];
    }
    tris.col(pid) = tri;
  }
  return 0;
}
cinolib::Trimesh<> &pmp_jie::CTriMesh::MeshRef()  {
  return this->trimesh_;
}


