//
// Created by jie zou on 2020/12/18.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_POLYGON_MESH_H_
#define DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_POLYGON_MESH_H_

#include <cinolib/meshes/meshes.h>

#include "base_mesh.h"
#include "tools.h"

namespace pmp_jie {

class CPMesh: public pmp_jie::Base {
 public:
  //io
  void read(const std::string &filename) override;
  void write(const std::string &prefix, const std::string &path) const;

  //mesh data
  size_t num_verts() const override;
  Eigen::Vector3d vert_data(const size_t &vid) const override;
  void update_vert(const size_t &vid, const Eigen::Vector3d &vert) override;
  std::vector<uint> poly_tri_verts(const size_t &pid, const size_t &vid) const;
  size_t angle_id(const size_t &pid, const size_t &vid) const;
  double vert_angle(const size_t &vid) const;
  double poly_angle_at_vert(const size_t &pid, const size_t &vid) const;
  double vert_1ring_min_length(const size_t &vid) const;

  //topo
  const std::vector<uint> & adj_v2v(const size_t &vid) const override;
  bool vert_is_boundary(const size_t &vid) const override;
  bool poly_is_corner(const size_t &pid) const;
  std::vector<uint> boundary_adj_v2v(const size_t &vid) const;
  int get_poly_ccw_vert(const size_t &pid, const size_t &vid) const;
  int get_poly_cw_vert(const size_t &pid, const size_t &vid) const;
  bool check_vert_is_local_interior_point(const size_t &vid, const Eigen::Vector3d &vert);
  bool check_vert_is_local_interior_point(const size_t &vid) const override;

  //opt data
  void update_optimal_vars();
  double optimal_var(const size_t &pid, const size_t &vid) const;

  //derivatives
  Eigen::Vector3d arccos_grad_poly_at_vert(const size_t &pid, const size_t &vid, const size_t &dvid) const;

 public:

  cinolib::Polygonmesh<> &mesh_ref() { return mesh_; }
  cinolib::Polygonmesh<> mesh_;

 private:
  Eigen::VectorXd opt_vars_;
  size_t nv_{0};
};

}

#endif //DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_POLYGON_MESH_H_
