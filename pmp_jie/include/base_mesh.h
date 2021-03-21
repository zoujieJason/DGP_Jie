//
// Created by jie zou on 2020/12/18.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_BASE_MESH_H_
#define DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_BASE_MESH_H_

#include <string>
#include <Eigen/Dense>

namespace pmp_jie {

class Base {
 public:
  virtual void read(const std::string &filename) = 0;
  virtual bool vert_is_boundary(const size_t &vid) const = 0 ;
  virtual size_t num_verts() const = 0;
  virtual Eigen::Vector3d vert_data(const size_t &vid) const = 0;
  virtual void update_vert(const size_t &vid, const Eigen::Vector3d &vert) = 0;
  //For parallel
  virtual const std::vector<uint> & adj_v2v(const size_t &vid) const = 0;
  //For untangle
  virtual bool check_vert_is_local_interior_point(const size_t &vid) const = 0;

};

}

#endif //DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_BASE_MESH_H_
