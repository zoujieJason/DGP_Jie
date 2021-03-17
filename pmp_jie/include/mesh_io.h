//
// Created by jie zou on 2020/7/23.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_MESH_IO_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_MESH_IO_H_

#include <string>
#include <cinolib/meshes/meshes.h>

namespace pmp_jie {

class CTriMesh {
 public:
  CTriMesh();

  int Load(const std::string &file);
  std::string GetLoadFilename() const;

  cinolib::Trimesh<> & MeshRef();

 protected:
  int AllocateVars();
  virtual int UsrAllocateVars()=0;

  size_t GetAngleId(const size_t &pid, const size_t &vid) const;
  int GetPolyCCWVert(const size_t &pid, const size_t &vid) const;
  double SumVertAngle(const size_t &vid) const;
  int GetBVertNeighbor(const size_t &vid, size_t &lvid, size_t &rvid) const;
  Eigen::Vector2d GetBVertNeighborLens(const size_t &vid) const;
  int GetTriangles(Eigen::MatrixXi &tris) const;

  //trimesh
  cinolib::Trimesh<> trimesh_;

  //mesh data
  size_t num_bod_verts_{0};
  size_t num_int_verts_{0};
  size_t num_angles_{0};
  size_t num_polys_{0};
  size_t num_verts_{0};

 private:
  int AllocateBaseVars();

};
}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_MESH_IO_H_
