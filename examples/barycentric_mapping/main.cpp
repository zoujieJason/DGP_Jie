//
// Created by jie zou on 2021/3/16.
//
#include "barycentric_mapping.h"

int main() {
  using namespace pmp_jie;

  std::string file   = DATA_PATH"/face.obj";
  std::string uv_res = DATA_PATH"/face_uv.obj";

  cinolib::Trimesh<> trimesh(file.c_str());
  const auto uv = BaryMapping(trimesh, UNIFORM, CircleBoundary(trimesh.get_boundary_vertices().size()));
  WriteUVToObj(uv, trimesh.vector_polys(), uv_res);

  return 1;
}
