//
// Created by jie zou on 2021/3/16.
//
#include "barycentric_mapping.h"

int main() {
  using namespace pmp_jie;

  std::string file   = DATA_PATH"/head.obj";
  std::string uv_res = DATA_PATH"/head_uv.obj";

  cinolib::Trimesh<> trimesh(file.c_str());
  const auto uv = barycentric_mapping(trimesh, UNIFORM, circle_boundary(trimesh.get_boundary_vertices().size()));
  write_uv_to_obj(uv, trimesh.vector_polys(), uv_res);

  return 1;
}
