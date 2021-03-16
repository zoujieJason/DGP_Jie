//
// Created by jie zou on 2021/3/16.
//
#include "least_squares_comformal_maps.h"
#include "tools.h"

int main() {
  using namespace pmp_jie;

  std::string file   = DATA_PATH"/face.obj";
  std::string uv_res = DATA_PATH"/lscm_face.obj";

  cinolib::Trimesh<> trimesh(file.c_str());
  const auto uv = lscm(trimesh);
  write_uv_to_obj(uv, trimesh.vector_polys(), uv_res);

  return 1;
}
