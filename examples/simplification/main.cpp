//
// Created by jie zou on 2021/3/16.
//


#include "mesh_QEM_simplification.h"

int main() {
  std::string file = DATA_PATH"/sphere.obj";
  std::string res  = DATA_PATH"/QEM_sphere.obj";

  pmp_jie::Mesh mesh;
  mesh.read(file.c_str());
  mesh.simplify(100);
  mesh.write(res.c_str());

  return 0;
}

