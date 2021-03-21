//
// Created by jie zou on 2021/3/16.
//

#include "mesh_optimization_Renka16.h"
#include "local_mesh_optimization.h"

int main() {
  std::string file = DATA_PATH"/square.obj";

  pmp_jie::CPMesh mesh;
  mesh.read(file);

  pmp_jie::LocalSurfaceMeshOptimizer optimizer(&mesh);
  //pmp_jie::optimize_vert_uniform_area_Renka16
  optimizer.optimize(100, pmp_jie::optimize_vert_angle_based_Renka16);

  mesh.write("Renka_", DATA_PATH);
  return 0;
}

