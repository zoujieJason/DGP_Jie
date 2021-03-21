//
// Created by jie zou on 2021/3/16.
//

#include "angle_based_optimization_Xu.h"
#include "local_mesh_optimization.h"

int main() {
  std::string file = DATA_PATH"/square.obj";

  pmp_jie::CPMesh mesh;
  mesh.read(file);

  pmp_jie::LocalSurfaceMeshOptimizer optimizer(&mesh);
  optimizer.optimize(100, pmp_jie::LM_optimize_vert_Xu06);

  mesh.write("xu_", DATA_PATH);
  return 0;
}

