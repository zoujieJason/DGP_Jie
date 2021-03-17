//
// Created by jie zou on 2021/3/16.
//

#include "mesh_denoising.h"

int main() {
  std::string file   = DATA_PATH"/Octahedron_n2.obj";

  pmp_jie::BilFilter bil_filter;
  bil_filter.Load(file);
  bil_filter.Smoothing();
  bil_filter.Write();
  return 0;
}

