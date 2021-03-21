//
// Created by jie zou on 2020/12/18.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_LOCAL_MESH_OPTIMIZATION_H_
#define DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_LOCAL_MESH_OPTIMIZATION_H_

#include <iostream>
#include <omp.h>

#include "base_mesh.h"
#include "function.h"
#include "timer.h"
#include "tools.h"

namespace pmp_jie {

class LocalSurfaceMeshOptimizer {
 public:
  LocalSurfaceMeshOptimizer(pmp_jie::Base *mesh);
  int optimize(const size_t &num_iters, double (*optimization_function)(pmp_jie::Base *, const size_t &));
  int parallel_optimize(const size_t &num_iters, double (*optimization_function)(pmp_jie::Base *, const size_t &));

  int local_untangle(const size_t &num_inner_iters, double (*optimization_function)(pmp_jie::Base *, const size_t &));

 private:
  double untangle_vert(const size_t &num_iters,
                       const size_t &vid,
                       double (*optimization_function)(pmp_jie::Base *,
                                                       const size_t &));
  bool greedy_graph_coloring(std::vector<std::vector<size_t>> &colored_vertices);

  pmp_jie::Base *mesh_ptr_;

};

}

#endif //DIGITAL_GEOMETRY_PROCESSING_DGP_INCLUDE_LOCAL_MESH_OPTIMIZATION_H_
