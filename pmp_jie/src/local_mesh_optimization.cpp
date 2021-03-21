//
// Created by jie zou on 2020/12/18.
//
#include "../include/local_mesh_optimization.h"
pmp_jie::LocalSurfaceMeshOptimizer::LocalSurfaceMeshOptimizer(pmp_jie::Base *mesh) {
  mesh_ptr_ = mesh;
}
int pmp_jie::LocalSurfaceMeshOptimizer::optimize(const size_t &num_iters,
                                                       double (*optimization_function)(pmp_jie::Base *,
                                                                                       const size_t &)) {
  Timer<> timer;
  timer.beginStage("Optimize.");
  size_t cnt(0);
  while(cnt < num_iters) {
    double residual(-1.0);
    for(size_t vid = 0; vid < mesh_ptr_->num_verts(); ++vid) {
      if(!mesh_ptr_->vert_is_boundary(vid)) {
        residual = std::max(residual, optimization_function(mesh_ptr_, vid));
      }
    }
    ++cnt;
    std::cout << "Iter( " << cnt << " ), Residual( " << residual << " ). " << std::endl;
    if(residual<1e-6&&residual>=0) {
      std::cout << "CONVERGENCE. Last residual( " << residual << " )" << std::endl;
      break;
    }
    if(residual==-1) {
      std::cout << "" << std::endl;
      break;
    }
  }
  timer.endStage();
  return cnt;
}
int pmp_jie::LocalSurfaceMeshOptimizer::parallel_optimize(const size_t &num_iters,
                                                                double (*optimization_function)(pmp_jie::Base *,
                                                                                                const size_t &)) {
  std::cout << "Parallel optimize: " << std::endl;
  std::vector<std::vector<size_t>> colored_vertices;
  if(!greedy_graph_coloring(colored_vertices)) {
    std::cout << "coloring failed." << std::endl;
  }
  const auto num_procs = omp_get_num_procs();
  std::cout << "Num procs( " << num_procs << " ). " << std::endl;
  Timer<> timer;
  timer.beginStage("Parallel optimize.");
  size_t cnt = 0;
  while(cnt < num_iters) {
    double residual(-1.0);
    for(const auto &vertices: colored_vertices) {
      omp_set_num_threads(num_procs);
//#pragma omp parallel for schedule(dynamic)
#pragma omp parallel for
      {
        for(size_t vid = 0; vid < vertices.size(); ++vid) {
#pragma omp critical
          {
            if(!mesh_ptr_->vert_is_boundary(vid)) {
              residual = std::max(residual, optimization_function(mesh_ptr_, vid));
            }
          }
        }
      }
    }
    ++cnt;
    std::cout << "Iter( " << cnt << " ), Residual( " << residual << " ). " << std::endl;
    if(residual<1e-6&&residual>=0) {
      std::cout << "CONVERGENCE. Last residual( " << residual << " )" << std::endl;
      break;
    }
  }
  timer.endStage();
  return cnt;
}
bool pmp_jie::LocalSurfaceMeshOptimizer::greedy_graph_coloring(std::vector<std::vector<size_t>> &colored_vertices) {
  std::vector<int> vertices_color(mesh_ptr_->num_verts(),-1);
  vertices_color[0] = 0;
  std::vector<bool> color_at_vert_not_available(mesh_ptr_->num_verts(), false);
  int max_num_colors = 0;
  for(size_t vid = 1; vid < mesh_ptr_->num_verts(); ++vid) {
    for(const auto &adj_vid: mesh_ptr_->adj_v2v(vid)) {
      if(vertices_color[adj_vid]!=-1) {
        color_at_vert_not_available[vertices_color[adj_vid]] = true;
      }
    }
    int cr = 0;
    for(;cr < mesh_ptr_->num_verts(); ++cr) {
      if(!color_at_vert_not_available[cr]) {
        break;
      }
    }
    vertices_color[vid] = cr;
    max_num_colors = std::max(cr,max_num_colors);

    for(const auto &adj_vid: mesh_ptr_->adj_v2v(vid)) {
      if(vertices_color[adj_vid]!=-1) {
        color_at_vert_not_available[vertices_color[adj_vid]] = false;
      }
    }
  }
  colored_vertices.resize(max_num_colors+1);
  for(size_t vid = 0; vid < mesh_ptr_->num_verts(); ++vid) {
    colored_vertices[vertices_color[vid]].push_back(vid);
  }
  std::set<int> vind;
  for(const auto &vertices: colored_vertices) {
    for(const auto &vid: vertices) {
      if(vind.count(vid)) {
        return false;
      }
      else {
        vind.insert(vid);
      }
    }
  }
  if(vind.size()!=mesh_ptr_->num_verts()) return false;
  return true;
}
int pmp_jie::LocalSurfaceMeshOptimizer::local_untangle(const size_t &num_inner_iters,
                                                             double (*optimization_function)(pmp_jie::Base *,
                                                                                             const size_t &)) {
  Timer<> timer;
  timer.beginStage("Untangle.");
  for(size_t vid = 0; vid < mesh_ptr_->num_verts(); ++vid) {
    if(!mesh_ptr_->vert_is_boundary(vid)&&!mesh_ptr_->check_vert_is_local_interior_point(vid)) {
      //optimize vertex domain.
      for(const auto &adj_vid: mesh_ptr_->adj_v2v(vid)) {
        untangle_vert(num_inner_iters, adj_vid, optimization_function);
      }
      //optimize vertex.
      untangle_vert(num_inner_iters, vid, optimization_function);
      if(mesh_ptr_->check_vert_is_local_interior_point(vid)) {
        std::cout << "Vert( " << vid << " ) untangle succeed." << std::endl;
      }
      else {
        std::cout << "Vert( " << vid << " ) untangle failed." << std::endl;
      }
    }
  }
  timer.endStage();
  return 0;
}
double pmp_jie::LocalSurfaceMeshOptimizer::untangle_vert(const size_t &num_iters,
                                                               const size_t &vid,
                                                               double (*optimization_function)(pmp_jie::Base *,
                                                                                               const size_t &)) {
  size_t cnt = 0;
  double sum_residual = 0.0;
  while(cnt < num_iters) {
    const auto residual = optimization_function(mesh_ptr_, vid);
    if(residual<1e-6) {
      return sum_residual;
    }
    ++cnt;
  }
  return sum_residual;
}
