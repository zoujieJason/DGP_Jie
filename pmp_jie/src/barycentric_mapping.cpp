//
// Created by jie zou on 2021/3/16.
//
#include "../include/barycentric_mapping.h"
Eigen::MatrixXd pmp_jie::barycentric_mapping(const cinolib::Trimesh<> &trimesh,
                                             const pmp_jie::WEIGHT_TYPE &w_type,
                                             const Eigen::MatrixXd &bnd_uv) {
  const auto num_verts = trimesh.num_verts();
  Eigen::MatrixXd UV;
  UV.setZero(num_verts, 2);

  //vmap
  std::map<int, int> vmap;
  size_t int_vid = 0;
  for(size_t vid = 0; vid < num_verts; ++vid) {
    if(!trimesh.vert_is_boundary(vid)) {
      vmap.insert(std::pair<size_t, size_t>(vid, int_vid));
      ++int_vid;
    }
  }

  const auto bnd_vind = trimesh.get_ordered_boundary_vertices();
  std::vector<uint> ordered_bnd_vind = bnd_vind;
  std::sort(ordered_bnd_vind.begin(), ordered_bnd_vind.end());

  //Fill matrix A.
  std::vector<Eigen::Triplet<double>> int_triplets;
  std::vector<Eigen::Triplet<double>> bnd_triplets;
  for(size_t vid = 0; vid < num_verts; ++vid) {
    if(trimesh.vert_is_boundary(vid)) continue;
    auto n_i = trimesh.adj_v2v(vid).size();
    Eigen::VectorXd w; w.setZero(n_i);
    for(size_t j = 0; j < n_i; ++j) {
      const auto adj_vid = trimesh.adj_v2v(vid)[j];
      const auto eid = trimesh.edge_id(vid, adj_vid);
      switch (w_type) {
        case UNIFORM: w(j) = 1.0/(double)n_i;
          break;
        case HARMONIC: {
          const auto len = trimesh.edge_length(eid);
          const auto &polys = trimesh.adj_e2p(eid);
          const auto laradius = trimesh.poly_angle_at_vert(polys[0], adj_vid) / 2.0;
          const auto raradius = trimesh.poly_angle_at_vert(polys[1], adj_vid) / 2.0;
          w(j) = (std::tan(laradius) + std::tan(raradius)) / len;
        } break;
        case COTLAPLACE: w(j) = trimesh.edge_cotangent_weight(eid); break;
      }
    }
    w /= w.sum();
    for(size_t j = 0; j < n_i; ++j) {
      const auto adj_vid = trimesh.adj_v2v(vid)[j];
      if (trimesh.vert_is_boundary(adj_vid)) {
        bnd_triplets.emplace_back(vmap.at(vid), linear_search<uint>(ordered_bnd_vind, adj_vid), w(j));
      }
      else {
        int_triplets.emplace_back(vmap.at(vid), vmap.at(adj_vid), w(j));
      }
    }
    int_triplets.emplace_back(vmap.at(vid), vmap.at(vid), -1.0);
  }

  const auto num_bnd_verts = bnd_vind.size();
  Eigen::MatrixXd circle;
  circle.setZero(num_bnd_verts, 2);
  for(size_t i = 0; i < num_bnd_verts; ++i) {
    circle.row(linear_search<uint>(ordered_bnd_vind, bnd_vind[i])) = bnd_uv.row(i);
    const size_t vid = bnd_vind[i];
    UV.row(vid) = bnd_uv.row(i);
//    trimesh.vert(bnd_vind[i]) = cinolib::vec3d(bnd_uv(i,0), bnd_uv(i,1), 0.0);
  }

  //Solve.
  assert(num_verts>=num_bnd_verts);
  const size_t num_itr_verts = num_verts - num_bnd_verts;
  Eigen::SparseMatrix<double> A(num_itr_verts, num_itr_verts);
  A.setFromTriplets(int_triplets.begin(), int_triplets.end());
  Eigen::SparseMatrix<double> B(num_itr_verts, num_bnd_verts);
  B.setFromTriplets(bnd_triplets.begin(), bnd_triplets.end());

  Eigen::MatrixXd result;
  opt::least_squares_solve_ATA(A,-B*circle,result);

  //Write to uv.
  int_vid = 0;
  for(size_t vid = 0; vid < num_verts; ++vid) {
    if(!trimesh.vert_is_boundary(vid)) {
      UV.row(vid) = result.row(int_vid);
//      trimesh.vert(vid) = cinolib::vec3d(result(int_vid,0), result(int_vid,1), 0.0);
      ++int_vid;
    }
  }
  return UV;
}
Eigen::MatrixXd pmp_jie::circle_boundary(const size_t &n, const double &scale) {
  Eigen::MatrixXd circle;
  circle.setZero(n, 2);
  const double angle(2.0 * M_PI/(double)n);
  for(size_t vid = 0; vid < n; ++vid) {
    const Eigen::Vector2d vuv(scale*std::cos(angle * (double)vid), scale*std::sin(angle * (double)vid));
    circle.row(vid) = vuv;
  }
  return circle;
}
