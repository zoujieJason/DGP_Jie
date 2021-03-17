//
// Created by jie zou on 2021/3/16.
//
#include "../include/least_squares_comformal_maps.h"

Eigen::MatrixXd pmp_jie::LSCM(const cinolib::Trimesh<> &trimesh) {
  const auto num_verts = trimesh.num_verts();
  const auto num_polys = trimesh.num_polys();
  const auto num_int_verts = num_verts - 2;

  const auto bnd    = trimesh.get_ordered_boundary_vertices();
  size_t v0id = bnd[0];
  size_t v1id = bnd[1];
  if(v0id > v1id) {
    pmp_jie::Swap(v0id, v1id);
  }
  const auto v0 = trimesh.vert(v0id);
  const auto v1 = trimesh.vert(v1id);
  const double scale = (v1-v0).length();
  std::vector<size_t> vec = {v0id, v1id};

  Eigen::VectorXd uv;
  uv.resize(2*num_verts);
  uv.setZero();
  uv(0) = 0.0;
  uv(2) = 0.0;
  uv(1) = scale;
  uv(3) = 0.0;
  std::vector<Eigen::Triplet<double>> triplets;
  for(size_t pid = 0; pid < num_polys; ++pid) {
    const auto tri = trimesh.poly_verts_id(pid);
    const Eigen::Map<const Eigen::Vector3d> vi(trimesh.vert(tri[0]).ptr());
    const Eigen::Map<const Eigen::Vector3d> vj(trimesh.vert(tri[1]).ptr());
    const Eigen::Map<const Eigen::Vector3d> vk(trimesh.vert(tri[2]).ptr());
    const double x1  = 0.0;
    const double y1  = 0.0;
    const double x2  = (vj - vi).norm();
    const double y2  = 0.0;
    const double x3  = (vk - vi).dot(vj - vi) / (vj - vi).norm();
    const double val = (vk - vi).norm();
    const double y3  = std::sqrt(val * val - x3 * x3);
    const double area = std::sqrt(std::abs(x2)*std::abs(y3));
    Eigen::MatrixXd coeffs;
    coeffs.resize(3,2);
    coeffs.row(0) = Eigen::Vector2d(x3 - x2, y3 - y2);
    coeffs.row(1) = Eigen::Vector2d(x1 - x3 ,y1 - y3);
    coeffs.row(2) = Eigen::Vector2d(x2 - x1, y2 - y1);
    for(size_t i = 0; i < 3; ++i) {
      const auto adj_vid = tri[i];
      auto curr_coeff = coeffs.row(i);
      if (adj_vid == v0id) {
        curr_coeff = -1.0 * curr_coeff;
        //real part
        triplets.emplace_back(pid, 0, curr_coeff(0) / area);
        triplets.emplace_back(pid + num_polys, 2, curr_coeff(0) / area);

        //image part
        triplets.emplace_back(pid, 2, -curr_coeff(1) / area);
        triplets.emplace_back(pid + num_polys, 0, curr_coeff(1) / area);
      } else if (adj_vid == v1id) {
        curr_coeff = -1.0 * curr_coeff;
        //real part
        triplets.emplace_back(pid, 1, curr_coeff(0) / area);
        triplets.emplace_back(pid + num_polys, 3, curr_coeff(0) / area);

        //image part
        triplets.emplace_back(pid, 3, -curr_coeff(1) / area);
        triplets.emplace_back(pid + num_polys, 1, curr_coeff(1) / area);
      } else {
        const int offset = LinearSearchOffset<size_t>(vec, adj_vid);
        assert(offset != -1);
        //real part
        triplets.emplace_back(pid, adj_vid + 4 - offset, curr_coeff(0) / area);
        triplets.emplace_back(pid + num_polys, num_int_verts + adj_vid + 4 - offset, curr_coeff(0) / area);

        //image part
        triplets.emplace_back(pid, num_int_verts + adj_vid + 4 - offset, -curr_coeff(1) / area);
        triplets.emplace_back(pid + num_polys, adj_vid + 4 - offset, curr_coeff(1) / area);
      }
    }
  }
  Eigen::SparseMatrix<double> M(2*num_polys,2*num_verts);
  M.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(M.block(0,4,2*num_polys,2*(num_verts-2)));
  const Eigen::VectorXd b = M.block(0,0,2*num_polys,4) * uv.block(0,0,4,1);
  uv.block(4,0,2*(num_verts-2),1) = solver.solve(b);

  Eigen::MatrixXd UV;
  UV.resize(num_verts, 2);
  UV.row(v0id) = Eigen::Vector2d(0.0,0.0);
  UV.row(v1id) = Eigen::Vector2d(scale,0.0);
  for(size_t vid = 0; vid < trimesh.num_verts(); ++vid) {
    const int offset = LinearSearchOffset<size_t>(vec, vid);
    if(offset==-1) continue;
    UV.row(vid) = Eigen::Vector2d(uv(vid+4-offset), uv(vid+4+num_verts-2-offset));
  }
  return UV;
}