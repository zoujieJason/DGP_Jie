//
// Created by jie zou on 2020/8/7.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_HALFEDGE_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_HALFEDGE_H_

#include "structures.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include <Eigen/Dense>

namespace pmp_jie {

class HMesh {
 public:
  HMesh();

  int Read(const std::string &file);
  int AddFace(const int &v1, const int &v2, const int &v3);

  size_t NumVerts() const;
  size_t NumFaces() const;
  bool IsBoundary(const size_t &vid) const;
  std::vector<std::weak_ptr<Halfedge>> Halfedges(const int &vid) const;
  std::vector<std::weak_ptr<Halfedge>> Verts(const int &fid) const;

 private:
  std::vector<Eigen::Vector3f> position_container_;
  std::vector<Eigen::Vector3f> normal_container_;
  std::vector<std::shared_ptr<Vertex>> vertices_container_;
  std::vector<std::shared_ptr<Halfedge>> halfedges_container_;
  std::vector<std::shared_ptr<Face>> faces_container_;
  std::map<int, std::shared_ptr<Edge>> edges_container_map_;

  std::weak_ptr<Halfedge> GetHalfedge(const std::weak_ptr<Edge>& edge, const int &from_vid, const int &to_vid) const;
  std::weak_ptr<Edge> FindEdge(const int &from_vid, const int &to_vid) const;
  std::weak_ptr<Edge> GetEdge(const int &from_vid, const int &to_vid);
  std::weak_ptr<Edge> AddEdge(const int &v1, const int &v2);
  std::weak_ptr<Face> AddTriangle(const int &fid, const int &vid);
  int ConnectTriangle(const std::weak_ptr<Halfedge> &halfedge12,
                      const std::weak_ptr<Halfedge> &halfedge23,
                      const std::weak_ptr<Halfedge> &halfedge31,
                      const int &fid);
};

}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_HALFEDGE_H_
