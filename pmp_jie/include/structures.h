//
// Created by jie zou on 2020/8/7.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_STRUCTURES_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_STRUCTURES_H_

#include <memory>

namespace pmp_jie {

class Vertex;
class Halfedge;
class Face;
class Edge;

class Vertex {
 public:
  Vertex(): idx(-1){}

  int idx;
  std::weak_ptr<Halfedge> halfedge;
};

class Face {
 public:
  Face(): idx(-1) {}

  int idx;
  std::weak_ptr<Halfedge> halfedge;
};

class Halfedge {
 public:
  Halfedge(): vertex_idx(-1), face_idx(-1) {}

  int vertex_idx;
  int face_idx;
  std::weak_ptr<Halfedge> opposite_halfedge;
  std::weak_ptr<Halfedge> next_halfedge;
  std::weak_ptr<Halfedge> prev_halfedge;
  std::weak_ptr<Edge> edge;

};

class Edge {
 public:
  Edge()= default;
  std::weak_ptr<Halfedge> from_halfedge;
  std::weak_ptr<Halfedge> to_halfedge;

  std::shared_ptr<Edge> next_edge;
};

}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_STRUCTURES_H_
