//
// Created by jie zou on 2019/11/7.
//

#ifndef PROGRESSIVEMESHES_MESH_H
#define PROGRESSIVEMESHES_MESH_H

#include "QEM_structures.h"

namespace pmp_jie {
class Mesh {
 public:
  int read(const char *filePath);
  int simplify(const int &resultNumber);
  int refine(const int &resultNumber);
  int write(const char *filePath) const;

  int nVertices() const;

 private:
  std::map<int, std::shared_ptr<Vertex>> verticesMap;
  std::set<std::shared_ptr<Triangle>> triangleSet;
  std::vector<int> permutation;
  std::vector<int> map;
  std::vector<VSplit> vsplits;

  int addTriangle(std::weak_ptr<Vertex> v1,
                  std::weak_ptr<Vertex> v2,
                  std::weak_ptr<Vertex> v3);
  int computeAllCollapse();
  int computeCollapseAtVertex(std::weak_ptr<Vertex> vertex);
  std::weak_ptr<Vertex> getMinimumCostVertex() const;
  int edgeCollapse();
  int vertexSplit();
  int deleteVertex(std::weak_ptr<Vertex> vertex);
  int deleteTriangle(std::weak_ptr<Triangle> tri);
  int updateCollapseAroundVertex(std::weak_ptr<Vertex> vertex);

  int setVsplit(std::weak_ptr<Vertex> vertex);
  void printPermutation() const;
};
}


#endif //PROGRESSIVEMESHES_MESH_H
