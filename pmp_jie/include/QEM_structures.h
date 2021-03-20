//
// Created by jie zou on 2019/11/7.
//

#ifndef PROGRESSIVEMESHES_STRUCTURES_H
#define PROGRESSIVEMESHES_STRUCTURES_H


#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <map>
#include <vector>
#include <wchar.h>
#include <assert.h>

#include <Eigen/Dense>

namespace pmp_jie {
class Vertex;
class Triangle;

class Vertex {
 public:
  Eigen::Vector3f position;
  int id;
  std::vector<std::weak_ptr<Vertex>> neighbor;
  std::vector<std::weak_ptr<Triangle>> face;
  float cost;
  std::weak_ptr<Vertex> collapseVertex;

  Vertex(const Eigen::Vector3f &pos, const int &idRef);

  int removeTriangle(std::weak_ptr<Triangle> tri);
  int removeNeighbor(std::weak_ptr<Vertex> vertex);
  int addTriangle(std::weak_ptr<Triangle> tri);
  int addNeighbor(std::weak_ptr<Vertex> vertex);
  int deleteNeighbor(std::weak_ptr<Vertex> vertex);
  int deleteTriangle(std::weak_ptr<Triangle> tri);

  bool hasTriangle(std::weak_ptr<Triangle> tri) const;
  bool hasNeighbor(std::weak_ptr<Vertex> vertex) const;
  Eigen::Matrix4f getMatQ() const;
  std::vector<std::weak_ptr<Triangle>> getTrianglesBetweenVlAndVr(std::weak_ptr<Vertex> vl,
                                                                  std::weak_ptr<Vertex> vr) const;

 private:
  std::weak_ptr<Triangle> getCwTriangleFromVertex(std::weak_ptr<Vertex> vertex) const;
  std::weak_ptr<Triangle> getCCwTriangleFromVertex(std::weak_ptr<Vertex> vertex) const;
};

class Triangle {
 public:
  std::weak_ptr<Vertex> v0;
  std::weak_ptr<Vertex> v1;
  std::weak_ptr<Vertex> v2;

  Eigen::Vector3f normal;

  Triangle(std::weak_ptr<Vertex> v0Ref, std::weak_ptr<Vertex> v1Ref, std::weak_ptr<Vertex> v2Ref);

  void computeNormal();
  int replaceVertex(std::weak_ptr<Vertex> vold, std::weak_ptr<Vertex> vnew);
  bool hasVertex(std::weak_ptr<Vertex> vertex) const;
  std::weak_ptr<Vertex> getNextVertex(std::weak_ptr<Vertex> vertex) const;
  std::weak_ptr<Vertex> getPervVertex(std::weak_ptr<Vertex> vertex) const;
};

class VSplit {
 public:
  Eigen::Vector3f vfPos;
  Eigen::Vector3f vtPos;
  Eigen::Vector4i indices;

  VSplit();

  void setVfPos(const Eigen::Vector3f &pos);
  void setVtPos(const Eigen::Vector3f &pos);
  void setVfId(const int &id);
  void setVtId(const int &id);
  void setVlId(const int &id);
  void setVrId(const int &id);
};
}
#endif //PROGRESSIVEMESHES_STRUCTURES_H
