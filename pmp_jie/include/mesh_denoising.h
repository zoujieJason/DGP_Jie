//
// Created by jie zou on 2020/8/18.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_MESH_SMOOTHING_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_MESH_SMOOTHING_H_

#include "mesh_io.h"
#include "tools.h"

namespace pmp_jie {

class BilFilter: public pmp_jie::CTriMesh {
 public:
  BilFilter();

  int Smoothing();
  int Write() const;

 private:
  int UsrAllocateVars() override;
  int InitNormals();
  int InitVerts();
  int NormalFilter();
  double SolveCurrNormalIteration();
  int Recovery();
  double SolveCurrVertIteration(const Eigen::MatrixXi &tris);

  Eigen::Vector3d GetPolyCenter(const size_t &pid) const;

  Eigen::MatrixXd normals_;
  Eigen::MatrixXd verts_;
  double sigma_{5.0};
  double epsilon_{1e-5};
  double volume_{0.0};

};
}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_MESH_SMOOTHING_H_
