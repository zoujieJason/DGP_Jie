//
// Created by jie zou on 2021/3/16.
//

#include "tools.h"
#include "abf_pp.h"

int main() {
  using namespace pmp_jie;

  std::string file   = DATA_PATH"/head.obj";
  std::string uv_res = DATA_PATH"/abf_head.obj";

  Eigen::VectorXd alpha;
  pmp_jie::ABFPara abf_para;
  abf_para.Load(file);
  abf_para.SetVerbose(true);
  abf_para.ToAngleSpace(alpha);

  pmp_jie::AngleRecons angle_recons;
  angle_recons.LoadMesh(abf_para.MeshRef());
  Eigen::MatrixXd verts;
  angle_recons.To2DSpace(alpha, verts);

  WriteUVToObj(verts, abf_para.MeshRef().vector_polys(), uv_res);

  return 1;
}
