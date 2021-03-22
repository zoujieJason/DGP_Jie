//
// Created by jie zou on 2021/3/16.
//
#include "../include/tools.h"
int pmp_jie::WriteUVToObj(const Eigen::MatrixXd &uv,
                          const std::vector<std::vector<uint>> &polys,
                          const std::string &filename) {
  std::ofstream os(filename, std::ofstream::out);
  if(!os.is_open()) {
    std::cerr << "OPEN " << filename << " FAILED." << std::endl;
    return __LINE__;
  }

  for(size_t vid = 0; vid < uv.rows(); ++vid) {
    os << "v " << uv.row(vid).x() << " " << uv.row(vid).y() << " 0\n";
  }
  for(const auto &poly: polys) {
    os << "f " << poly[0]+1 << " " << poly[1]+1 << " " <<  poly[2]+1 << "\n";
  }
  os.close();
  return 1;
}
std::string pmp_jie::GetFilePath(const std::string &string) {
  size_t pos = string.find_last_of('/');
  if(pos>=string.size()) return "./";
  return string.substr(0,pos+1);
}
std::string pmp_jie::GetFileName(const std::string &string, const bool &with_extension) {
  size_t pos = string.find_last_of('/');
  std::string substring = (pos>=string.size()) ? string : string.substr(pos+1);
  if(with_extension) return substring;
  pos = substring.find_last_of('.');
  return substring.substr(0,pos);
}
std::string pmp_jie::GetFileExtension(const std::string &string) {
  size_t pos = string.find_last_of('.');
  if(pos>=string.size()) return "";
  return string.substr(pos+1);
}
double pmp_jie::GaussValue(const double &x, const double &sigma) {
  const double val = sqrt(2 * M_PI) * sigma;
  return exp(-1.0 * x * x / sigma / sigma) / val;
}
Eigen::Matrix3d pmp_jie::GetPolyGradientMatrix(const cinolib::Trimesh<> &trimesh,
                                               const size_t &pid,
                                               const bool &inv_sign) {
  Eigen::Matrix3d Q = pmp_jie::GetPolyMatrix(trimesh,pid);

  const Eigen::Vector3d v1 = Q.col(0);
  Q.col(0) = Q.col(1) - v1;
  Q.col(1) = Q.col(2) - v1;
  Q.col(2) = Q.col(0).cross(Q.col(1));
  Q.col(2) /= Q.col(2).norm();

  if(inv_sign) {
    return Q.inverse();
  }
  else {
    return Q;
  }
}
Eigen::Matrix3d pmp_jie::GetPolyMatrix(const cinolib::Trimesh<> &trimesh, const size_t &pid) {
  Eigen::Matrix3d V;
  V.setZero();
  int num_col = 0;
  for(const auto &vidx: trimesh.adj_p2v(pid)) {
    const Eigen::Map<const Eigen::Vector3d> position(trimesh.vert(vidx).ptr());
    V.col(num_col) = position;
    ++num_col;
  }
  return V;
}
bool pmp_jie::Vec3dIsNan(const Eigen::Vector3d &vert) {
  if(std::isnan(vert.x()) || std::isnan(vert.y()) || std::isnan(vert.z())) {
    return true;
  }
  return false;
}
bool pmp_jie::Vec3dIsInf(const Eigen::Vector3d &vert) {
  if(std::isinf(vert.x()) || std::isinf(vert.y()) || std::isinf(vert.z())) {
    return true;
  }
  return false;
}
Eigen::Vector2d pmp_jie::PlanRotate(const Eigen::Vector2d &vec,
                                    const double &rotate_angle,
                                    const bool &counter_click_wise) {
  Eigen::Matrix2d R;
  const auto cos_a = std::cos(rotate_angle);
  const auto sin_a = std::sin(rotate_angle);
  if(counter_click_wise) {
    R << cos_a, -sin_a, sin_a, cos_a;
  }
  else {
    R << cos_a, sin_a, -sin_a, cos_a;
  }
  return R*vec;
}
bool pmp_jie::IsDegenerate(const Eigen::Vector3d &vert) {
  if(vert(0)==0 && vert(1)==0 && vert(2)==0) return true;
  if(std::isnan(vert(0))) return true;
  if(std::isnan(vert(1))) return true;
  if(std::isnan(vert(2))) return true;
  if(std::isinf(vert(0))) return true;
  if(std::isinf(vert(1))) return true;
  if(std::isinf(vert(2))) return true;
  return false;
}
