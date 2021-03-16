//
// Created by jie zou on 2021/3/16.
//
#include "../include/tools.h"
int pmp_jie::write_uv_to_obj(const Eigen::MatrixXd &uv,
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
