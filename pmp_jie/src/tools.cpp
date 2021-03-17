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
