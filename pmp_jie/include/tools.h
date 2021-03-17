//
// Created by jie zou on 2021/3/16.
//

#ifndef DGP_JIE_PMP_JIE_INCLUDE_TOOLS_H_
#define DGP_JIE_PMP_JIE_INCLUDE_TOOLS_H_

#include <vector>
#include <iostream>
#include <fstream>

#include <Eigen/Dense>

namespace pmp_jie {
template<class search_t>
int LinearSearch(const std::vector<search_t> &vec, const search_t &target) {
  for(size_t i = 0; i < vec.size(); ++i) {
    if(vec[i] == target) {
      return i;
    }
  }
  return -1;
}
template<class search_t>
int LinearSearchOffset(const std::vector<search_t> &vec, const search_t &target) {
  for(size_t i = 0; i < vec.size(); ++i) {
    if(target == vec[i]) {
      return -1;
    }
    else if(target < vec[i]){
      return i;
    }
  }
  return vec.size();
}

int WriteUVToObj(const Eigen::MatrixXd &uv,
                 const std::vector<std::vector<uint>> &polys,
                 const std::string &filename);
template<class search_t>
void Swap(search_t &lhs, search_t &rhs) {
  search_t tmp = lhs;
  lhs = rhs;
  rhs = tmp;
}
std::string GetFilePath(const std::string &string);
std::string GetFileName(const std::string &string, const bool &with_extension);
std::string GetFileExtension(const std::string &string);
double GaussValue(const double &x, const double &sigma);
}
#endif //DGP_JIE_PMP_JIE_INCLUDE_TOOLS_H_
