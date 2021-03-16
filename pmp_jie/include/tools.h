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
int linear_search(const std::vector<search_t> &vec, const search_t &target) {
  for(size_t i = 0; i < vec.size(); ++i) {
    if(vec[i] == target) {
      return i;
    }
  }
  return -1;
}

int write_uv_to_obj(const Eigen::MatrixXd &uv,
                    const std::vector<std::vector<uint>> &polys,
                    const std::string &filename);
}
#endif //DGP_JIE_PMP_JIE_INCLUDE_TOOLS_H_
