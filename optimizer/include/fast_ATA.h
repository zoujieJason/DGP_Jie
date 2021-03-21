//
// Created by jie zou on 2020/12/18.
//

#ifndef DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_FAST_ATA_H_
#define DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_FAST_ATA_H_

#include <stdint.h>
#include <algorithm>
#include <iostream>

#include <Eigen/Sparse>

namespace opt {
//! @param AAT must be a sorted csc with correct patten
//! @param method: which_part: -1 lower, 0 full, 1 upper, -1 is not
//! supported currently.
void fast_AAT(const Eigen::SparseMatrix<double> &A,
              Eigen::SparseMatrix<double> &AAT,
              int which_part = 0);

void fast_ADAT(const Eigen::SparseMatrix<double> &A,
               const Eigen::VectorXd &D,
               Eigen::SparseMatrix<double> &AAT,
               int which_part = 0);

}

#endif //DIGITAL_GEOMETRY_PROCESSING_SOLVER_INCLUDE_FAST_ATA_H_
