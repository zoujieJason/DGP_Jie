//
// Created by jie zou on 2021/3/16.
//

#ifndef DGP_JIE_OPTIMIZER_INCLUDE_LINEAR_OPTIMIZER_H_
#define DGP_JIE_OPTIMIZER_INCLUDE_LINEAR_OPTIMIZER_H_

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace opt {
int least_squares_solve_ATA(const Eigen::SparseMatrix<double> &A,
                            const Eigen::MatrixXd &b,
                            Eigen::MatrixXd &results);
}

#endif //DGP_JIE_OPTIMIZER_INCLUDE_LINEAR_OPTIMIZER_H_
