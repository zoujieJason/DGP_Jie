//
// Created by jie zou on 2020/12/18.
//
#include "../include/fast_ATA.h"
void opt::fast_AAT(const Eigen::SparseMatrix<double> &A, Eigen::SparseMatrix<double> &AAT, int which_part) {
  using namespace std;

  if(which_part == -1) {
    cerr << "not support currently." << endl;
    return;
  }

  std::fill(AAT.valuePtr(), AAT.valuePtr() + AAT.nonZeros(), 0);
  assert(A.rows() == AAT.rows() && A.rows() == AAT.cols());
  const int n = A.cols();
  for(int i = 0;i < n;++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it_col(A,i); it_col; ++it_col){
      const int col_idx_of_AAT = it_col.index();
      const int nz_beg = AAT.outerIndexPtr()[col_idx_of_AAT];
      const int nz_end = AAT.outerIndexPtr()[col_idx_of_AAT+1];
      for (Eigen::SparseMatrix<double>::InnerIterator it_row(A,i); it_row; ++it_row){
        if(it_row.index() <= col_idx_of_AAT){
          const int nz_of_AAT =
              lower_bound(&AAT.innerIndexPtr()[nz_beg], AAT.innerIndexPtr()+nz_end, it_row.index()) - AAT.innerIndexPtr();
          AAT.valuePtr()[nz_of_AAT] += it_row.value()*it_col.value();
        }
      }
    }
  }

  // get another part
  for(int col_i = 0; col_i < AAT.cols() && which_part == 0; ++col_i) {
    for(Eigen::SparseMatrix<double>::InnerIterator it_col(AAT,col_i); it_col; ++it_col) {
      const int row_i = it_col.index();
      if(row_i > col_i) {
        if(it_col.value() != 0) {
          cerr << "error: " << it_col.value() << endl;
        }
        const int nz_end = AAT.outerIndexPtr()[row_i+1];
        const int offset = lower_bound(&AAT.innerIndexPtr()[AAT.outerIndexPtr()[row_i]],
                                       AAT.innerIndexPtr()+nz_end, col_i)-AAT.innerIndexPtr();
        it_col.valueRef() = AAT.valuePtr()[offset];
      }
    }
  }
}
void opt::fast_ADAT(const Eigen::SparseMatrix<double> &A,
                         const Eigen::VectorXd &D,
                         Eigen::SparseMatrix<double> &AAT,
                         int which_part) {
  using namespace std;

  if(which_part == -1) {
    cerr << "not support currently." << endl;
    return;
  }

  std::fill(AAT.valuePtr(), AAT.valuePtr() + AAT.nonZeros(), 0);
  assert(A.rows() == AAT.rows() && A.rows() == AAT.cols());
  const int n = A.cols();
  for(int i = 0;i < n;++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it_col(A,i); it_col; ++it_col){
      const int col_idx_of_AAT = it_col.index();
      const int nz_beg = AAT.outerIndexPtr()[col_idx_of_AAT];
      const int nz_end = AAT.outerIndexPtr()[col_idx_of_AAT+1];
      for (Eigen::SparseMatrix<double>::InnerIterator it_row(A,i); it_row; ++it_row){
        if(it_row.index() <= col_idx_of_AAT){
          const int nz_of_AAT =
              lower_bound(&AAT.innerIndexPtr()[nz_beg], AAT.innerIndexPtr()+nz_end, it_row.index()) - AAT.innerIndexPtr();
          AAT.valuePtr()[nz_of_AAT] += it_row.value()*it_col.value()*D[it_row.col()];
        }
      }
    }
  }

  // get another part
  for(int col_i = 0; col_i < AAT.cols() && which_part == 0; ++col_i) {
    for(Eigen::SparseMatrix<double>::InnerIterator it_col(AAT,col_i); it_col; ++it_col) {
      const int row_i = it_col.index();
      if(row_i > col_i) {
        if(it_col.value() != 0) {
          cerr << "error: " << it_col.value() << endl;
        }
        const int nz_end = AAT.outerIndexPtr()[row_i+1];
        const int offset = lower_bound(&AAT.innerIndexPtr()[AAT.outerIndexPtr()[row_i]],
                                       AAT.innerIndexPtr()+nz_end, col_i)-AAT.innerIndexPtr();
        it_col.valueRef() = AAT.valuePtr()[offset];
      }
    }
  }
}
