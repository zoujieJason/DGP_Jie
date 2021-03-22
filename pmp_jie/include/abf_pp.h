//
// Created by jie zou on 2020/7/23.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_ABF_PP_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_ABF_PP_H_

#include "angle_reconstruct.h"
#include "linear_solver.h"
#include "mesh_io.h"
#include "tools.h"
#include "timer.h"

namespace pmp_jie {
class ABFPara: public CTriMesh {
 public:
  ABFPara();

  int SetVerbose(const bool &verbose);
  int ToAngleSpace(Eigen::VectorXd &alpha);
  int LinearABF(Eigen::VectorXd &alpha);

 protected:
  virtual int UsrAllocateVars() override;
  virtual int InitAlpha();
  virtual int InitBeta();
  virtual int InitLambda();
  virtual int InitWeight();
  int InitGamma();


  int OrganizeInverseLambda();
  int OrganizeInverseLambdaStar();
  int UpdateJorden();
  virtual int UpdateDerivatives();
  int UpdateDeLambda();
  double UpdateVars(const Eigen::VectorXd &delta_alpha, const Eigen::VectorXd &delta_lambda, const double &ratio);
  virtual int DeAlpha();
  virtual int DeLambdaCtri();
  virtual int DeLambdaCplan();
  virtual int DeLambdaClen();
  virtual size_t SetCtri(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const;
  virtual size_t SetCplan(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const;
  virtual size_t SetClen(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const;
  //TODO: Debug linear ABF.
  size_t SetLinearClen(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const;
  size_t SetLinearSolverConstants(Eigen::VectorXd &b) const;

  double ProdSinVoplus1(const size_t &vid) const;
  double ProdSinVominus1(const size_t &vid) const;

  virtual double SolveCurrIteration();
  virtual double CalStepLength(const Eigen::VectorXd &delta);

  double DFVal() const;
  double FVal() const;
  virtual double TargetEnergy() const;
  virtual double Ctri() const;
  virtual double Cplan() const;
  virtual double Clen() const;

  //matrix
  Eigen::SparseMatrix<double> invL_;
  Eigen::SparseMatrix<double> invL_star_;
  Eigen::SparseMatrix<double> Jt_;

  //angles
  Eigen::VectorXd alpha_;
  Eigen::VectorXd beta_;
  Eigen::VectorXd gamma_;
  Eigen::VectorXd optimal_angle_weight_;

  //lambda
  Eigen::VectorXd lctri_;
  Eigen::VectorXd lcplan_;
  Eigen::VectorXd lclen_;

  //derivatives
  Eigen::VectorXd dalpha_;
  Eigen::VectorXd dlctri_;
  Eigen::VectorXd dlcplan_;
  Eigen::VectorXd dlclen_;
  Eigen::VectorXd dlambda_;

  //coefficients
  double epsilon_{1e-5};
  bool verbose_{false};
  size_t num_ctri_{0};
  size_t num_cplan_{0};
  size_t num_clen_{0};
  size_t num_constrains_{0};

 private:
  double GetBoundaryVertAngle(const size_t &vid);
};
}
#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_ABF_PP_H_
