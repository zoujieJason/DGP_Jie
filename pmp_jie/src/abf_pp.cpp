//
// Created by jie zou on 2020/7/23.
//
#include "../include/abf_pp.h"

pmp_jie::ABFPara::ABFPara()=default;
int pmp_jie::ABFPara::UsrAllocateVars() {
  num_ctri_ = num_polys_;
  num_cplan_ = num_int_verts_;
  num_clen_ = num_int_verts_;
  num_constrains_ = num_ctri_+num_cplan_+num_clen_;
  InitBeta();
  InitWeight();
  InitAlpha();
  InitLambda();
  OrganizeInverseLambda();
  OrganizeInverseLambdaStar();
  UpdateJorden();
  UpdateDerivatives();
  return 0;
}
int pmp_jie::ABFPara::InitAlpha() {
  alpha_.setZero(num_angles_);
  alpha_ = beta_;
  return 0;
}
int pmp_jie::ABFPara::InitBeta() {
  beta_.resize(num_angles_);
  beta_.setZero();
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      double ratio = 1.0;
      if(!trimesh_.vert_is_boundary(adj_vid)) {
        const double sum_angle = SumVertAngle(adj_vid);
        ratio = 2 * M_PI / sum_angle;
      }
      const auto aid = GetAngleId(pid, adj_vid);
      beta_(aid) = trimesh_.poly_angle_at_vert(pid, adj_vid) * ratio;
      beta_(aid) = std::max(beta_(aid), 3.0 * M_PI / 180.0);
      beta_(aid) = std::min(beta_(aid), 175.0 * M_PI / 180.0);
    }
  }
  return 0;
}
int pmp_jie::ABFPara::InitWeight() {
  optimal_angle_weight_.setZero(num_angles_);
  for(size_t i = 0; i < beta_.size(); ++i) {
    optimal_angle_weight_(i) = 1.0 / (beta_(i) * beta_(i));
  }
  return 0;
}
int pmp_jie::ABFPara::OrganizeInverseLambda() {
  std::vector<Eigen::Triplet<double>> triplets;
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    //set triplet for each angle of triangles.
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      const auto aid = GetAngleId(pid, adj_vid);
      triplets.emplace_back(aid, aid, 1.0 / (2.0 * optimal_angle_weight_(aid)));
    }
  }
  invL_.resize(num_angles_, num_angles_ );
  invL_.setZero();
  invL_.setFromTriplets(triplets.begin(), triplets.end());
  return 0;
}
int pmp_jie::ABFPara::UpdateJorden() {
  std::vector<Eigen::Triplet<double>> triplets;
  SetCtri(0,triplets);
  SetCplan(num_ctri_, triplets);
  SetClen(num_ctri_+num_cplan_, triplets);
  Jt_.resize(num_angles_,num_constrains_);
  Jt_.setZero();
  Jt_.setFromTriplets(triplets.begin(), triplets.end());
  return 0;
}
size_t pmp_jie::ABFPara::SetCtri(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const {
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    for(size_t i = 0; i < 3; ++i) {
      triplets.emplace_back(3 * pid + i, curr_cols + pid, 1.0);
    }
  }
  return num_polys_;
}
size_t pmp_jie::ABFPara::SetCplan(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const {
  size_t int_vid(0);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
        triplets.emplace_back(GetAngleId(adj_pid, vid), curr_cols + int_vid, 1.0);
      }
      ++int_vid;
    }
  }
  return int_vid;
}
size_t pmp_jie::ABFPara::SetClen(const size_t &curr_cols, std::vector<Eigen::Triplet<double>> &triplets) const {
  size_t int_vid(0);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      const auto prod_voplus1_sin  = ProdSinVoplus1(vid);
      const auto prod_vominus1_sin = -1.0 * ProdSinVominus1(vid);
      for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
        const auto voplus1 = GetPolyCCWVert(adj_pid, vid);
        assert(voplus1 != -1);
        const auto voplus1_radian = alpha_(GetAngleId(adj_pid, voplus1));
        triplets.emplace_back(GetAngleId(adj_pid, voplus1), curr_cols + int_vid, prod_voplus1_sin * std::cos(voplus1_radian) / std::sin(voplus1_radian));

        const auto vominus1 = GetPolyCCWVert(adj_pid, voplus1);
        assert(vominus1 != -1);
        const auto vominus1_radian = alpha_(GetAngleId(adj_pid, vominus1));
        triplets.emplace_back(GetAngleId(adj_pid, vominus1), curr_cols + int_vid, prod_vominus1_sin * std::cos(vominus1_radian) / std::sin(vominus1_radian));
      }
      ++int_vid;
    }
  }
  return int_vid;
}
double pmp_jie::ABFPara::ProdSinVoplus1(const size_t &vid) const {
  double prod(1.0);
  for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
    const auto voplus1 = GetPolyCCWVert(adj_pid, vid);
    assert(voplus1 != -1);
    const auto voplus1_radian = alpha_(GetAngleId(adj_pid, voplus1));
    prod *= std::sin(voplus1_radian);
  }
  return prod;
}
double pmp_jie::ABFPara::ProdSinVominus1(const size_t &vid) const {
  double prod(1.0);
  for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
    const auto voplus1 = GetPolyCCWVert(adj_pid, vid);
    assert(voplus1 != -1);
    const auto vominus1 = GetPolyCCWVert(adj_pid, voplus1);
    assert(vominus1 != -1);
    const auto vominus1_radian = alpha_(GetAngleId(adj_pid, vominus1));
    prod *= std::sin(vominus1_radian);
  }
  return prod;
}
int pmp_jie::ABFPara::UpdateDerivatives() {
  DeAlpha();
  DeLambdaCtri();
  DeLambdaCplan();
  DeLambdaClen();
  UpdateDeLambda();
  return 0;
}
int pmp_jie::ABFPara::DeAlpha() {
  dalpha_.setZero(num_angles_);
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      const auto aid   = GetAngleId(pid, adj_vid);
      const double energy   = 2.0 * (alpha_(aid) - beta_(aid)) * optimal_angle_weight_(aid);
      dalpha_(aid) -= energy;
      dalpha_(aid) -= lctri_(pid);
    }
  }
  size_t int_vid = 0;
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      const auto prod_next_sin     = ProdSinVoplus1(vid);
      const auto prod_prev_sin     = ProdSinVominus1(vid);
      for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
        const auto voplus1           = GetPolyCCWVert(adj_pid, vid);
        const auto vominus1          = GetPolyCCWVert(adj_pid, voplus1);
        const double voplus1_radian  =  alpha_(GetAngleId(adj_pid, voplus1));
        const double vominus1_radian =  alpha_(GetAngleId(adj_pid, vominus1));
        dalpha_(GetAngleId(adj_pid, vid))      -= lcplan_(int_vid);
        dalpha_(GetAngleId(adj_pid, voplus1))  -= lclen_(int_vid) * prod_next_sin  * std::cos(voplus1_radian) / std::sin(voplus1_radian);
        dalpha_(GetAngleId(adj_pid, vominus1)) += lclen_(int_vid) * prod_prev_sin  * std::cos(vominus1_radian) / std::sin(vominus1_radian);
      }
      ++int_vid;
    }
  }
  return 0;
}
int pmp_jie::ABFPara::DeLambdaCtri() {
  dlctri_.setZero(num_ctri_);
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    double sum_angles(0.0);
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      sum_angles += alpha_(GetAngleId(pid, adj_vid));
    }
    const double diff = sum_angles - M_PI;
    dlctri_(pid) -= diff;
  }
  return 0;
}
int pmp_jie::ABFPara::InitLambda() {
  lctri_.setZero(num_ctri_);
  lcplan_.setZero(num_cplan_);
  lclen_.setZero(num_clen_);
  dlambda_.setZero(num_constrains_);
  return 0;
}
int pmp_jie::ABFPara::DeLambdaCplan() {
  dlcplan_.setZero(num_cplan_);
  size_t int_vid(0);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      double sum_angles(0.0);
      for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
        sum_angles += alpha_(GetAngleId(adj_pid, vid));
      }
      const double offset = sum_angles - 2 * M_PI;
      dlcplan_(int_vid) -= offset;
      ++int_vid;
    }
  }
  assert(int_vid == num_int_verts_);
  return int_vid;
}
int pmp_jie::ABFPara::DeLambdaClen() {
  dlclen_.setZero(num_clen_);
  size_t int_vid(0);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      const double clen = ProdSinVoplus1(vid) - ProdSinVominus1(vid);
      dlclen_(int_vid) -= clen;
      ++int_vid;
    }
  }
  assert(int_vid == num_int_verts_);
  return int_vid;
}
int pmp_jie::ABFPara::ToAngleSpace(Eigen::VectorXd &alpha) {
  double curr_dfval = DFVal();
  double curr_fval = FVal();
  double delta(0.0);
  std::cout << "INIT( dfval: " << curr_dfval << " fval: " << curr_fval << " )"<< std::endl;
  size_t it_cnt(0);
  do {
    delta = SolveCurrIteration();
    OrganizeInverseLambda();
    OrganizeInverseLambdaStar();
    UpdateJorden();
    UpdateDerivatives();
    curr_dfval = DFVal();
    curr_fval = FVal();
    if(verbose_) {
      std::cout << "iteration(" << it_cnt << "): ";
      std::cout << "    dfval: " << curr_dfval;
      std::cout << "         delta: " << delta;
      std::cout << "         fval : " << curr_fval <<std::endl;
    }
    ++it_cnt;
  } while (curr_dfval > epsilon_ && it_cnt < 100 && delta > epsilon_);
  if(curr_dfval <= epsilon_) {
    std::cout << "Gradient converge." << std::endl;
  }
  if(delta <= epsilon_) {
    std::cout << "Delta converge." << std::endl;
  }
  std::cout << "LAST ITERATION( dfval: " << curr_dfval << " delta: " << delta << " )" << std::endl;
  alpha.setZero(num_angles_);
  alpha = alpha_;
  return 0;
}
double pmp_jie::ABFPara::DFVal() const {
  return dalpha_.lpNorm<1>() + dlctri_.lpNorm<1>() + dlcplan_.lpNorm<1>() + dlclen_.lpNorm<1>();
}
double pmp_jie::ABFPara::FVal() const {
  return TargetEnergy()+Ctri()+Cplan()+Clen();
}
double pmp_jie::ABFPara::TargetEnergy() const {
  double sum_energy(0.0);
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    //set triplet for each angle of triangles.
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      const auto curr_aid = GetAngleId(pid, adj_vid);
      const auto alpha_radian = alpha_(curr_aid);
      const auto beta_radian  = beta_(curr_aid);
      sum_energy += (alpha_radian - beta_radian) * (alpha_radian - beta_radian) * optimal_angle_weight_(curr_aid);
    }
  }
  return sum_energy;
}
double pmp_jie::ABFPara::Ctri() const {
  double sum_energy(0.0);
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    double sum_angles = 0.0;
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      sum_angles += alpha_(GetAngleId(pid, adj_vid));
    }
    sum_energy += lctri_(pid) * (sum_angles - M_PI);
  }
  return sum_energy;
}
double pmp_jie::ABFPara::Cplan() const {
  double sum_energy(0.0);
  size_t int_vid(0);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      double sum_angles = 0.0;
      for(const auto &adj_pid: trimesh_.adj_v2p(vid)) {
        sum_angles += alpha_(GetAngleId(adj_pid, vid));
      }
      sum_energy += lcplan_(int_vid) * (sum_angles - 2 * M_PI);
      ++int_vid;
    }
  }
  return sum_energy;
}
double pmp_jie::ABFPara::Clen() const {
  double sum_energy(0.0);
  size_t int_vid(0);
  for(size_t vid = 0; vid < num_verts_; ++vid) {
    if(!trimesh_.vert_is_boundary(vid)) {
      const auto prod_voplus1  = ProdSinVoplus1(vid);
      const auto prod_vominus1 = ProdSinVominus1(vid);
      sum_energy += lclen_(int_vid) * (prod_voplus1 - prod_vominus1);
      ++int_vid;
    }
  }
  return sum_energy;
}
double pmp_jie::ABFPara::SolveCurrIteration() {
  //J_t = [J1, J2, J3]
  //decompose J_2 = [J1, J2]
  Eigen::SparseMatrix<double> J1t_;
  J1t_.resize(num_angles_, num_ctri_);
  J1t_ = Jt_.leftCols(num_ctri_);
  Eigen::SparseMatrix<double> J2t_;
  J2t_.resize(num_angles_, num_cplan_ + num_clen_);
  J2t_ = Jt_.rightCols(num_cplan_ + num_clen_);

  //[Lambda-1* J*^t]
  //[    J*    J** ]
  const auto J_star_t    = (J2t_.transpose() * invL_ * J1t_).transpose();
  const auto J_star_star = J2t_.transpose() * invL_ * J2t_;

  //b1 = dalpha
  //b2 = dlambda = [dlambda_ctri, dlambda_cplan, dlambda_clen]
  //
  const auto b_star   = Jt_.transpose() * invL_ * dalpha_ - dlambda_;
  const auto b_star_1 = b_star.head(dlctri_.size());
  const auto b_star_2 = b_star.tail(dlcplan_.size() + dlclen_.size());

  Eigen::VectorXd delta_lambda1, delta_lambda2;
  opt::LinearBlockSolve(invL_star_, J_star_t, J_star_star, b_star_1, b_star_2, delta_lambda1, delta_lambda2);

  Eigen::VectorXd delta_lambda;
  delta_lambda.resize(delta_lambda1.size() + delta_lambda2.size());
  delta_lambda.head(delta_lambda1.size()) = delta_lambda1;
  delta_lambda.tail(delta_lambda2.size()) = delta_lambda2;
  const auto delta_alpha = invL_ * (dalpha_ - Jt_ * delta_lambda);

  const double ratio = CalStepLength(delta_alpha);
  return UpdateVars(delta_alpha, delta_lambda, ratio);
}
int pmp_jie::ABFPara::OrganizeInverseLambdaStar() {
  std::vector<Eigen::Triplet<double>> triplets;
  for(size_t pid = 0; pid < num_polys_; ++pid) {
    double sum(0.0);
    for(const auto &adj_vid: trimesh_.adj_p2v(pid)) {
      const auto aid = GetAngleId(pid, adj_vid);
      sum += 1.0 / (2.0 * optimal_angle_weight_(aid));
    }
    triplets.emplace_back(pid, pid, 1.0 / sum);
  }
  invL_star_.resize( num_polys_, num_polys_ );
  invL_star_.setZero();
  invL_star_.setFromTriplets(triplets.begin(), triplets.end());
  return 0;
}
int pmp_jie::ABFPara::UpdateDeLambda() {
  dlambda_.setZero(num_constrains_);
  dlambda_.head(num_ctri_) = dlctri_;
  dlambda_.middleRows(num_ctri_, num_cplan_) = dlcplan_;
  dlambda_.tail(num_clen_) = dlclen_;
  return 0;
}
double pmp_jie::ABFPara::CalStepLength(const Eigen::VectorXd &delta) {
  double ratio(1.0);
  for(size_t aid = 0; aid < alpha_.size(); ++aid) {
    if(alpha_(aid) + delta(aid) < 10.0 * epsilon_) {
      double r1 = -.5 * (alpha_(aid) - 10.0 * epsilon_)/delta(aid);
      ratio = std::min(ratio, r1) ;
      optimal_angle_weight_(aid) *= 1.2;
    } else if(alpha_(aid) + delta(aid) > M_PI - 10.0 * epsilon_) {
//	   double r1 = .5*(M_PI - alpha_(aid)+10.0 * epsilon_)/dalpha_(aid);
//	   ratio = std::min(ratio, r1);
      optimal_angle_weight_(aid) *= 1.2;
    }
  }
  return ratio;
}
double pmp_jie::ABFPara::UpdateVars(const Eigen::VectorXd &delta_alpha, const Eigen::VectorXd &delta_lambda, const double &ratio) {
  alpha_  = alpha_ + ratio * delta_alpha;
  lctri_  = lctri_ + ratio * delta_lambda.head(num_ctri_);
  lcplan_ = lcplan_ + ratio * delta_lambda.middleRows(num_ctri_, num_cplan_);
  lclen_  = lclen_ + ratio * delta_lambda.tail(num_clen_);
  return ratio * (delta_alpha.lpNorm<1>() + delta_lambda.lpNorm<1>());
}
int pmp_jie::ABFPara::SetVerbose(const bool &verbose) {
  verbose_ = verbose;
  return 0;
}

