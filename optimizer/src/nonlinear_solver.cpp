//
// Created by jie zou on 2020/12/14.
//
#include "../include/nonlinear_solver.h"
opt::LM_Gauss_Newton::LM_Gauss_Newton(const opt::Function &func) {
  const auto nx(func.dim_of_x()), nf(func.dim_of_f());
  JT_.resize(nx,nf);
  g_.setZero(nx);
  h_.setZero(nx);
}
int opt::LM_Gauss_Newton::solve(double *x_ptr, opt::Function &func, int &num_iter) {
  const double min_max_mu[] = {1e-5, 1e6};
  const double eps(1e-8);
  double mu(1e-3), nu(2.0);
  Eigen::Map<Eigen::VectorXd> x0(x_ptr, JT_.rows());
  Eigen::VectorXd f0(JT_.cols());
  func.val(&x0[0],&f0[0]);
  double Fx = func.energy(&x0[0]);
  for(size_t i = 0; i < num_iter; ++i) {
    compute_grad(&x0[0],func,&f0[0]);
    const auto cond1 = g_.lpNorm<Eigen::Infinity>();

    std::cout << "current iter(" << i << "). ";

    if(cond1< eps) {
      std::cout << "Condition 1 (|g|_inf) reached." << std::endl;
      num_iter = i;
      break;
    }
    solve_step(&h_[0], mu);
    const auto cond2 = h_.norm() - eps*(x0.norm()+eps);

    std::cout << "var: "
              << "mu("      << mu << "), "
              << "nu("      << nu << "), "
              << "|g|_inf(" << cond1 << "), "
              << "|h|("     << cond2 << "). ";

    if(cond2 < 0) {
      x0 = x0 + h_;
      std::cout << "Condition 2 (|h|) reached. " << std::endl;
      num_iter = i;
      break;
    }
    const double Fx_new          = func.energy(&Eigen::VectorXd(x0 + h_)[0]);
    const double actual_residual = Fx - Fx_new;
    const double linear_residual = 0.5*h_.transpose()*(mu*h_-g_);
    const double varrho          = actual_residual/linear_residual;
    if(varrho > 0) {
      std::cout << "Acceptable iter. " << std::endl;
      x0 = x0 + h_;
      mu *= std::max(1.0/3.0,1.0-std::pow(2.0*varrho-1,3));
      nu = 2;
      Fx = Fx_new;
    }
    else {
      std::cout << "Unacceptable iter. " << std::endl;
      mu *= nu;
      nu *= 2;
    }
    if(mu<min_max_mu[0]) mu = min_max_mu[0];
    if(mu>min_max_mu[1]) mu = min_max_mu[1];
  }
  return 0;
}
int opt::LM_Gauss_Newton::compute_grad(const double *x_ptr, opt::Function &func, double *f_ptr) {
  func.jac(x_ptr, JT_);
  func.val(x_ptr, f_ptr);

  Eigen::Map<Eigen::VectorXd> f(f_ptr, JT_.cols());
  g_ = JT_*f;
  return 0;
}
int opt::LM_Gauss_Newton::get_diag_ptr() {
  if(!JTJ_.nonZeros())
    return 0;
  const int num_diag = JTJ_.rows();
  diag_ptr_.resize(num_diag);
  for(int i = 0; i < num_diag; ++i) {
    diag_ptr_[i] = &JTJ_.coeffRef(i, i);
  }
  return 1;
}
int opt::LM_Gauss_Newton::solve_step(double *h_ptr, const double &mu) {
  clock_t beg = clock();
  if(JTJ_.nonZeros())
    opt::fast_AAT(JT_, JTJ_);
  else {
    JTJ_ = JT_ * JT_.transpose();
    get_diag_ptr();
    slv_.analyzePattern(JTJ_);
  }
  std::cout << "# JTJ time: " << double(clock()-beg)/CLOCKS_PER_SEC << std::endl;

  for(auto &ptr: diag_ptr_) {
    (*ptr) += mu;
  }

  slv_.factorize(JTJ_);
  if(slv_.info()!= Eigen::Success) {
    std::cerr << "solve_step err: matrix decomposition failed!" << std::endl;
    return __LINE__;
  }
  Eigen::Map<Eigen::VectorXd> h(h_ptr,JTJ_.rows());
  h = slv_.solve(-g_);

  return 0;
}
opt::More_LM_Gauss_Newton::More_LM_Gauss_Newton(const opt::Function &func) : LM_Gauss_Newton(func) {
  diag_.setZero(func.dim_of_x());
}
int opt::More_LM_Gauss_Newton::solve(double *x_ptr, opt::Function &func, int &num_iter) {
  const double min_max_mu[] = {1e-5, 1e6};
  const double eps(1e-8);
  double mu(1e-3), nu(2.0);
  Eigen::Map<Eigen::VectorXd> x0(x_ptr, JT_.rows());
  Eigen::VectorXd f0(JT_.cols());
  func.val(&x0[0],&f0[0]);
  double Fx = func.energy(&x0[0]);
  for(size_t i = 0; i < num_iter; ++i) {
    compute_grad(&x0[0],func,&f0[0]);
    const auto cond1 = g_.lpNorm<Eigen::Infinity>();

    std::cout << "current iter(" << i << "). ";

    if(cond1< eps) {
      std::cout << "Condition 1 (|g|_inf) reached." << std::endl;
      num_iter = i;
      break;
    }
    if(solve_step(&h_[0], mu)) {
      mu *= 4.0;
      continue;
    }

    const auto cond2 = h_.norm() - eps*(x0.norm()+eps);

    std::cout << "var: "
              << "mu("      << mu << "), "
              << "|g|_inf(" << cond1 << "), "
              << "|h|("     << cond2 << "). ";

    if(cond2 < 0) {
      x0 = x0 + h_;
      std::cout << "Condition 2 (|h|) reached. " << std::endl;
      num_iter = i;
      break;
    }
    const double Fx_new          = func.energy(&Eigen::VectorXd(x0 + h_)[0]);
    const double actual_residual = Fx - Fx_new;
    const double linear_residual = 0.5*h_.transpose()*(mu*h_-g_);
    const double varrho          = actual_residual/linear_residual;
    if(varrho > 0) {
      std::cout << "Acceptable iter. " << std::endl;
      x0 = x0 + h_;
      mu *= std::max(1.0/3.0,1.0-std::pow(2.0*varrho-1,3));
      nu = 2;
      Fx = Fx_new;
    }
    else {
      std::cout << "Unacceptable iter. " << std::endl;
      mu *= nu;
      nu *= 2;
    }
    if(mu<min_max_mu[0]) mu = min_max_mu[0];
    if(mu>min_max_mu[1]) mu = min_max_mu[1];
  }
  return 0;
}
int opt::More_LM_Gauss_Newton::solve_step(double *h_ptr, const double &mu) {
  clock_t beg = clock();
  if(JTJ_.nonZeros())
    opt::fast_AAT(JT_, JTJ_);
  else {
    JTJ_ = JT_ * JT_.transpose();
    get_diag_ptr();
    slv_.analyzePattern(JTJ_);
  }
  std::cout << "# JTJ time: " << double(clock()-beg)/CLOCKS_PER_SEC << std::endl;

  compute_diag(&diag_[0]);
  {
    const Eigen::VectorXd diag = mu*diag_;
    set_JTJ_diag(&diag[0]);
  }

  slv_.factorize(JTJ_);
  if(slv_.info()!= Eigen::Success) {
    std::cerr << "solve_step err: matrix decomposition failed!" << std::endl;
    return __LINE__;
  }
  Eigen::Map<Eigen::VectorXd> h(h_ptr,JTJ_.rows());
  h = slv_.solve(-g_);
  return 0;
}
void opt::More_LM_Gauss_Newton::set_JTJ_diag(const double *diag_ptr) {
//  Eigen::Map<const Eigen::VectorXd> diag(diag_ptr, diag_.size());
  for(size_t i = 0; i < diag_ptr_.size(); ++i) {
    (*diag_ptr_[i]) = diag_ptr[i];
  }
}
void opt::More_LM_Gauss_Newton::compute_diag(double *diag) {
  Eigen::VectorXd D1 = Eigen::VectorXd::Zero(diag_.size());
  for(size_t ci = 0; ci < JT_.cols(); ++ci) {
    for(Eigen::SparseMatrix<double>::InnerIterator it(JT_,ci); it; ++it){
      D1[it.index()] += it.value() * it.value();
    }
//        for(size_t nzi = JT_.outerIndexPtr()[ci]; nzi < JT_.outerIndexPtr()[ci+1]; ++nzi) {
//            D1[JT_.innerIndexPtr()[nzi]] += JT_.valuePtr()[nzi]*JT_.valuePtr()[nzi];
//        }
  }
  D1 = D1.cwiseSqrt();
  // STRANGE: should be squared, DTD, according to
  // the book, but removing the sqrt leads to slow
  // convergence.
  for(ptrdiff_t i = 0; i < JT_.rows(); ++i) {
    if(diag[i] < D1[i])
      diag[i] = D1[i];
  }
}
opt::Dog_Leg::Dog_Leg(const opt::Function &func) {
  const auto nx(func.dim_of_x()), nf(func.dim_of_f());
  JT_.resize(nx,nf);
  g_.setZero(nx);
  h_.setZero(nx);
}
int opt::Dog_Leg::solve(double *x_ptr, opt::Function &func, int &num_iter) {
  const double eps(1e-8);
  Eigen::Map<Eigen::VectorXd> x0(x_ptr, JT_.rows());
  Eigen::VectorXd f0(JT_.cols());
  func.val(&x0[0],&f0[0]);
  double delta(1.0), Fx(func.energy(&x0[0])), alpha(0.0);
  for(size_t i = 0; i < num_iter; ++i) {
    compute_grad(&x0[0],func,&f0[0]);
    const auto cond1 = g_.lpNorm<Eigen::Infinity>();
    const auto cond2 = f0.lpNorm<Eigen::Infinity>();
    const auto cond3 = delta - eps*(x0.norm()+eps);

    alpha = g_.squaredNorm()/(JT_.transpose()*g_).squaredNorm();

    std::cout << "current iter(" << i << "). ";

    std::cout << "|g|_inf( " << cond1 << "), "
              << "|f|_inf( " << cond2 << "), "
              << "delta( " << cond3 << "), ";

    if(cond1<eps) {
      std::cout << "Condition 1 (|g|_inf) reached." << std::endl;
      break;
    }
    if(cond2<eps) {
      std::cout << "Condition 2 (|f|_inf) reached." << std::endl;
      break;
    }
    if(cond3<=0) {
      std::cout << "Condition 3 (delta) reached." << std::endl;
      break;
    }
    Eigen::VectorXd h_sd = -alpha*g_;
    Eigen::VectorXd h_gn(h_.size());
    solve_step(&h_gn[0], 0.0);
    const double linear_residual = dog_leg_step(alpha, &h_gn[0], &h_[0], delta, Fx);
    const double Fx_new          = func.energy(&Eigen::VectorXd(x0 + h_)[0]);
    const double actual_residual = Fx - Fx_new;
    const double ratio = actual_residual/linear_residual;
    const double cond4 = h_.norm()-eps*(x0.norm()-eps);
    std::cout << "|h|_dl( " << cond4 << " ). ";
    if(cond4<=0) {
      x0 = x0 + h_;
      std::cout << "Condition 4 (|h_dl|) reached." << std::endl;
      break;
    }
    //update delta
    if(ratio > 0.75) {
      std::cout << "\033[32m Increase delta. \033[0m";
      delta = std::max(delta, 3.0 * h_.norm());
    }
    else if(ratio < 0.25) {
      std::cout << "\033[31m Reduce delta. \033[0m";
      delta /= 2.0;
    }
    else {
      std::cout << "\033[30m Delta unchanged. \033[0m";
    }
    if(ratio > 0) {
      std::cout << "\033[32m Acceptable iter. \033[0m" << std::endl;
      x0 = x0 + h_;
      Fx = Fx_new;
    }
    else {
      std::cout << "\033[31m Unacceptable iter. \033[0m" << std::endl;
    }
  }
  return 0;
}
int opt::Dog_Leg::compute_grad(const double *x_ptr, opt::Function &func, double *f_ptr) {
  func.jac(x_ptr, JT_);
  func.val(x_ptr, f_ptr);

  Eigen::Map<Eigen::VectorXd> f(f_ptr, JT_.cols());
  g_ = JT_*f;
  return 0;
}
int opt::Dog_Leg::solve_step(double *h_ptr, const double &mu) {
  clock_t beg = clock();
  if(JTJ_.nonZeros())
    opt::fast_AAT(JT_, JTJ_);
  else {
    JTJ_ = JT_ * JT_.transpose();
    slv_.analyzePattern(JTJ_);
  }
  std::cout << "# JTJ time: " << double(clock()-beg)/CLOCKS_PER_SEC << std::endl;

  slv_.factorize(JTJ_);
  if(slv_.info()!= Eigen::Success) {
    std::cerr << "solve_step err: matrix decomposition failed!" << std::endl;
    return __LINE__;
  }
  Eigen::Map<Eigen::VectorXd> h_gn(h_ptr, h_.size());
  h_gn = slv_.solve(-g_);
  return 0;
}
double opt::Dog_Leg::dog_leg_step(const double &alpha,
                                       const double *h_gn_ptr,
                                       double *h_ptr,
                                       const double &delta,
                                       const double &Fx) {
  Eigen::Map<const Eigen::VectorXd> h_gn(h_gn_ptr, h_.size());
  Eigen::VectorXd h_sd = -g_;
  Eigen::Map<Eigen::VectorXd> h(h_ptr, h_.size());

  if(h_gn.norm()<=delta) {
    h = h_gn;
    return Fx;
  }
  else if((alpha*h_sd).norm()>=delta) {
    h = delta * h_sd.normalized();
    return (delta*(2.0*((alpha*g_).norm())-delta))/(2.0*alpha);
  }
  else {
    const auto a = alpha*h_sd;
    const auto b_a = h_gn-a;
    const double c = a.dot(b_a);
    double beta(0.0);
    if(c<=0) {
      beta =(-c+std::sqrt(c*c+b_a.squaredNorm()*(delta*delta-a.squaredNorm())))/b_a.squaredNorm();
    }
    else {
      beta = (delta*delta-a.squaredNorm())/(c+std::sqrt(c*c+b_a.squaredNorm()*(delta*delta-a.squaredNorm())));
    }
    return 0.5*alpha*(1.0-beta*beta)*g_.squaredNorm()+beta*(2.0-beta)*Fx;
  }
}
