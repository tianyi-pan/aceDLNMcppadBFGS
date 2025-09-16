/** Include **/
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#include <random> // for generating samples from standard normal distribution


#include "aceDLNMcppadClass.hpp"
#include "aceDLNMeigenClass.hpp"

// #include <LBFGSB.h>
// using namespace LBFGSpp;
// TODO: use https://github.com/yixuan/LBFGSpp



struct LAMLResult {
    double fn;
    Eigen::VectorXd gradient;
};


CppAD::ADFun<double> LAMLTape(Model& modelobj, Modelcppad& modelcppadobj){
    // if(!modelcppadobj.ifhastape) {
      int kE = modelobj.kE;
      int kw = modelobj.kw;
      int kwopt = modelobj.kwopt;
      int kbetaR = modelobj.kbetaR;
      int kbetaF = modelobj.kbetaF;
      int p = modelobj.p;

      Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
      Eigen::VectorXd R_betaR = modelobj.betaR;
      Eigen::VectorXd R_betaF = modelobj.betaF;
      Eigen::VectorXd R_phi = modelobj.phi;

      double R_log_theta = modelobj.log_theta;
      double R_log_smoothing_f = modelobj.log_smoothing_f;
      double R_log_smoothing_w = modelobj.log_smoothing_w;
      Eigen::VectorXd R_logsmoothing = modelobj.logsmoothing;

      Vec alpha_f = R_alpha_f.cast<Scalar>();
      Vec betaR = R_betaR.cast<Scalar>();
      Vec betaF = R_betaF.cast<Scalar>();
      Vec phi = R_phi.cast<Scalar>();

      Scalar log_theta = R_log_theta;
      Scalar log_smoothing_f = R_log_smoothing_f;
      Scalar log_smoothing_w = R_log_smoothing_w;
      Vec logsmoothing = R_logsmoothing.cast<Scalar>();

      // alpha_f, phi, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing


      Vec at(kE+kbetaR+kbetaF + kwopt + 3 + p);
      at << alpha_f, phi, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing;


      // start reverse mode
      Vec result(1);
      CppAD::ADFun<double> gr;

      CppAD::Independent(at);
      modelcppadobj.setAlphaF(at.segment(0, kE));
      modelcppadobj.setPhi(at.segment(kE, kwopt));
      modelcppadobj.setBetaR(at.segment(kE+kwopt, kbetaR));
      modelcppadobj.setBetaF(at.segment(kE+kwopt+kbetaR, kbetaF));
      modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+kwopt));
      modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+kwopt+1));
      modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+kwopt+2));
      modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+kwopt+3, p));


      modelcppadobj.derivative_coef();
      modelcppadobj.derivative_he();
      result(0) = modelcppadobj.logdetH05();

      gr.Dependent(at, result);
      // gr.optimize();
      // modelcppadobj.gr = gr;
      return gr;
      // END reverse



      // Forward mode 1
      // std::vector<Scalar> at_std(at.data(), at.data() + (kE+kbetaR+kbetaF + kw-1 + 3 + p));

      // CppAD::Independent(at_std);

      // std::vector<Scalar> result(1);
      // Vec at_eigen = Eigen::Map<Vec>(at_std.data(), (kE+kbetaR+kbetaF + kw-1 + 3 + p));
      // modelcppadobj.setAlphaF(at_eigen.segment(0, kE));
      // modelcppadobj.setPhi(at_eigen.segment(kE, kw-1));
      // modelcppadobj.setBetaR(at_eigen.segment(kE+kw-1, kbetaR));
      // modelcppadobj.setBetaF(at_eigen.segment(kE+kw-1+kbetaR, kbetaF));
      // modelcppadobj.setLogTheta(at_eigen(kE+kbetaR+kbetaF+kw-1));
      // modelcppadobj.setLogSmoothingF(at_eigen(kE+kbetaR+kbetaF+kw-1+1));
      // modelcppadobj.setLogSmoothingW(at_eigen(kE+kbetaR+kbetaF+kw-1+2));
      // modelcppadobj.setLogsmoothing(at_eigen.segment(kE+kbetaR+kbetaF+kw-1+3, p));

      // modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      // result[0] = modelcppadobj.logdetH05();

      // CppAD::ADFun<double> gr(at_std, result);
      // END forward 1


      // Forward mode 2
      // CppAD::Independent(at);

      // Vec result(1);
      // modelcppadobj.setAlphaF(at.segment(0, kE));
      // modelcppadobj.setPhi(at.segment(kE, kw-1));
      // modelcppadobj.setBetaR(at.segment(kE+kw-1, kbetaR));
      // modelcppadobj.setBetaF(at.segment(kE+kw-1+kbetaR, kbetaF));
      // modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+kw-1));
      // modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+kw-1+1));
      // modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+kw-1+2));
      // modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+kw-1+3, p));


      // modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      // result(0) = modelcppadobj.logdetH05();

      // CppAD::ADFun<double> gr(at, result);
      // END forward 2


      // modelcppadobj.ifhastape = true;





    // }
}

LAMLResult LAML(Model& modelobj, Modelcppad& modelcppadobj) {
    LAMLResult result;
    double u_LAML;
    int kE = modelobj.kE;
    int kw = modelobj.kw;
    int kwopt = modelobj.kwopt;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int p = modelobj.p;
    modelobj.derivative_coef();
    modelobj.derivative_he();
    modelobj.derivative_full();
    modelobj.NegativeLogLikelihood();

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd phi = modelobj.phi;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;

    Eigen::VectorXd gr_s_u_vec = modelobj.gr_s_u_vec;
    Eigen::VectorXd gr_s_par_vec = modelobj.gr_s_par_vec;
    Eigen::MatrixXd he_s_u_mat = modelobj.he_s_u_mat;
    Eigen::MatrixXd he_s_par_u_mat = modelobj.he_s_par_u_mat;


    double log_theta = modelobj.log_theta;
    double log_smoothing_f = modelobj.log_smoothing_f;
    double log_smoothing_w = modelobj.log_smoothing_w;
    Eigen::VectorXd logsmoothing = modelobj.logsmoothing;

    // First derivative of LAML

    CppAD::ADFun<double> cppadgr = LAMLTape(modelobj, modelcppadobj);



    Eigen::VectorXd at0(kE+kbetaR+kbetaF + kwopt + 3 + p);
    at0 << alpha_f, phi, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing;

    // // reverse mode
    // Eigen::VectorXd g_LAML(kE + kwopt + 3 + kbetaR+kbetaF + p);
    // g_LAML.setZero();
    // g_LAML = cppadgr.Jacobian(at0);
    // // END reverse mode


    // using reverse mode Reverse(1).
    // Should be faster than Jacobian when dim(output) is 1.
    std::vector<double> x(at0.data(), at0.data() + at0.size());
    // 1) Evaluate f at x (zero-order forward)
    cppadgr.Forward(0, x);
    // 2) Seed the single output with 1.0 (since output is scalar)
    std::vector<double> w(1, 1.0);
    // 3) One reverse sweep gives the whole gradient
    std::vector<double> grad = cppadgr.Reverse(1, w);
    // 4) Map to Eigen
    Eigen::Map<const Eigen::VectorXd> g_LAML_map(grad.data(), static_cast<Eigen::Index>(grad.size()));
    Eigen::VectorXd g_LAML = g_LAML_map;
    // END using reverse mode Reverse(1)


    // forward mode 1, 2
    // std::vector<Scalar> at0_std(at0.data(), at0.data() + (kE+kbetaR+kbetaF + kw-1 + 3 + p));
    // std::vector<double> dx((kE+kbetaR+kbetaF + kw-1 + 3 + p), 0.0);
    // std::vector<double> g_LAML_std((kE+kbetaR+kbetaF + kw-1 + 3 + p));
    //  std::vector<double> dy(1);
    // for (size_t i = 0; i < (kE+kbetaR+kbetaF + kw-1 + 3 + p); ++i) {
    //     dx[i] = 1.0;
    //     dy = modelcppadobj.gr.Forward(1, dx);
    //     g_LAML_std[i] = dy[0];
    //     dx[i] = 0.0;
    // }

    // Eigen::VectorXd g_LAML = Eigen::Map<Eigen::VectorXd>(g_LAML_std.data(), g_LAML_std.size());
    // END forward mode 1, 2




    u_LAML = modelobj.logdetH05() + modelobj.NegLogL - modelobj.n/2.0 * log(2*3.141592653589793238462643383279);


    // In R: grad[-(1:(kE+kw-1))] - H.full[-(1:(kE+kw-1)),(1:(kE+kw-1))] %*% as.vector(solve(H.alpha, grad[(1:(kE+kw-1))]))
    Eigen::VectorXd g1 = g_LAML.segment(0, kE+kwopt+kbetaR+kbetaF) + gr_s_u_vec;
    Eigen::VectorXd g2 = g_LAML.segment(kE+kwopt+kbetaR+kbetaF, 3+p) + gr_s_par_vec;
    Eigen::VectorXd gr = g2 - he_s_par_u_mat * he_s_u_mat.ldlt().solve(g1);

    result.fn = u_LAML;
    result.gradient = gr;
    return result;
}




// [[Rcpp::export]]
List aceDLNMbuild(const Eigen::VectorXd R_y,
                   const Eigen::MatrixXd R_B_inner,
                   const Eigen::VectorXd R_knots_f,
                   const Eigen::MatrixXd R_Sw_large,
                   const Eigen::MatrixXd R_Sf,
                   const Eigen::MatrixXd R_Dw,
                   const Eigen::MatrixXd R_Xrand,
                   const Eigen::MatrixXd R_Xfix,
                   const Eigen::MatrixXd R_Zf,
                   const Eigen::VectorXd R_Xoffset,
                   const Eigen::VectorXd R_r,
                   const Eigen::MatrixXd R_K,
                   const Eigen::VectorXd R_a,
                   Eigen::VectorXd R_alpha_f,
                   Eigen::VectorXd R_phi,
                   double R_log_theta,
                   double R_log_smoothing_f,
                   double R_log_smoothing_w,
                   Eigen::VectorXd R_betaR,
                   Eigen::VectorXd R_betaF,
                   Eigen::VectorXd R_logsmoothing) {
    // convert
    Vec y = R_y.cast<Scalar>();
    Mat B_inner = R_B_inner.cast<Scalar>();
    Vec knots_f = R_knots_f.cast<Scalar>();
    Mat Sw_large = R_Sw_large.cast<Scalar>();
    Mat Sf = R_Sf.cast<Scalar>();
    Mat Dw = R_Dw.cast<Scalar>();
    Mat Xrand = R_Xrand.cast<Scalar>();
    Mat Xfix = R_Xfix.cast<Scalar>();
    Mat Zf = R_Zf.cast<Scalar>();
    Vec Xoffset = R_Xoffset.cast<Scalar>();
    Vec r = R_r.cast<Scalar>();
    Mat K = R_K.cast<Scalar>();
    Vec a = R_a.cast<Scalar>();
    Vec alpha_f = R_alpha_f.cast<Scalar>();
    Vec phi = R_phi.cast<Scalar>();
    Scalar log_theta = R_log_theta;
    Scalar log_smoothing_f = R_log_smoothing_f;
    Scalar log_smoothing_w = R_log_smoothing_w;
    Vec betaR = R_betaR.cast<Scalar>();
    Vec betaF = R_betaF.cast<Scalar>();
    Vec logsmoothing = R_logsmoothing.cast<Scalar>();


    Modelcppad* modelcppadobj_ptr = new Modelcppad(y, B_inner, knots_f, Sw_large, Sf, Dw,
                                    Xrand, Xfix, Zf, Xoffset, r, K, a,
                                    alpha_f, phi, log_theta, log_smoothing_f, log_smoothing_w,
                                    betaR, betaF, logsmoothing);

    Rcpp::XPtr<Modelcppad> ptrcppad(modelcppadobj_ptr);

    Model* modelobj_ptr = new Model(R_y, R_B_inner, R_knots_f, R_Sw_large, R_Sf, R_Dw,
                                R_Xrand, R_Xfix, R_Zf, R_Xoffset, R_r, R_K, R_a,
                                R_alpha_f, R_phi, R_log_theta, R_log_smoothing_f, R_log_smoothing_w,
                                R_betaR, R_betaF, R_logsmoothing);
    Rcpp::XPtr<Model> ptr(modelobj_ptr);


    return List::create(Named("address.eigen") = ptr,
                        Named("address.cppad") = ptrcppad);
}






// [[Rcpp::export]]
List aceDLNMopt(SEXP ptr,
                SEXP ptrcppad,
                Eigen::VectorXd R_alpha_f,
                Eigen::VectorXd R_phi,
                double R_log_theta,
                double R_log_smoothing_f,
                double R_log_smoothing_w,
                Eigen::VectorXd R_betaR,
                Eigen::VectorXd R_betaF,
                Eigen::VectorXd R_logsmoothing,
                bool verbose) {


    Rcpp::XPtr<Model> modelobj_ptr(ptr);
    Rcpp::XPtr<Modelcppad> modelcppadobj_ptr(ptrcppad);

    Model& modelobj = *modelobj_ptr;
    Modelcppad& modelcppadobj = *modelcppadobj_ptr;

    modelobj.setAlphaF(R_alpha_f);
    modelobj.setPhi(R_phi);
    modelobj.setBetaR(R_betaR);
    modelobj.setBetaF(R_betaF);
    modelobj.setLogTheta(R_log_theta);
    modelobj.setLogSmoothingF(R_log_smoothing_f);
    modelobj.setLogSmoothingW(R_log_smoothing_w);
    modelobj.setLogsmoothing(R_logsmoothing);



    // Vec alpha_f = R_alpha_f.cast<Scalar>();
    // Vec phi = R_phi.cast<Scalar>();
    // Scalar log_theta = R_log_theta;
    // Scalar log_smoothing_f = R_log_smoothing_f;
    // Scalar log_smoothing_w = R_log_smoothing_w;
    // Vec betaR = R_betaR.cast<Scalar>();
    // Vec betaF = R_betaF.cast<Scalar>();
    // Vec logsmoothing = R_logsmoothing.cast<Scalar>();

    // modelobj.setAlphaF(alpha_f);
    // modelobj.setPhi(phi);
    // modelobj.setBetaR(betaR);
    // modelobj.setBetaF(betaF);
    // modelobj.setLogTheta(log_theta);
    // modelobj.setLogSmoothingF(log_smoothing_f);
    // modelobj.setLogSmoothingW(log_smoothing_w);
    // modelobj.setLogsmoothing(logsmoothing);

    // Inner opt
    Inner(modelobj, verbose);
    // get gr of LAML
    LAMLResult LAMLresult;
    LAMLresult = LAML(modelobj, modelcppadobj);
    // if(hasNaN(LAMLresult.gradient)) {
    //   // rebuild the tape
    //   modelcppadobj.ifhastape = false;
    //   LAMLresult = LAML(modelobj, modelcppadobj);
    // }
    return List::create(Named("LAML.fn") = LAMLresult.fn,
                        Named("LAML.gradient") = LAMLresult.gradient,
                        Named("alpha_f.mod") = modelobj.alpha_f,
                        Named("phi.mod") = modelobj.phi,
                        Named("betaR.mod") = modelobj.betaR,
                        Named("betaF.mod") = modelobj.betaF,
                        Named("address") = modelobj_ptr
                        );
}




// [[Rcpp::export]]
List aceDLNMCI(const Eigen::VectorXd R_y,
                  const Eigen::MatrixXd R_B_inner,
                  const Eigen::VectorXd R_knots_f,
                  const Eigen::MatrixXd R_Sw_large,
                  const Eigen::MatrixXd R_Sf,
                  const Eigen::MatrixXd R_Dw,
                  const Eigen::MatrixXd R_Xrand,
                  const Eigen::MatrixXd R_Xfix,
                  const Eigen::MatrixXd R_Zf,
                  const Eigen::VectorXd R_Xoffset,
                  const Eigen::VectorXd R_r,
                  const Eigen::MatrixXd R_K,
                  const Eigen::VectorXd R_a,
                  Eigen::VectorXd R_alpha_f,
                  Eigen::VectorXd R_phi,
                  double R_log_theta,
                  double R_log_smoothing_f,
                  double R_log_smoothing_w,
                  Eigen::VectorXd R_betaR,
                  Eigen::VectorXd R_betaF,
                  Eigen::VectorXd R_logsmoothing,
                  const int Rci,
                  const int rseed,
                  bool ifeta,
                  bool delta,
                  bool verbose) {


  // construct model
  Model modelobj(R_y, R_B_inner, R_knots_f, R_Sw_large, R_Sf, R_Dw,
                  R_Xrand, R_Xfix, R_Zf, R_Xoffset, R_r, R_K, R_a,
                  R_alpha_f, R_phi, R_log_theta, R_log_smoothing_f, R_log_smoothing_w,
                  R_betaR, R_betaF, R_logsmoothing);


  int kw = modelobj.kw;
  int kwopt = modelobj.kwopt;
  int kE = modelobj.kE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;

  int paraSize = kE+kwopt+kbetaR+kbetaF;
  int paraSizefull;

  // hessian
  Eigen::MatrixXd R_he;
  Eigen::VectorXd R_alpha_w(kw);

  // Vectors for sampling
  Eigen::VectorXd R_phi_sample(kwopt);
  Eigen::VectorXd R_phiKa_sample(kw-1);
  Eigen::VectorXd R_alpha_w_sample(kw);
  Eigen::VectorXd R_alpha_f_sample(kE);
  Eigen::VectorXd R_betaR_sample(kbetaR);
  Eigen::VectorXd R_betaF_sample(kbetaF);

  // Matrices to save results
  Eigen::MatrixXd phi_sample_mat(Rci, kwopt);
  Eigen::MatrixXd alpha_w_sample_mat(Rci, kw);
  Eigen::MatrixXd alpha_f_sample_mat(Rci, kE);
  Eigen::MatrixXd betaR_sample_mat(Rci, kbetaR);
  Eigen::MatrixXd betaF_sample_mat(Rci, kbetaF);

  int n = R_y.size();
  // components for eta
  Eigen::VectorXd R_E;
  Eigen::VectorXd eta_sample;
  Eigen::MatrixXd eta_sample_mat;
  if(ifeta) {
    eta_sample.resize(n);
    eta_sample_mat.resize(Rci, n);
  }

  // Mode of phi
  Eigen::VectorXd R_phi_mod = R_phi;
  // Generate phi
  double R_alpha_w_C_denominator;
  Eigen::VectorXd R_phi_long(kw); // phi_long = c(1, phi)
  R_phi_long(0) = 1.0;
  Eigen::VectorXd R_phiKa = R_K * R_phi + R_a;

  // d alpha_f / d phi
  Eigen::MatrixXd R_deriv_g(kw,kwopt);
  Eigen::MatrixXd R_deriv_g_large(kw,kw);
  Eigen::MatrixXd R_DiagMat_tmp(kw,kw);
  R_DiagMat_tmp.setZero();
  double R_diag_tmp;

  Eigen::MatrixXd R_deriv(paraSize+1, paraSize);
  R_deriv.setZero();

  paraSizefull = paraSize;

  // Joint
  R_he = modelobj.he_s_u_mat;
  Eigen::VectorXd R_u_mod(paraSizefull);
  // Hessian
  // cholesky of inverse Hessian
  Eigen::MatrixXd R_he_u_L(paraSize, paraSize);
  Eigen::MatrixXd R_he_u_L_inv(paraSizefull, paraSize);
  Eigen::VectorXd zjoint(paraSize);
  Eigen::VectorXd samplejoint(paraSizefull);


  R_u_mod << R_alpha_f, R_phi, R_betaR, R_betaF;

  // cholesky of inverse Hessian
  R_he_u_L = R_he.llt().matrixL();
  if(delta) {
    R_he_u_L_inv = R_deriv * (invertL(R_he_u_L)).transpose();
  } else {
    R_he_u_L_inv = (invertL(R_he_u_L)).transpose();
  }

  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(rseed);
  std::normal_distribution<> dist(0, 1);


  for(int i = 0; i < Rci; i++)
  {
    // Jointly sample
    for (int j = 0; j < paraSize; j++) {
      zjoint(j) = dist(gen);
    }
    samplejoint = R_u_mod + R_he_u_L_inv * zjoint;
    // get alpha_f
    R_alpha_f_sample = samplejoint.segment(0, kE);

    if(delta) {
      R_alpha_w_sample = samplejoint.segment(kE, kw);
      R_betaR_sample = samplejoint.segment(kE+kw, kbetaR);
      R_betaF_sample = samplejoint.segment(kE+kw+kbetaR, kbetaF);
    } else {
      // get phi
      R_phi_sample = samplejoint.segment(kE, kwopt);
      R_phiKa_sample = R_K * R_phi_sample + R_a;
      // get alpha_w
      for (int j = 0; j < (kw - 1); j++) {
        R_phi_long(j + 1) = R_phiKa_sample(j);
      }
      R_alpha_w_C_denominator = sqrt(R_phi_long.dot(R_Dw * R_phi_long));
      R_alpha_w_sample = R_phi_long / R_alpha_w_C_denominator;
      // get betaR
      R_betaR_sample = samplejoint.segment(kE+kwopt, kbetaR);
      // get betaF
      R_betaF_sample = samplejoint.segment(kE+kwopt+kbetaR, kbetaF);

      // save
      phi_sample_mat.row(i) = R_phi_sample.transpose();
    }

    if(ifeta) {
      R_E = R_B_inner * R_alpha_w_sample;
      for (int ii = 0; ii < n; ii++) {
        eta_sample(ii) = BsplinevecCon(R_E(ii), R_knots_f, 4, R_Zf).dot(R_alpha_f_sample) + R_Xfix.row(ii).dot(R_betaF_sample) + R_Xrand.row(ii).dot(R_betaR_sample) + R_Xoffset(ii);
      }
      eta_sample_mat.row(i) = eta_sample.transpose();
    }

    // save
    alpha_w_sample_mat.row(i) = R_alpha_w_sample.transpose();
    alpha_f_sample_mat.row(i) = R_alpha_f_sample.transpose();
    betaR_sample_mat.row(i) = R_betaR_sample.transpose();
    betaF_sample_mat.row(i) = R_betaF_sample.transpose();
  }
  if(delta) {
    return List::create(Named("alpha_w_sample") = alpha_w_sample_mat,
                        Named("alpha_f_sample") = alpha_f_sample_mat,
                        Named("betaR_sample") = betaR_sample_mat,
                        Named("betaF_sample") = betaF_sample_mat,
                        Named("eta_sample_mat") = eta_sample_mat,
                        Named("Hessian_inner") = R_he);
  } else {
    return List::create(Named("phi_sample") = phi_sample_mat,
                      Named("alpha_w_sample") = alpha_w_sample_mat,
                      Named("alpha_f_sample") = alpha_f_sample_mat,
                      Named("betaR_sample") = betaR_sample_mat,
                      Named("betaF_sample") = betaF_sample_mat,
                      Named("eta_sample_mat") = eta_sample_mat,
                      Named("Hessian_inner") = R_he);
  }

}



// [[Rcpp::export]]
List ConditionalAIC(const Eigen::VectorXd R_y,
                      const Eigen::MatrixXd R_B_inner,
                      const Eigen::VectorXd R_knots_f,
                      const Eigen::MatrixXd R_Sw_large,
                      const Eigen::MatrixXd R_Sf,
                      const Eigen::MatrixXd R_Dw,
                      const Eigen::MatrixXd R_Xrand,
                      const Eigen::MatrixXd R_Xfix,
                      const Eigen::MatrixXd R_Zf,
                      const Eigen::VectorXd R_Xoffset,
                      const Eigen::VectorXd R_r,
                      const Eigen::MatrixXd R_K,
                      const Eigen::VectorXd R_a,
                      Eigen::VectorXd R_alpha_f,
                      Eigen::VectorXd R_phi,
                      double R_log_theta,
                      double R_log_smoothing_f,
                      double R_log_smoothing_w,
                      Eigen::VectorXd R_betaR,
                      Eigen::VectorXd R_betaF,
                      Eigen::VectorXd R_logsmoothing) {

  // construct model
  Model modelobj(R_y, R_B_inner, R_knots_f, R_Sw_large, R_Sf, R_Dw,
                  R_Xrand, R_Xfix, R_Zf, R_Xoffset, R_r,
                  R_K, R_a,
                  R_alpha_f, R_phi, R_log_theta, R_log_smoothing_f, R_log_smoothing_w,
                  R_betaR, R_betaF, R_logsmoothing);
  modelobj.prepare_AIC();
  // hessian
  Eigen::MatrixXd R_he;
  R_he = modelobj.he_s_u_mat;
  // I
  Eigen::MatrixXd R_I;
  R_I = modelobj.I_mat;
  //
  Eigen::MatrixXd mat_AIC = R_he.ldlt().solve(R_I);

  double l = modelobj.NegLogL_l;
  double edf = (2 * mat_AIC - mat_AIC * mat_AIC).trace();
  double AIC = 2.0*l + 2.0*edf;

  return List::create(Named("AIC") = AIC,
                      Named("l") = -1.0*l,
                      Named("edf") = edf);
}



// TODO: revise the code ... THE CURRENT CODE BELOW ARE WRONG!
// ************** PART 3: g(mu) = DL term + linear term. WITHOUT other smooth terms **********************

class Model_nosmooth {
  // The DLNM model

private:
  // DATA
  const Vec& y; // Response
  const Mat& Sw; // penalty matrix for w(l)
  const Mat& Sf; // penalty matrix for f(E)
  const Mat& B_inner;
  const Vec& knots_f; // knots for f(E) B-spline
  const Mat& Dw;  // \int w(l)^2 dl = 1

  const Mat& Xfix; // fixed effects
  const Mat& Zf; // for point contraint in f(E)

  const Vec& Xoffset; // offset

public:
  int n;
  int kw;
  int kE;
  int kbetaF;

  // PARAMETERS
  Vec alpha_f;
  Vec phi;
  Scalar log_theta;
  Scalar log_smoothing_f;
  Scalar log_smoothing_w;

  Vec betaF; // parameters for fixed effects

  // Components generated
  Scalar theta;
  Scalar smoothing_f;
  Scalar smoothing_w;
  Vec phi_long;
  Scalar alpha_w_C_denominator;
  Vec alpha_w_C;
  Vec alpha_w_C_pen;
  Mat Bf_matrix;
  Vec E;
  Vec eta;
  Vec eta_remaining; // remaining terms = Xfix * betaF
  Vec mu; // log(mu) = eta + eta_remaining + offset
  Scalar NegLogL; // NegativeLogLikelihood value


  // Components for derivatives

  Mat dlogmu_df_mat;
  Mat dlogmu_dbetaF_mat;
  Mat dlogmu_dw_mat;
  Vec dlogdensity_dmu_vec;
  Mat dmu_df_mat;
  Mat dmu_dbetaF_mat;
  Mat dmu_dw_mat;
  Mat dw_dphi_mat;
  Vec gr_alpha_w_vec;
  Vec d2logdensity_dmudmu_vec;
  std::vector<Mat> d2mu_dfdf_list;
  std::vector<Mat> d2logmu_dwdw_list;
  std::vector<Mat> d2mu_dbetaFdbetaF_list;
  std::vector<Mat> d2mu_dwdw_list;
  std::vector<Mat> d2w_dphidphi_list;
  std::vector<Mat> d2logmu_dfdw_list;
  std::vector<Mat> d2mu_dfdw_list;
  Mat he_alpha_w_mat;
  Mat he_alpha_f_alpha_w_mat;
  Scalar dlogdensity_dtheta_scalar;
  // Scalar d2logdensity_dthetadtheta_scalar;
  Vec d2logdensity_dmudtheta_vec;

  // gradient and hessian for updating alpha_f and betaF
  Vec gr_inner_vec;
  Mat he_inner_mat;

  // full gradient
  Vec gr_alpha_f_vec;
  Vec gr_betaF_vec;
  Vec gr_phi_vec;
  Scalar gr_log_smoothing_f_scalar;
  Scalar gr_log_smoothing_w_scalar;
  Scalar gr_log_theta_scalar;

  Vec gr_s_u_vec;
  Vec gr_s_par_vec;

  // full hessian
  Mat he_alpha_f_mat;
  Mat he_betaF_mat;
  Mat he_phi_mat;
  Mat he_alpha_f_phi_mat;
  Mat he_alpha_f_betaF_mat;
  Mat he_phi_betaF_mat;
  // Scalar he_log_smoothing_f_scalar;
  // Scalar he_log_smoothing_w_scalar;
  // Scalar he_log_theta_scalar;
  Vec he_alpha_f_log_smoothing_f_vec;
  Vec he_phi_log_smoothing_w_vec;
  Vec he_alpha_f_log_theta_vec;
  Vec he_phi_log_theta_vec;
  Vec he_betaF_log_theta_vec;

  Mat he_s_u_mat;
  Mat he_s_par_u_mat;

  // To compute AIC
  Scalar NegLogL_l; // NegativeLogLikelihood without penalty
  // matrix for I (hessian of log likelihood without penalty)
  Mat I_alpha_f_mat;
  Mat I_phi_mat;
  Mat I_alpha_w_mat;
  Mat I_mat;

  // results for profile likelihood
  Eigen::VectorXd PL_gradient;
  Eigen::MatrixXd PL_hessian;
  int converge; // 0: converge. 99: not converge

  // Constructor
  Model_nosmooth(const Vec& y_,
        const Mat& B_inner_,
        const Vec& knots_f_,
        const Mat& Sw_,
        const Mat& Sf_,
        const Mat& Dw_,
        const Mat& Xfix_,
        const Mat& Zf_,
        const Vec& Xoffset_,
        Vec& alpha_f_,
        Vec& phi_,
        Scalar log_theta_,
        Scalar log_smoothing_f_,
        Scalar log_smoothing_w_,
        Vec& betaF_) :
    y(y_), B_inner(B_inner_), knots_f(knots_f_), Sw(Sw_), Sf(Sf_), Dw(Dw_), Xfix(Xfix_), Zf(Zf_), Xoffset(Xoffset_),
    alpha_f(alpha_f_), phi(phi_), log_theta(log_theta_), log_smoothing_f(log_smoothing_f_), log_smoothing_w(log_smoothing_w_), betaF(betaF_) {

      n = y.size(); // sample size
      kw = phi.size() + 1;
      kE = alpha_f.size(); // It is different from kE.mgcv (alpha_f.size()) in mgcv.
      kbetaF = betaF.size();

      theta = exp(log_theta);
      smoothing_f = exp(log_smoothing_f);
      smoothing_w = exp(log_smoothing_w);

      phi_long.resize(kw); // phi_long = c(1, phi)
      phi_long(0) = 1.0;
      for (int j = 0; j < (kw - 1); j++) {
        phi_long(j + 1) = phi(j);
      }
      alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));
      alpha_w_C = phi_long / alpha_w_C_denominator;
      alpha_w_C_pen = phi / alpha_w_C_denominator;

      E = B_inner * alpha_w_C;

      Bf_matrix.resize(n, kE);
      eta.resize(n);
      eta_remaining.resize(n);
      mu.resize(n);
      Vec Bf;
      for (int i = 0; i < n; i++) {
        Bf = BsplinevecCon(E(i), knots_f, 4, Zf);
        Bf_matrix.row(i) = Bf;
        eta(i) = Bf.dot(alpha_f);
        eta_remaining(i) = Xfix.row(i).dot(betaF);
        mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
      }

      // Initialize the derivative components and NegativeLogLikelihood
      dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
      d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi

      gr_s_u_vec.resize(kw+kE-1+kbetaF);
      he_s_u_mat.resize(kw+kE-1+kbetaF, kw+kE-1+kbetaF);
      gr_s_par_vec.resize(3);
      he_s_par_u_mat.resize(3, kw+kE-1+kbetaF);

      gr_inner_vec.resize(kE+kbetaF);
      he_inner_mat.resize(kE+kbetaF, kE+kbetaF);

      derivative_coef();
      derivative_he();
      derivative_full();
      NegativeLogLikelihood();

      // Initialize PL_nosmooth
      PL_gradient.resize(kw-1);
      PL_hessian.resize(kw-1, kw-1);
    }

  // Functions to set parameters
  void setAlphaF(const Vec alpha_f_) {
    alpha_f = alpha_f_;

    for (int i = 0; i < n; i++) {
      eta(i) = Bf_matrix.row(i).dot(alpha_f);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }

  void setPhi(const Vec phi_) {
    phi = phi_;

    // re-generate
    for (int j = 0; j < (kw - 1); j++) {
      phi_long(j + 1) = phi(j);
    }
    alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));
    alpha_w_C = phi_long / alpha_w_C_denominator;
    alpha_w_C_pen = phi / alpha_w_C_denominator;

    E = B_inner * alpha_w_C;
    Vec Bf;
    for (int i = 0; i < n; i++) {
      Bf = BsplinevecCon(E(i), knots_f, 4, Zf);
      Bf_matrix.row(i) = Bf;
      eta(i) = Bf.dot(alpha_f);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }

    dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
    d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi
  }

  void setBetaF(const Vec betaF_) {
    betaF = betaF_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }

  void setLogTheta(const Scalar log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const Scalar log_smoothing_f_) {
    log_smoothing_f = log_smoothing_f_;
    smoothing_f = exp(log_smoothing_f);
  }

  void setLogSmoothingW(const Scalar log_smoothing_w_) {
    log_smoothing_w = log_smoothing_w_;
    smoothing_w = exp(log_smoothing_w);
  }

  // get private members
  Mat getB_inner () {
    return B_inner;
  }

  Mat getDw () {
    return Dw;
  }

  Vec getknots_f () {
    return knots_f;
  }

  // Function to update derivatives.
  // RUN the function derivative_coef(), derivative_he() and derivative_full() after update parameters.
  // update derivatives related to spline coefficients alpha_f and phi
  void derivative_coef() {
    dlogmu_dw_mat = dlogmu_dw();
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();
    dmu_dbetaF_mat = dmu_dbetaF();
    dmu_dw_mat = dmu_dw();
    gr_alpha_w_vec = gr_alpha_w();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
    d2mu_dfdf_list = d2mu_dfdf();
    d2mu_dbetaFdbetaF_list = d2mu_dbetaFdbetaF();
    d2logmu_dwdw_list = d2logmu_dwdw();
    d2mu_dwdw_list = d2mu_dwdw();
    he_alpha_w_mat = he_alpha_w();
    d2logmu_dfdw_list = d2logmu_dfdw();
    d2mu_dfdw_list = d2mu_dfdw();
    he_alpha_f_alpha_w_mat = he_alpha_f_alpha_w();
    dlogdensity_dtheta_scalar = dlogdensity_dtheta();
    // d2logdensity_dthetadtheta_scalar = d2logdensity_dthetadtheta();
    d2logdensity_dmudtheta_vec = d2logdensity_dmudtheta();

    // obtain gradient
    gr_alpha_f_vec = gr_alpha_f();
    gr_betaF_vec = gr_betaF();
    gr_phi_vec = gr_phi();
    // obtain hessian
    he_alpha_f_mat = he_alpha_f();
    he_betaF_mat = he_betaF();
    he_phi_mat = he_phi();
    he_alpha_f_phi_mat = he_alpha_f_phi();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_phi_betaF_mat = he_phi_betaF();
  }

  // update full gradient and hessian of alpha_f, phi and betaF
  void derivative_he () {
    gr_s_u_vec << gr_alpha_f_vec, gr_phi_vec, gr_betaF_vec;

    he_s_u_mat.setZero();
    // he_s_mat = [he_alpha_f_mat, he_alpha_f_phi_mat, he_alpha_f_log_theta_vec, he_alpha_f_log_smoothing_f_vec, 0;
    //             he_alpha_f_phi_mat.transpose(), he_phi_mat, he_phi_log_theta_vec, 0, he_phi_log_smoothing_w_vec;
    //             he_alpha_f_log_theta_vec.transpose(), he_phi_log_theta_vec.transpose(), he_log_theta_scalar, 0, 0;
    //             he_alpha_f_log_smoothing_f_vec.transpose(), 0, 0, he_log_smoothing_f_scalar, 0;
    //             0, he_phi_log_smoothing_w_vec.transpose(), 0, 0, he_log_smoothing_w_scalar]
    he_s_u_mat.block(0, 0, kE, kE)  = he_alpha_f_mat;
    he_s_u_mat.block(0, kE, kE, kw-1) = he_alpha_f_phi_mat;
    he_s_u_mat.block(kE, 0, kw-1, kE) = he_alpha_f_phi_mat.transpose();
    he_s_u_mat.block(kE, kE, kw-1, kw-1) = he_phi_mat;

    he_s_u_mat.block(kE+kw-1, kE+kw-1, kbetaF, kbetaF) = he_betaF_mat;
    he_s_u_mat.block(0,kE+kw-1,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_s_u_mat.block(kE+kw-1,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();

    he_s_u_mat.block(kE,kE+kw-1,kw-1,kbetaF) = he_phi_betaF_mat;
    he_s_u_mat.block(kE+kw-1,kE,kbetaF,kw-1) = he_phi_betaF_mat.transpose();

    // make it symmetric. Comment out...
    // he_s_u_mat = (he_s_u_mat + he_s_u_mat.transpose())/2.0;
  }

  // update derivatives related to overdispersion and smoothing parameters
  // Full derivative for LAML
  void derivative_full () {
    // obtain full gradient

    gr_log_smoothing_f_scalar = gr_log_smoothing_f();
    gr_log_smoothing_w_scalar = gr_log_smoothing_w();
    gr_log_theta_scalar = gr_log_theta();

    // u represents spline coefficient alpha_f and phi and betaF
    // par represents overdispersion and smoothing parameters

    gr_s_par_vec << gr_log_theta_scalar, gr_log_smoothing_f_scalar, gr_log_smoothing_w_scalar;

    // obtain full hessian
    // he_log_smoothing_f_scalar = he_log_smoothing_f();
    // he_log_smoothing_w_scalar = he_log_smoothing_w();
    // he_log_theta_scalar = he_log_theta();
    he_alpha_f_log_smoothing_f_vec = he_alpha_f_log_smoothing_f();
    he_phi_log_smoothing_w_vec = he_phi_log_smoothing_w();
    he_alpha_f_log_theta_vec = he_alpha_f_log_theta();
    he_phi_log_theta_vec = he_phi_log_theta();
    he_betaF_log_theta_vec = he_betaF_log_theta();

    he_s_par_u_mat.setZero();


    he_s_par_u_mat.row(0) << he_alpha_f_log_theta_vec.transpose(), he_phi_log_theta_vec.transpose(), he_betaF_log_theta_vec.transpose();
    he_s_par_u_mat.block(1, 0, 1, kE) = he_alpha_f_log_smoothing_f_vec.transpose();
    he_s_par_u_mat.block(2, kE, 1, kw-1) = he_phi_log_smoothing_w_vec.transpose();
  }



  // update variables related to alpha_f.
  // Used only in updating alpha_f.
  void derivative_f () {
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();
    dmu_dbetaF_mat = dmu_dbetaF();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
    d2mu_dfdf_list = d2mu_dfdf();
    d2mu_dbetaFdbetaF_list = d2mu_dbetaFdbetaF();

    gr_alpha_f_vec = gr_alpha_f();
    gr_betaF_vec = gr_betaF();

    he_alpha_f_mat = he_alpha_f();
    he_betaF_mat = he_betaF();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();

    gr_inner_vec << gr_alpha_f_vec, gr_betaF_vec;

    he_inner_mat.setZero();
    he_inner_mat.block(0,0,kE,kE) = he_alpha_f_mat;
    he_inner_mat.block(kE,kE,kbetaF,kbetaF) = he_betaF_mat;

    he_inner_mat.block(0,kE,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_inner_mat.block(kE,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
  }


  // functions for NegativeLogLikelihood
  void NegativeLogLikelihood() {

    Scalar loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lanczos_lgamma(y(i) + theta) - lanczos_lgamma(theta) - lanczos_lgamma(y(i) + 1) -
                                   theta * log(1 + mu(i)/theta) +
                                   y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }
    // part 1: DLNM
    // Smooth Penalty
    loglik += -0.5 * smoothing_w * alpha_w_C_pen.dot(Sw * alpha_w_C_pen) - 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f);
    // Scale
    loglik += (kw-1-1) / 2.0 * log_smoothing_w + (kE-1) / 2.0 * log_smoothing_f;

    // part 2: Remaining smooth terms
    // EMPTY

    NegLogL = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  // functions for NegativeLogLikelihood WITHOUT penalty for AIC
  void NegativeLogLikelihood_l() {

    Scalar loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lanczos_lgamma(y(i) + theta) - lanczos_lgamma(theta) - lanczos_lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }

    NegLogL_l = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  void prepare_AIC () {
    NegativeLogLikelihood_l();
     // hessian of log likelihood without penalty
    I_alpha_f_mat = I_alpha_f();
    I_alpha_w_mat = I_alpha_w();
    I_phi_mat = I_phi();
    I_mat = he_s_u_mat;
    I_mat.block(0, 0, kE, kE)  = I_alpha_f_mat;
    I_mat.block(kE, kE, kw-1, kw-1) = I_phi_mat;
  }

  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Vec dlogdensity_dmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu^2
  Vec d2logdensity_dmudmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }
  // d log(exponential family density) / d theta
  Scalar dlogdensity_dtheta () {
    Scalar out = 0.0;
    // std::cout << "x" << 3.5 << std::endl;
    // std::cout << "lgamma1st" << CppAD::Value(lgamma1st(3.5)) << std::endl;

    // TO DO: optimize it. Use property of gamma function...
    for (int i = 0; i < n; i++) {
      out += log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    }
    return out;
  }
  // d^2 log(exponential family density) / d theta^2
  // Scalar d2logdensity_dthetadtheta () {
  //   Scalar out = 0.0;

  //   for (int i = 0; i < n; i++) {
  //     out += 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
  //   }
  //   return out;
  // }
  Vec d2logdensity_dmudtheta () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = (y(i) - mu(i)) / pow(theta+mu(i), 2);
    }
    return out;
  }



  // 2. mean model
  // d log(mu) / d alpha_f
  Mat dlogmu_df () {
    return Bf_matrix;
  }
  // d mu / d alpha_f
  Mat dmu_df () {
    Mat out(n, kE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_df_mat.row(i) * mu(i);
    }
    return out;
  }
  // d log(mu) / d alpha_w
  Mat dlogmu_dw () {
    Mat out(n, kw);
    Vec Bf1st;
    for (int i = 0; i < n; i++) {
      Bf1st = BsplinevecCon1st(E(i), knots_f, 4, Zf);
      out.row(i) = B_inner.row(i) * (Bf1st.dot(alpha_f));
    }
    return out;
  }
  // d mu / d alpha_w
  Mat dmu_dw () {
    Mat out(n, kw);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dw_mat.row(i) * mu(i);
    }
    return out;
  }

  // d log(mu) / d betaF
  Mat dlogmu_dbetaF () {
    return Xfix;
  }
  // d mu / d betaF
  Mat dmu_dbetaF () {
    Mat out(n, kbetaF);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaF_mat.row(i) * mu(i);
    }
    return out;
  }


  // d^2 mu / d alpha_f^2
  std::vector<Mat> d2mu_dfdf () {
    std::vector<Mat> out;
    for (int i = 0; i < n; i++) {
      out.push_back(mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }
    return out;
  }
  // d^2 mu / d betaF^2
  std::vector<Mat> d2mu_dbetaFdbetaF () {
    std::vector<Mat> out;
    for (int i = 0; i < n; i++) {
      out.push_back(mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));
    }
    return out;
  }
  // d^2 log(mu) / d alpha_w^2
  std::vector<Mat> d2logmu_dwdw () {
    std::vector<Mat> out;
    for (int i = 0; i < n; i++) {
      Vec Bf2nd = BsplinevecCon2nd(E(i), knots_f, 4, Zf);
      out.push_back((Bf2nd.dot(alpha_f)) * B_inner.row(i).transpose() * B_inner.row(i));
    }
    return out;
  }
  // d^2 mu / d alpha_w^2
  std::vector<Mat> d2mu_dwdw () {
    std::vector<Mat> out;
    for (int i = 0; i < n; i++) {
      out.push_back(mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i) * d2logmu_dwdw_list.at(i));
    }
    return out;
  }
  // d^2 log(mu) / d alpha_f d alpha_w
  std::vector<Mat> d2logmu_dfdw () {
    std::vector<Mat> out;
    for (int i = 0; i < n; i++) {
      out.push_back(BsplinevecCon1st(E(i), knots_f, 4, Zf) * B_inner.row(i));
    }
    return out;
  }
  // d^2 mu / d alpha_f d alpha_w
  std::vector<Mat> d2mu_dfdw () {
    std::vector<Mat> out;
    for (int i = 0; i < n; i++) {
      out.push_back(mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i)*d2logmu_dfdw_list.at(i));
    }
    return out;
  }



  // 3. Re-parameterization
  // d alpha_w / d phi
  Mat dw_dphi () {
    // deriv_g <- diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw) - phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    // deriv_g <- deriv_g[1:(kw), 2:kw]

    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));

    // diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw)
    Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D(kw);
    D.diagonal().setConstant(1.0/alpha_w_C_denominator);
    Mat Ddense = D.toDenseMatrix();

    // phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    Mat deriv_g2 = (phi_long * (Dw * phi_long).transpose()) * pow(alpha_w_C_denominator, -3);

    Mat deriv_g = Ddense - deriv_g2;
    // Remove the first column
    Mat out = deriv_g.block(0, 1, kw, kw - 1);
    return out;
  }


  std::vector<Mat> d2w_dphidphi () {
    std::vector<Mat> out;
    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long)) = tmp^(1/2);
    Scalar tmp1 = pow(alpha_w_C_denominator, -3); // pow(tmp, -1.5)
    Scalar tmp2 = pow(alpha_w_C_denominator, -5); // pow(tmp, -2.5)
    Vec Dwphi = Dw * phi_long;
    Mat outlarge(kw, kw);
    for (int s = 0; s < kw; s++) {
      if (s == 0) {
        outlarge = -1.0 * tmp1 * Dw + 3.0 * tmp2 * Dwphi * Dwphi.transpose();
      } else {
        Mat m1(kw, kw);
        m1.setZero();
        m1.row(s) = Dwphi.transpose()*tmp1;
        m1.col(s) = m1.col(s) + Dwphi*tmp1;
        Mat m2 = -1.0 * tmp1* Dw + 3.0 * tmp2 * Dwphi * Dwphi.transpose();
        outlarge = -1.0 * m1 + m2; // or m1 - m2 or m2 - m1
      }
      out.push_back(outlarge.block(1, 1, kw-1, kw-1));
    }
    return out;
  }

  // *** GRADIENT ***
  Vec gr_alpha_f () {
    Vec out = - dmu_df_mat.transpose() * dlogdensity_dmu_vec + smoothing_f * Sf * alpha_f;
    return out;
  }
  Vec gr_betaF () {
    Vec out = - dmu_dbetaF_mat.transpose() * dlogdensity_dmu_vec;
    return out;
  }
  Vec gr_alpha_w () {
    Vec gr_pen_w = smoothing_w * Sw * alpha_w_C_pen;
    Vec gr_pen_w_long(kw);
    gr_pen_w_long(0) = 0.0;
    for (int j = 0; j < (kw - 1); j++) {
      gr_pen_w_long(j + 1) = gr_pen_w(j);
    }

    Vec out = - dmu_dw_mat.transpose() * dlogdensity_dmu_vec + gr_pen_w_long;
    return out;
  }
  Vec gr_phi () {
    Vec out = dw_dphi_mat.transpose() * gr_alpha_w_vec;
    return out;
  }
  Scalar gr_log_smoothing_f () {
    return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * (kE-1);
  }
  Scalar gr_log_smoothing_w () {
    return 0.5 * smoothing_w * alpha_w_C_pen.dot(Sw * alpha_w_C_pen) - 0.5 * (kw-1-1);
  }
  Scalar gr_log_theta () {
    return -1.0 * theta * dlogdensity_dtheta_scalar;
  }


  // *** Hessian ***
  Mat he_alpha_f () {
    Mat out1(kE, kE);
    Mat out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * d2mu_dfdf_list.at(i);
    }
    return - out1 - out2 + smoothing_f*Sf;
  }
  Mat I_alpha_f () {
    Mat out1(kE, kE);
    Mat out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * d2mu_dfdf_list.at(i);
    }
    return - out1 - out2;
  }
  Mat he_betaF () {
    Mat out1(kbetaF, kbetaF);
    Mat out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * d2mu_dbetaFdbetaF_list.at(i);
    }
    return - out1 - out2;
  }
  Mat he_alpha_w () {
    Mat out1(kw, kw);
    Mat out2(kw, kw);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * d2mu_dwdw_list.at(i);
    }
    Mat Sw_large(kw, kw);
    Sw_large.setZero();
    Sw_large.block(1, 1, kw-1, kw-1) = Sw;
    return - out1 - out2 + smoothing_w*Sw_large;
  }

  Mat I_alpha_w () {
    Mat out1(kw, kw);
    Mat out2(kw, kw);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * d2mu_dwdw_list.at(i);
    }
    return - out1 - out2;
  }

  Mat he_phi () {
    Mat out1 = dw_dphi_mat.transpose() * he_alpha_w_mat * dw_dphi_mat;
    Mat out2(kw-1, kw-1);
    out2.setZero();
    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }

  Mat I_phi () {
    Mat out1 = dw_dphi_mat.transpose() * I_alpha_w_mat * dw_dphi_mat;
    Mat out2(kw-1, kw-1);
    out2.setZero();
    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }

  Mat he_alpha_f_alpha_w () {
    Mat out1(kE, kw);
    Mat out2(kE, kw);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * d2mu_dfdw_list.at(i);
    }
    return - out1 - out2;
  }

  Mat he_alpha_f_phi () {
    Mat out = he_alpha_f_alpha_w_mat * dw_dphi_mat;
    return out;
  }
  Mat he_alpha_f_betaF () {
    Mat out1(kE, kbetaF);
    Mat out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }
  Mat he_phi_betaF () {
    Mat out1(kw-1, kbetaF);
    Mat out2(kw-1, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
    }
    return - out1 - out2;
  }

  // Scalar he_log_smoothing_f () {
  //   return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f);
  // }
  // Scalar he_log_smoothing_w () {
  //   return 0.5 * smoothing_w * alpha_w_C_pen.dot(Sw * alpha_w_C_pen);
  // }
  // Scalar he_log_theta () {
  //   return -1.0*theta*theta * d2logdensity_dthetadtheta_scalar - theta * dlogdensity_dtheta_scalar;
  // }

  Vec he_alpha_f_log_smoothing_f () {
    return smoothing_f * Sf * alpha_f;
  }
  Vec he_phi_log_smoothing_w () {

    Vec he_alpha_w_C_pen_log_smoothing_w = smoothing_w * Sw * alpha_w_C_pen;

    Vec he_alpha_w_log_smoothing_w(kw);
    he_alpha_w_log_smoothing_w(0) = 0.0;
    for (int i = 1; i < kw; i++)
    {
      he_alpha_w_log_smoothing_w(i) = he_alpha_w_C_pen_log_smoothing_w(i-1);
    }

    return dw_dphi_mat.transpose() * he_alpha_w_log_smoothing_w;
  }

  Vec he_alpha_f_log_theta () {
    // he_alpha_f_theta = dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Vec he_phi_log_theta () {
    // he_alpha_w_theta = dmu_dw_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*theta * dw_dphi_mat.transpose() * ( dmu_dw_mat.transpose() * d2logdensity_dmudtheta_vec );
  }
  Vec he_betaF_log_theta () {
    return -1.0*dmu_dbetaF_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }




  // *********** LAML ***********
  Scalar logdetH05() {
    // Scalar out = 0.5 * log(he_s_u_mat.determinant());
    // return out;

    Scalar logdetH05 = 0.0;
    Eigen::PartialPivLU<Mat> lu(he_s_u_mat);
    Mat LU = lu.matrixLU();
    // Scalar c = lu.permutationP().determinant(); // -1 or 1
    Scalar lii;
    for (int i = 0; i < LU.rows(); i++) {
      lii = LU(i,i);
      logdetH05 += log(CppAD::abs(lii));
    }
    return logdetH05/2.0;
  }

};







// optimize alpha_f and betaF for a given phi, log_smoothing_f, log_smoothing_w, and log_theta
// Use newton method with eigenvalue modification
void PL_nosmooth(Model_nosmooth& modelobj, bool verbose){
    int maxitr = 50, itr = 0;
    const double eps = 1e-05;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 50, stephalve = 0;
    int resetitr = 0, maxreset = 1; // if step is nan, reset coefficients as 0. only once.
    int additr = 50; // allow further iterations after resetting coefficients as 0

    Vec alpha_f = modelobj.alpha_f;
    Vec betaF = modelobj.betaF;
    Eigen::VectorXd R_alpha_f = convertToDouble(alpha_f);
    Eigen::VectorXd R_betaF = convertToDouble(betaF);

    int kE = modelobj.kE;
    int kw = modelobj.kw;
    int kbetaF = modelobj.kbetaF;
    int converge = 0;

    int paraSize = kE+kbetaF;
    // Optimize ALPHA_F
    Scalar u;
    Scalar u_tmp;

    Mat H;
    Vec g;


    // update steps
    Eigen::VectorXd step(paraSize);
    step.setZero();



    // Cast results
    Eigen::MatrixXd R_H(paraSize, paraSize);
    Eigen::VectorXd R_g(paraSize);

    R_g.setZero();
    R_H.setZero();

    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = R_H * eigen_in; // Easy version
      eigen_out.noalias() += R_H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, paraSize, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smalles eigenvalue

    // eigen decomposition
    Eigen::VectorXd eigvals(paraSize);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(paraSize);
    invabseigvals.setZero();
    // double eigval;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(R_H,false); // Only values, not vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(R_H,true); // Both values and vectors

    double delta = 1.0; // for New Q-Newton method
    double R_g_norm;

    for (int i = 0; i < paraSize; i++) R_g(i) = 1. + eps;

    while (itr < maxitr) {
      itr++;
      modelobj.derivative_f();

      g = modelobj.gr_inner_vec;
      R_g = convertToDouble(g);

      R_g_norm = R_g.norm();

      if (R_g_norm < eps) break;

      modelobj.NegativeLogLikelihood();
      u = modelobj.NegLogL;

      H = modelobj.he_inner_mat;
      R_H = convertToDouble(H);

      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        eigvec.compute(R_H); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();
        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < paraSize; iii++) eigvals(iii) += delta*R_g_norm;
        }

        // for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / max(abs(eigvals(i)), mineig); // flip signs
        for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * R_g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = R_H.llt().solve(R_g);
        step = R_H.ldlt().solve(R_g);
      }

      // check nan in step
      // Really needed
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          R_alpha_f.setZero(); // reset alpha_f
          R_betaF.setZero();
          alpha_f = R_alpha_f.cast<Scalar>();
          betaF = R_betaF.cast<Scalar>();
          modelobj.setAlphaF(alpha_f);
          modelobj.setBetaF(betaF);
          if(verbose) std::cout << "reset alpha_f and betaF as 0" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          converge = 99;
          break;
        }
      }

      R_alpha_f -= step.segment(0, kE);
      R_betaF -= step.segment(kE, kbetaF);

      alpha_f = R_alpha_f.cast<Scalar>();
      betaF = R_betaF.cast<Scalar>();

      modelobj.setAlphaF(alpha_f);
      modelobj.setBetaF(betaF);

      modelobj.NegativeLogLikelihood();
      u_tmp = modelobj.NegLogL;

      // halving if objective function increase.
      stephalve = 0;
      while ((CppAD::Value(u_tmp) > CppAD::Value(u) + 1e-8) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;

        R_alpha_f += step.segment(0, kE);
        R_betaF += step.segment(kE, kbetaF);
        alpha_f = R_alpha_f.cast<Scalar>();
        betaF = R_betaF.cast<Scalar>();
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaF(betaF);

        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }


      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(CppAD::Value(u_tmp)) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration

        R_alpha_f += step.segment(0, kE);
        R_betaF += step.segment(kE, kbetaF);
        alpha_f = R_alpha_f.cast<Scalar>();
        betaF = R_betaF.cast<Scalar>();
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaF(betaF);

        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }

      stephalve = 0;
      // if (stephalve > 0) std::cout << "Performed " << stephalve << " iterations of step-halving." << std::endl;
      if (std::isnan(CppAD::Value(u_tmp))) {
        // Step-halving didn't work
        // std::cout << "AlphaF: Step-halving failed with nan function value. Returning failure." << std::endl;
        converge = 99;
        break;
      }
    }
    if(itr == maxitr){
      // std::cout << "AlphaF: Newton method for updating alpha fails" << std::endl;
      converge = 99;
    }

    if(verbose) std::cout << "-- AlphaF Gradient Max: " << R_g.maxCoeff() << std::endl;

    modelobj.derivative_coef();
    modelobj.NegativeLogLikelihood();

    Vec gr_PL = modelobj.gr_phi_vec;
    Mat he_PL(kw-1, kw-1);

    // OLD:
    // Mat mat1 = modelobj.he_alpha_f_mat.ldlt().solve(modelobj.he_alpha_f_phi_mat);
    // he_PL = modelobj.he_phi_mat - mat1.transpose() * modelobj.he_alpha_f_phi_mat;

    // NEW:
    modelobj.derivative_f();
    Mat mat_tmp_PL(kE+kbetaF,kw-1);
    mat_tmp_PL.block(0, 0, kE, kw-1) = modelobj.he_alpha_f_phi_mat;
    mat_tmp_PL.block(kE, 0, kbetaF, kw-1) = modelobj.he_phi_betaF_mat.transpose();
    Mat mat1 = modelobj.he_inner_mat.ldlt().solve(mat_tmp_PL);
    he_PL = modelobj.he_phi_mat - mat1.transpose() * mat_tmp_PL;

    modelobj.PL_gradient = convertToDouble(gr_PL);
    modelobj.PL_hessian = convertToDouble(he_PL);
    modelobj.converge = converge; // 0: converge. 99: not converge
}


// function to be differentiate
Scalar logdetH05_nosmooth(Vec& alpha_f,
                 Vec& phi,
                 Vec& betaF,
                 Scalar& log_theta,
                 Scalar& log_smoothing_f,
                 Scalar& log_smoothing_w,
                 Model_nosmooth& modelobj){

    modelobj.setAlphaF(alpha_f);
    modelobj.setPhi(phi);
    modelobj.setBetaF(betaF);
    modelobj.setLogTheta(log_theta);
    modelobj.setLogSmoothingF(log_smoothing_f);
    modelobj.setLogSmoothingW(log_smoothing_w);

    modelobj.derivative_coef();
    modelobj.derivative_he();
    return modelobj.logdetH05();
}


double LAML_fn_nosmooth(Model_nosmooth& modelobj) {
    Scalar out;
    out = modelobj.logdetH05() +  modelobj.NegLogL - modelobj.n/2.0 * log(2*3.141592653589793238462643383279);
    return CppAD::Value(out);
}

struct LAMLResult_nosmooth {
    double fn;
    Eigen::VectorXd gradient;
};


void Inner_nosmooth(Model_nosmooth& modelobj, bool verbose) {
    // newton method
    int maxitr = 50, itr = 0;
    const double eps = 1e-05;
    const double largereps = 1e-03;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    // double mineig = 1e-02; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 50, stephalve = 0;
    int maxErangehalve = 20;
    int stephalve_inner = 0;
    int resetitr = 0, maxreset = 1; // if step is nan or always diverge, reset coefficients as 0. Only reset once
    int additr = 50; // allow further 20 iteraions after resetting coeffiicents as 0.

    // check non-moving step following https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
    int maxnonmovingsteps = 5; // Maximum number of iterations for which we will tolerate no movement. 5 in awstringer1/varcomptest
    double stepeps = 1e-12; // If max(abs(step)) < stepeps then we say the iteration resulted in no movement.
    int stepcounter = 0; // Count the number of non-moving steps

    Vec phi = modelobj.phi;
    Eigen::VectorXd R_phi = convertToDouble(phi);

    int kE = modelobj.kE;
    int kw = modelobj.kw;

    // catch double PL_nosmooth.fn
    double s;
    double s_tmp;

    // update steps
    Eigen::VectorXd step(kw-1);
    step.setZero();

    // Cast results
    Eigen::MatrixXd R_H(kw-1, kw-1);
    Eigen::VectorXd R_g(kw-1);

    R_g.setZero();
    R_H.setZero();

    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = R_H * eigen_in; // Easy version
      eigen_out.noalias() += R_H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, kw-1, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smalles eigenvalue

    // eigen decomposition
    Eigen::VectorXd eigvals(kw-1);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(kw-1);
    invabseigvals.setZero();
    // double eigval;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(R_H,false); // Only values, not vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(R_H,true); // Both values and vectors

    // check range of E
    Vec phi_long(kw);
    phi_long(0) = 1.0;
    Vec alpha_w_C;
    Vec E;
    Eigen::VectorXd R_E;
    int kmin;
    int kmax;
    int krange;
    int krangemin = 4;
    // int krangemin = 3;
    bool krangewarning = false;
    Mat B_inner = modelobj.getB_inner();
    Mat Dw = modelobj.getDw();
    Vec knots_f = modelobj.getknots_f();

    double delta = 1.0; // for New Q-Newton method
    double R_g_norm;

    // initialize gradient
    for (int i = 0; i < (kw-1); i++) R_g(i) = 1. + eps;

    if(verbose) std::cout << "* Start optimize profile likelihood" << std::endl;

    // start newton's method
    while (itr < maxitr) {
      itr++;
      // modelobj.NegativeLogLikelihood();

      resetcon_label: // reset phi as all zero. If the current phi always leads to the divergence in updating alpha_f
      // update alpha_f
      PL_nosmooth(modelobj, verbose);
      R_g = modelobj.PL_gradient;
      R_g_norm = R_g.norm();
      if (R_g_norm < eps) break;

      s = CppAD::Value(modelobj.NegLogL);
      R_H = modelobj.PL_hessian;


      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        // std::cout << "smallest_eigval" << smallest_eigval << std::endl;
        eigvec.compute(R_H); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();
        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < (kw-1); iii++) eigvals(iii) += delta*R_g_norm;
        }
        // for (int i = 0; i < (kw-1); i++) invabseigvals(i) = 1. / max(abs(eigvals(i)), mineig); // flip signs
        for (int i = 0; i < (kw-1); i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        // std::cout << "invabseigvals max" << invabseigvals.maxCoeff() << std::endl;
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * R_g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = R_H.llt().solve(R_g);
        step = R_H.ldlt().solve(R_g);
      }
      // check nan in step
      // NOT really needed. checking here to align with alpha_f.
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          R_phi.setZero(); // reset alpha_f
          phi = R_phi.cast<Scalar>();
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0 because of nan step" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          break;
        }
      }

      R_phi -= step;
      phi = R_phi.cast<Scalar>();

      // ****** start checking range of E
      for (int j = 0; j < (kw - 1); j++) {
        phi_long(j + 1) = phi(j);
      }
      alpha_w_C = phi_long / sqrt(phi_long.dot(Dw * phi_long));
      E = B_inner * alpha_w_C;
      R_E = convertToDouble(E);
      kmin = knotindex((Scalar) R_E.minCoeff(), knots_f);
      kmax = knotindex((Scalar) R_E.maxCoeff(), knots_f);
      krange = kmax - kmin;
      stephalve = 0;
      while ((krange < krangemin) & (stephalve < maxErangehalve)){
        stephalve++;
        step /= 2.;
        R_phi += step;
        phi = R_phi.cast<Scalar>();

        for (int j = 0; j < (kw - 1); j++) {
          phi_long(j + 1) = phi(j);
        }
        alpha_w_C = phi_long / sqrt(phi_long.dot(Dw * phi_long));
        E = B_inner * alpha_w_C;
        R_E = convertToDouble(E);
        kmin = knotindex((Scalar) R_E.minCoeff(), knots_f);
        kmax = knotindex((Scalar) R_E.maxCoeff(), knots_f);
        krange = kmax - kmin;
        if(verbose) std::cout << "E range krange" << krange << std::endl;
        // std::cout << "halving" << std::endl;
        // std::cout << "R_phi" << R_phi << std::endl;
        // std::cout << "alpha_f" << convertToDouble(modelobj.alpha_f) << std::endl;
        // Example: getLAML(c(3, 3, 3)) ## fails when kE = kw = 20. Nt = 1000. wl <- function(l) dnorm(l, mean = 10, sd = 10)/wl_de
      }
      if (stephalve >= maxErangehalve) {
        krangewarning = true;
        if(verbose) std::cout << "E range krange: " << krange << " < " << krangemin << std::endl;
        // std::cout << "Range of weighted exposure is small. Consider increasing kE and resetting starting values." << std::endl;
      }
      // finish checking ******
      modelobj.setPhi(phi);

      PL_nosmooth(modelobj, verbose);

      // halving if the optimization for alpha_f fails
      stephalve = 0;
      while ((modelobj.converge != 0) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;
        R_phi += step;
        phi = R_phi.cast<Scalar>();
        modelobj.setPhi(phi);
        PL_nosmooth(modelobj, verbose);
      }
      if (modelobj.converge != 0) {
        if (resetitr < maxreset){
          resetitr++;
          R_phi.setZero(); // reset alpha_f
          phi = R_phi.cast<Scalar>();
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0 because of divergence" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          std::cout << "Optimization for alpha_f fails" << std::endl;
          break;
        }
      }

      s_tmp = CppAD::Value(modelobj.NegLogL);

      // halving if objective function increase.
      stephalve = 0;
      while ((s_tmp > s + 1e-8) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;
        R_phi += step;
        phi = R_phi.cast<Scalar>();
        modelobj.setPhi(phi);
        PL_nosmooth(modelobj, verbose);
        // when dealing with increase: halving if the optimization for alpha_f fails
        stephalve_inner = 0;
        while ((modelobj.converge != 0) & (stephalve_inner < maxstephalve)){
          stephalve_inner++;
          step /= 2.;
          R_phi += step;
          phi = R_phi.cast<Scalar>();
          modelobj.setPhi(phi);
          PL_nosmooth(modelobj, verbose);
        }
        if (modelobj.converge != 0) {
          if (resetitr < maxreset){
            resetitr++;
            R_phi.setZero(); // reset alpha_f
            phi = R_phi.cast<Scalar>();
            modelobj.setPhi(phi);
            if(verbose) std::cout << "reset phi as 0 because of divergence" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            std::cout << "Optimization for alpha_f fails" << std::endl;
            break;
          }
        }
        s_tmp = CppAD::Value(modelobj.NegLogL);
      }

      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(s_tmp) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration
        R_phi += step;
        phi = R_phi.cast<Scalar>();
        modelobj.setPhi(phi);
        PL_nosmooth(modelobj, verbose);
        // when dealing with NaN: halving if the optimization for alpha_f fails
        stephalve_inner = 0;
        while ((modelobj.converge != 0) & (stephalve_inner < maxstephalve)){
          stephalve_inner++;
          step /= 2.;
          R_phi += step;
          phi = R_phi.cast<Scalar>();
          modelobj.setPhi(phi);
          PL_nosmooth(modelobj, verbose);
        }
        if (modelobj.converge != 0) {
          if (resetitr < maxreset){
            resetitr++;
            R_phi.setZero(); // reset alpha_f
            phi = R_phi.cast<Scalar>();
            modelobj.setPhi(phi);
            if(verbose) std::cout << "reset phi as 0 because of divergence" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            break; // break the current innner while. The next code to run: if (std::isnan(s_tmp)) {...}
          }
        }
        s_tmp = CppAD::Value(modelobj.NegLogL);
      }

      if (std::isnan(s_tmp)) {
        if (resetitr < maxreset){
          resetitr++;
          R_phi.setZero(); // reset alpha_f
          phi = R_phi.cast<Scalar>();
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0 because of nan function" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          // Step-halving didn't work
          std::cout << "Phi: Step-halving failed with nan function value. Returning failure." << std::endl;
          break;
        }
      }
      // Count the number of iterations where we didn't move; if too many, we got stuck.
      // a part of the code follows https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
      if (step.lpNorm<Eigen::Infinity>() < stepeps) {
        stepcounter++;
        if (stepcounter > maxnonmovingsteps) {
          if (resetitr < maxreset){
            resetitr++;
            R_phi.setZero(); // reset alpha_f
            phi = R_phi.cast<Scalar>();
            modelobj.setPhi(phi);
            if(verbose) std::cout << "reset phi as 0 because of non-moving" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            std::cout << "The algorithm hasn't moved for " << stepcounter << " steps; terminating. Please check the answer." << std::endl;
            break;
          }
        }
      } else {
        stepcounter = 0; // if move
      }

      // The last reset
      if((itr == maxitr) & (R_g.norm() >= largereps)) {
        if (resetitr < maxreset){
          resetitr++;
          R_phi.setZero(); // reset alpha_f
          phi = R_phi.cast<Scalar>();
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0. The last one" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          if (verbose) {
            std::cout << "Newton method for updating weight function might fail. Profile Likelihood Gradient Max: " << R_g.maxCoeff() << std::endl;
          } else {
            // report less information if verbose is false
            if (R_g.maxCoeff() >= 1e-2) {
             std::cout << "Newton method for updating weight function might fail. Profile Likelihood Gradient Max: " << R_g.maxCoeff() << std::endl;
            }
          }
          break;
        }
      }
      // initilized count
      stephalve = 0;
    }
    if(krangewarning) std::cout << "Range of weighted exposure is small. Consider increasing kE or resetting starting values. Profile Likelihood Gradient Max: " << R_g.maxCoeff() << std::endl;
    if(verbose) std::cout << "* Finish middle opt. Profile Likelihood Gradient Max: " << R_g.maxCoeff() << std::endl;
}


LAMLResult_nosmooth LAML_nosmooth(Model_nosmooth& modelobj) {
    LAMLResult_nosmooth result;
    Scalar u_LAML;
    int kE = modelobj.kE;
    int kw = modelobj.kw;
    int kbetaF = modelobj.kbetaF;
    modelobj.derivative_coef();
    modelobj.derivative_he();
    modelobj.derivative_full();
    modelobj.NegativeLogLikelihood();

    Vec alpha_f = modelobj.alpha_f;
    Vec phi = modelobj.phi;
    Vec betaF = modelobj.betaF;

    Vec gr_s_u_vec = modelobj.gr_s_u_vec;
    Vec gr_s_par_vec = modelobj.gr_s_par_vec;
    Mat he_s_u_mat = modelobj.he_s_u_mat;
    Mat he_s_par_u_mat = modelobj.he_s_par_u_mat;
    Eigen::VectorXd gr_s_u_vec_R = convertToDouble(gr_s_u_vec);
    Eigen::VectorXd gr_s_par_vec_R = convertToDouble(gr_s_par_vec);
    Eigen::MatrixXd he_s_u_mat_R = convertToDouble(he_s_u_mat);
    Eigen::MatrixXd he_s_par_u_mat_R = convertToDouble(he_s_par_u_mat);

    Scalar log_theta = modelobj.log_theta;
    Scalar log_smoothing_f = modelobj.log_smoothing_f;
    Scalar log_smoothing_w = modelobj.log_smoothing_w;

    // First derivative of LAML
    u_LAML = 0.0;
    // Cast results
    Eigen::VectorXd R_g_LAML(kE+kw+2+kbetaF);
    R_g_LAML.setZero();

    // g_LAML = gradient(logdetH05_nosmooth, wrt(alpha_f, phi, betaF, log_theta, log_smoothing_f, log_smoothing_w),
    //                              at(alpha_f, phi, betaF, log_theta, log_smoothing_f, log_smoothing_w, modelobj),
    //                              u_LAML);
    u_LAML += modelobj.NegLogL - modelobj.n/2.0 * log(2*3.141592653589793238462643383279);



    // In R: grad[-(1:(kE+kw-1))] - H.full[-(1:(kE+kw-1)),(1:(kE+kw-1))] %*% as.vector(solve(H.alpha, grad[(1:(kE+kw-1))]))
    Eigen::VectorXd g1 = R_g_LAML.segment(0, kE+kw-1+kbetaF) + gr_s_u_vec_R;
    Eigen::VectorXd g2 = R_g_LAML.segment(kE+kw-1+kbetaF, 3) + gr_s_par_vec_R;
    Eigen::VectorXd gr = g2 - he_s_par_u_mat_R * he_s_u_mat_R.ldlt().solve(g1);


    result.fn = CppAD::Value(u_LAML);
    result.gradient = gr;
    return result;
}




// [[Rcpp::export]]
List aceDLNMopt_nosmooth(const Eigen::VectorXd R_y,
                            const Eigen::MatrixXd R_B_inner,
                            const Eigen::VectorXd R_knots_f,
                            const Eigen::MatrixXd R_Sw,
                            const Eigen::MatrixXd R_Sf,
                            const Eigen::MatrixXd R_Dw,
                            const Eigen::MatrixXd R_Xfix,
                            const Eigen::MatrixXd R_Zf,
                            const Eigen::VectorXd R_Xoffset,
                            Eigen::VectorXd R_alpha_f,
                            Eigen::VectorXd R_phi,
                            double R_log_theta,
                            double R_log_smoothing_f,
                            double R_log_smoothing_w,
                            Eigen::VectorXd R_betaF,
                            bool verbose) {
    // convert
    Vec y = R_y.cast<Scalar>();
    Mat B_inner = R_B_inner.cast<Scalar>();
    Vec knots_f = R_knots_f.cast<Scalar>();
    Mat Sw = R_Sw.cast<Scalar>();
    Mat Sf = R_Sf.cast<Scalar>();
    Mat Dw = R_Dw.cast<Scalar>();
    Mat Xfix = R_Xfix.cast<Scalar>();
    Mat Zf = R_Zf.cast<Scalar>();
    Vec Xoffset = R_Xoffset.cast<Scalar>();
    Vec alpha_f = R_alpha_f.cast<Scalar>();
    Vec phi = R_phi.cast<Scalar>();
    Scalar log_theta = R_log_theta;
    Scalar log_smoothing_f = R_log_smoothing_f;
    Scalar log_smoothing_w = R_log_smoothing_w;
    Vec betaF = R_betaF.cast<Scalar>();

    // construct model
    Model_nosmooth modelobj(y, B_inner, knots_f, Sw, Sf, Dw, Xfix, Zf, Xoffset, alpha_f, phi, log_theta, log_smoothing_f, log_smoothing_w, betaF);


     // Inner opt
    Inner_nosmooth(modelobj, verbose);
    // get gr of LAML
    LAMLResult_nosmooth LAMLresult;
    LAMLresult = LAML_nosmooth(modelobj); // true: fn and gr
    return List::create(Named("LAML.fn") = LAMLresult.fn,
                        Named("LAML.gradient") = LAMLresult.gradient,
                        Named("alpha_f.mod") = convertToDouble(modelobj.alpha_f),
                        Named("phi.mod") = convertToDouble(modelobj.phi),
                        Named("betaF.mod") = convertToDouble(modelobj.betaF));
}




// [[Rcpp::export]]
List aceDLNMCI_nosmooth(const Eigen::VectorXd R_y,
                           const Eigen::MatrixXd R_B_inner,
                           const Eigen::VectorXd R_knots_f,
                           const Eigen::MatrixXd R_Sw,
                           const Eigen::MatrixXd R_Sf,
                           const Eigen::MatrixXd R_Dw,
                           const Eigen::MatrixXd R_Xfix,
                           const Eigen::MatrixXd R_Zf,
                           const Eigen::VectorXd R_Xoffset,
                           Eigen::VectorXd R_alpha_f,
                           Eigen::VectorXd R_phi,
                           double R_log_theta,
                           double R_log_smoothing_f,
                           double R_log_smoothing_w,
                           Eigen::VectorXd R_betaF,
                           const int Rci,
                           const int rseed,
                           bool ifeta,
                           bool delta,
                           bool verbose) {
  // convert
  Vec y = R_y.cast<Scalar>();
  Mat B_inner = R_B_inner.cast<Scalar>();
  Vec knots_f = R_knots_f.cast<Scalar>();
  Mat Sw = R_Sw.cast<Scalar>();
  Mat Sf = R_Sf.cast<Scalar>();
  Mat Dw = R_Dw.cast<Scalar>();
  Mat Xfix = R_Xfix.cast<Scalar>();
  Mat Zf = R_Zf.cast<Scalar>();
  Vec Xoffset = R_Xoffset.cast<Scalar>();
  Vec alpha_f = R_alpha_f.cast<Scalar>();
  Vec phi = R_phi.cast<Scalar>();
  Scalar log_theta = R_log_theta;
  Scalar log_smoothing_f = R_log_smoothing_f;
  Scalar log_smoothing_w = R_log_smoothing_w;
  Vec betaF = R_betaF.cast<Scalar>();

  // construct model
  Model_nosmooth modelobj(y, B_inner, knots_f, Sw, Sf, Dw, Xfix, Zf, Xoffset, alpha_f, phi, log_theta, log_smoothing_f, log_smoothing_w, betaF);

  int kw = modelobj.kw;
  int kE = modelobj.kE;
  int kbetaF = modelobj.kbetaF;
  int paraSize = kE+kw-1+kbetaF;
  int paraSizefull;

  // hessian
  Eigen::MatrixXd R_he;
  Eigen::VectorXd R_alpha_w(kw);

  // Vectors for sampling
  Eigen::VectorXd R_phi_sample(kw-1);
  Eigen::VectorXd R_alpha_w_sample(kw);
  Eigen::VectorXd R_alpha_f_sample(kE);
  Eigen::VectorXd R_betaF_sample(kbetaF);

  // Matrices to save results
  Eigen::MatrixXd phi_sample_mat(Rci, kw-1);
  Eigen::MatrixXd alpha_w_sample_mat(Rci, kw);
  Eigen::MatrixXd alpha_f_sample_mat(Rci, kE);
  Eigen::MatrixXd betaF_sample_mat(Rci, kbetaF);

  int n = y.size();
  // components for eta
  Vec E;
  Vec eta_sample;
  Eigen::VectorXd R_eta_sample;
  Eigen::MatrixXd eta_sample_mat;
  if(ifeta) {
    eta_sample.resize(n);
    R_eta_sample.resize(n);
    eta_sample_mat.resize(Rci, n);
  }

  // Mode of phi
  Eigen::VectorXd R_phi_mod = R_phi;
  // Generate phi
  double R_alpha_w_C_denominator;
  Eigen::VectorXd R_phi_long(kw); // phi_long = c(1, phi)
  R_phi_long(0) = 1.0;



  // d alpha_f / d phi
  Eigen::MatrixXd R_deriv_g(kw,kw-1);
  Eigen::MatrixXd R_deriv_g_large(kw,kw);
  Eigen::MatrixXd R_DiagMat_tmp(kw,kw);
  R_DiagMat_tmp.setZero();
  double R_diag_tmp;

  Eigen::MatrixXd R_deriv(paraSize+1, paraSize);
  R_deriv.setZero();
  if(delta) {
    // DELTA METHOD
    paraSizefull = paraSize+1;
    // deriv_g <- diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw) - phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    // deriv_g <- deriv_g[1:(kw), 2:kw]
    // Var_alpha_w <- deriv_g %*% Var_phi %*% t(deriv_g)
    for (int j = 0; j < (kw - 1); j++) {
      R_phi_long(j + 1) = R_phi(j);
    }

    R_diag_tmp = 1/sqrt(R_phi_long.dot(R_Dw * R_phi_long));
    R_alpha_w = R_phi_long * R_diag_tmp;

    for (int j = 0; j < kw; j++) {
      R_DiagMat_tmp(j, j) = R_diag_tmp;
    }
    R_deriv_g_large = R_DiagMat_tmp - R_phi_long * R_phi_long.transpose() * R_Dw * pow(R_phi_long.transpose() * R_Dw * R_phi_long, -3/2);
    R_deriv_g = R_deriv_g_large.block(0, 1, kw, (kw-1));

    for (int j = 0; j < paraSize; j++) {
      if (j < kE) {
        R_deriv(j,j) = 1.0;
      }
      if (j >= (kE + kw - 1)) {
        R_deriv(j+1,j) = 1.0;
      }
    }
    R_deriv.block(kE, kE, kw, kw-1) = R_deriv_g;
  } else {
    paraSizefull = paraSize;
  }



  // Joint
  R_he = convertToDouble(modelobj.he_s_u_mat);
  Eigen::VectorXd R_u_mod(paraSizefull);
  // Hessian
  // cholesky of inverse Hessian
  Eigen::MatrixXd R_he_u_L(paraSize, paraSize);
  Eigen::MatrixXd R_he_u_L_inv(paraSizefull, paraSize);
  Eigen::VectorXd zjoint(paraSize);
  Eigen::VectorXd samplejoint(paraSizefull);

  if(delta) {
    R_u_mod << R_alpha_f, R_alpha_w, R_betaF;
  } else {
    R_u_mod << R_alpha_f, R_phi, R_betaF;
  }


  // cholesky of inverse Hessian
  R_he_u_L = R_he.llt().matrixL();
  if(delta) {
    R_he_u_L_inv = R_deriv * (invertL(R_he_u_L)).transpose();
  } else {
    R_he_u_L_inv = (invertL(R_he_u_L)).transpose();
  }



  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(rseed);
  std::normal_distribution<> dist(0, 1);


  for(int i = 0; i < Rci; i++)
  {
    // Jointly sample
    for (int j = 0; j < paraSize; j++) {
      zjoint(j) = dist(gen);
    }
    samplejoint = R_u_mod + R_he_u_L_inv * zjoint;
    // get alpha_f
    R_alpha_f_sample = samplejoint.segment(0, kE);

    if(delta) {
      R_alpha_w_sample = samplejoint.segment(kE, kw);
      R_betaF_sample = samplejoint.segment(kE+kw, kbetaF);
    } else {
      // get phi
      R_phi_sample = samplejoint.segment(kE, kw-1);
      // get alpha_w
      for (int j = 0; j < (kw - 1); j++) {
        R_phi_long(j + 1) = R_phi_sample(j);
      }
      R_alpha_w_C_denominator = sqrt(R_phi_long.dot(R_Dw * R_phi_long));
      R_alpha_w_sample = R_phi_long / R_alpha_w_C_denominator;
      // get betaF
      R_betaF_sample = samplejoint.segment(kE+kw-1, kbetaF);

      // save
      phi_sample_mat.row(i) = R_phi_sample.transpose();
    }

    if(ifeta) {
      E = B_inner * R_alpha_w_sample.cast<Scalar>();
      for (int ii = 0; ii < n; ii++) {
        eta_sample(ii) = BsplinevecCon(E(ii), knots_f, 4, Zf).dot(R_alpha_f_sample.cast<Scalar>()) + Xfix.row(ii).dot(R_betaF_sample.cast<Scalar>()) + Xoffset(ii);
      }
      R_eta_sample = convertToDouble(eta_sample);
      eta_sample_mat.row(i) = R_eta_sample.transpose();
    }

    // save
    alpha_w_sample_mat.row(i) = R_alpha_w_sample.transpose();
    alpha_f_sample_mat.row(i) = R_alpha_f_sample.transpose();
    betaF_sample_mat.row(i) = R_betaF_sample.transpose();
  }

  if(delta) {
    return List::create(Named("alpha_w_sample") = alpha_w_sample_mat,
                        Named("alpha_f_sample") = alpha_f_sample_mat,
                        Named("betaF_sample") = betaF_sample_mat,
                        Named("eta_sample_mat") = eta_sample_mat,
                        Named("Hessian_inner") = R_he);
  } else {
    return List::create(Named("phi_sample") = phi_sample_mat,
                      Named("alpha_w_sample") = alpha_w_sample_mat,
                      Named("alpha_f_sample") = alpha_f_sample_mat,
                      Named("betaF_sample") = betaF_sample_mat,
                      Named("eta_sample_mat") = eta_sample_mat,
                      Named("Hessian_inner") = R_he);
  }

}







// [[Rcpp::export]]
List ConditionalAIC_nosmooth(const Eigen::VectorXd R_y,
                           const Eigen::MatrixXd R_B_inner,
                           const Eigen::VectorXd R_knots_f,
                           const Eigen::MatrixXd R_Sw,
                           const Eigen::MatrixXd R_Sf,
                           const Eigen::MatrixXd R_Dw,
                           const Eigen::MatrixXd R_Xfix,
                           const Eigen::MatrixXd R_Zf,
                           const Eigen::VectorXd R_Xoffset,
                           Eigen::VectorXd R_alpha_f,
                           Eigen::VectorXd R_phi,
                           double R_log_theta,
                           double R_log_smoothing_f,
                           double R_log_smoothing_w,
                           Eigen::VectorXd R_betaF) {
  // convert
  Vec y = R_y.cast<Scalar>();
  Mat B_inner = R_B_inner.cast<Scalar>();
  Vec knots_f = R_knots_f.cast<Scalar>();
  Mat Sw = R_Sw.cast<Scalar>();
  Mat Sf = R_Sf.cast<Scalar>();
  Mat Dw = R_Dw.cast<Scalar>();
  Mat Xfix = R_Xfix.cast<Scalar>();
  Mat Zf = R_Zf.cast<Scalar>();
  Vec Xoffset = R_Xoffset.cast<Scalar>();
  Vec alpha_f = R_alpha_f.cast<Scalar>();
  Vec phi = R_phi.cast<Scalar>();
  Scalar log_theta = R_log_theta;
  Scalar log_smoothing_f = R_log_smoothing_f;
  Scalar log_smoothing_w = R_log_smoothing_w;
  Vec betaF = R_betaF.cast<Scalar>();

  // construct model
  Model_nosmooth modelobj(y, B_inner, knots_f, Sw, Sf, Dw, Xfix, Zf, Xoffset, alpha_f, phi, log_theta, log_smoothing_f, log_smoothing_w, betaF);
  modelobj.prepare_AIC();
  // hessian
  Eigen::MatrixXd R_he;
  R_he = convertToDouble(modelobj.he_s_u_mat);
  // I
  Eigen::MatrixXd R_I;
  R_I = convertToDouble(modelobj.I_mat);
  //
  Eigen::MatrixXd mat_AIC = R_he.ldlt().solve(R_I);

  double l = CppAD::Value(modelobj.NegLogL_l);
  double edf = (2 * mat_AIC - mat_AIC * mat_AIC).trace();
  double AIC = 2.0*l + 2.0*edf;

  return List::create(Named("AIC") = AIC,
                      Named("l") = -1.0*l,
                      Named("edf") = edf);
}

