// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// full_looper
List full_looper(arma::mat locations, arma::mat yy, arma::mat knots, arma::vec previous_knot_res, arma::uword maxiters, arma::SpMat<double> XX_old, arma::SpMat<double> sigmastar_old, arma::SpMat<double> mu, double a_pi, double b_pi, double a_g, double nu, arma::vec max_dist, arma::vec res_dist, double init_log_marginal_lik, double y_ssq, arma::uword numthreads, int pi_method, int g_method, double r2_old, double m0_size, int Kernel_type_c);
RcppExport SEXP _MSSS_full_looper(SEXP locationsSEXP, SEXP yySEXP, SEXP knotsSEXP, SEXP previous_knot_resSEXP, SEXP maxitersSEXP, SEXP XX_oldSEXP, SEXP sigmastar_oldSEXP, SEXP muSEXP, SEXP a_piSEXP, SEXP b_piSEXP, SEXP a_gSEXP, SEXP nuSEXP, SEXP max_distSEXP, SEXP res_distSEXP, SEXP init_log_marginal_likSEXP, SEXP y_ssqSEXP, SEXP numthreadsSEXP, SEXP pi_methodSEXP, SEXP g_methodSEXP, SEXP r2_oldSEXP, SEXP m0_sizeSEXP, SEXP Kernel_type_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type yy(yySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type previous_knot_res(previous_knot_resSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxiters(maxitersSEXP);
    Rcpp::traits::input_parameter< arma::SpMat<double> >::type XX_old(XX_oldSEXP);
    Rcpp::traits::input_parameter< arma::SpMat<double> >::type sigmastar_old(sigmastar_oldSEXP);
    Rcpp::traits::input_parameter< arma::SpMat<double> >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type a_pi(a_piSEXP);
    Rcpp::traits::input_parameter< double >::type b_pi(b_piSEXP);
    Rcpp::traits::input_parameter< double >::type a_g(a_gSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type max_dist(max_distSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type res_dist(res_distSEXP);
    Rcpp::traits::input_parameter< double >::type init_log_marginal_lik(init_log_marginal_likSEXP);
    Rcpp::traits::input_parameter< double >::type y_ssq(y_ssqSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type numthreads(numthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type pi_method(pi_methodSEXP);
    Rcpp::traits::input_parameter< int >::type g_method(g_methodSEXP);
    Rcpp::traits::input_parameter< double >::type r2_old(r2_oldSEXP);
    Rcpp::traits::input_parameter< double >::type m0_size(m0_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type Kernel_type_c(Kernel_type_cSEXP);
    rcpp_result_gen = Rcpp::wrap(full_looper(locations, yy, knots, previous_knot_res, maxiters, XX_old, sigmastar_old, mu, a_pi, b_pi, a_g, nu, max_dist, res_dist, init_log_marginal_lik, y_ssq, numthreads, pi_method, g_method, r2_old, m0_size, Kernel_type_c));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSSS_full_looper", (DL_FUNC) &_MSSS_full_looper, 22},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSSS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
