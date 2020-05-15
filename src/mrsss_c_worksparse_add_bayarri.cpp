#include <RcppArmadillo.h>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
using namespace std;
#include <cmath>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]



double cube_rootfind(arma::vec coefs)
{
  arma::mat companion(3,3);
  companion.fill(0);
  companion(1,0)=1;
  companion(2,1)=1;
  companion(0,2)=-coefs(0)/coefs(3);
  companion(1,2)=-coefs(1)/coefs(3);
  companion(2,2)=-coefs(2)/coefs(3);
  
  double a;
  double b;
  double c;
  
  c=coefs(0)/coefs(3);
  b=coefs(1)/coefs(3);
  a=coefs(2)/coefs(3);
  
  double Q;
  double R;
  double Q3;
  
  
  Q=(pow(a,2)-3.*b)/9.;
  R=(2*pow(a,3)-9*a*b+27.*c)/54.;
  Q3=pow(Q,3);
  double disc;
  disc=pow(R,2)-Q3;
  
  double A;
  double B;
  double returner;
  if(disc>=0.){
    if(R>=0) A=-cbrt(R+sqrt(disc));
    else A=-cbrt(R-sqrt(disc));
    
    if(A==0.) B=0.;
    else B=Q/(A);
    returner=(A+B)-a/3.;
  }
  else
  {
    returner=1;
    Rcpp::Rcout <<  "2 roots, something awful happened"<< std::endl;
  }
  
  
  /*
  
  arma::Col<arma::cx_double> eigs(3);
  eigs=arma::eig_gen( companion ) ;
  
 
  returner=0;
  for(int i=0;i<(int)eigs.n_rows;++i)
  {
    returner=eigs(i).real();
    if(abs(eigs(i).imag())<.00000001)
      {
        break;
      }
  }
*/

  return(returner);
}


/*create knots based on a location and a distance*/

arma::mat create_knots1(arma::mat knot,double res_dist)
{
  arma::mat returner(2*knot.n_cols,knot.n_cols);
  returner.fill(0);
  for(int i=0;i<(int)returner.n_rows;++i)
  {
    for(int j=0;j<(int)returner.n_cols;++j)
    {
      /*    returner(i,j)=knot(0,j)+pow(-1,i+1)*/
      returner(i,j)=knot(0,j)+pow(-1,i+1)*pow(pow(-1,round((i+1)/2)),j)*.5*res_dist;
    }
  }
  return(returner);
}
/*This function is equilvalent to rdist in R, and is slower than rdist, but not by much
this could be improved using sparse matrix stuff but I can't figure out how
*/

arma::mat euclidean_cdist(arma::mat A, arma::mat B) {
  
  arma::colvec An =  sum(square(A), 1);
  arma::colvec Bn =  sum(square(B), 1);
  
  arma::mat C = -2* (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  return (sqrt(C)); 
}



/*Now we write the create_column function */

arma::sp_mat create_column1(arma::mat knot,
                         int resolution,
                         arma::mat locations,
                         arma::vec max_dist,
                         double nu,
                         int Kernel_type_c)
{
  arma::mat returner(locations.n_rows,1);
  returner.fill(0);
  returner=euclidean_cdist(locations,knot);
  double ll=0;
  ll=floor(((double)locations.n_cols)/2.0)+2.0;
  /*Bezier stuff*/
  if(Kernel_type_c==1)
  {
  
    for(int j=0;j<(int)returner.n_elem;++j)
    {
      if(returner(j)>max_dist(resolution-1))
      {
        returner(j)=0;
      }
      else
      {
        returner(j)=pow((1-pow(returner(j)/max_dist(resolution-1),2)),nu);
      }
    }
  }
  /*Wendland stuff*/
  if(Kernel_type_c==2)
  {
    
    for(int j=0;j<(int)returner.n_elem;++j)
    {
      if(returner(j)>max_dist(resolution-1))
      {
        returner(j)=0;
      }
      else
      {

        returner(j)=pow((1-returner(j)/max_dist(resolution-1)),ll)*(1+(ll+1)*returner(j)/max_dist(resolution-1));
      }
    }
  }
  return(arma::sp_mat(returner));
}

/*used in below function*/
/*stackoverflow.com/questions/40222092/access-and-modify-the-non-zero-elements-of-sparse-matrix-of-class-armasp-mat-u */
arma::umat get_locations(arma::sp_mat& B)
{
  
  // Make const iterator
  arma::sp_mat::const_iterator start = B.begin();
  arma::sp_mat::const_iterator end   = B.end();
  
  // Calculate number of points
  int n = std::distance(start, end);
  
  // Kill process if no values are found (very sparse matrix)
  if (n <= 0) { Rcpp::stop("No values found!"); }
  
  // Build a location storage matrix
  arma::umat locs(2, n);
  
  // Create a vector to store each row information in. (Row, Col)
  arma::uvec temp(2);
  
  // Start collecting locations
  arma::sp_mat::const_iterator it = start; 
  for(int i = 0; i < n; ++i)
  {
    temp(0) = it.row();
    temp(1) = it.col();
    locs.col(i) = temp;
    ++it; // increment
  }
  
  return locs;
}

/*this is slower than the above function for some reason I don't really get it*/
/*maybe reach out to stackexchange, but not priority*/

arma::SpMat<double> create_column2(arma::mat knot,
                                   int resolution,
                                   arma::mat locations,
                                   arma::vec max_dist,
                                   double nu,
                                   arma::SpMat<double> parent_col)
{
  arma::umat imp_ind=get_locations(parent_col);
  arma::mat imp_loc=locations.rows(imp_ind.row(0));
  arma::mat temp=euclidean_cdist(imp_loc,knot);
  temp=pow((1-pow(temp/max_dist(resolution-1),2)),nu);
  /*  fix fractional powers*/
  temp.replace(arma::datum::nan,0);
  temp=clamp(temp,0,temp.max());
  arma::SpMat<double> returner(imp_ind,temp.col(0),parent_col.n_rows,1);
  return(returner);
}



/*Code sourced from BAS within bayesreg.c, quite obfuscated*/ 

double log_laplace_2F1(double a, double b, double c, double z)
{
  
  double ghat1, ghat,logint,sigmahat, D, A, B, C;
  
  // 2F1(a,b,c, Z)      =        Gamma(c)
  //                        --------------------- * integral
  //			     Gamma(b) Gamma(c - b)
  
  // integral = \int_0^\infty  g^{b -1} (1 + g)^{c - b - 1}(1 + (1 - z)g)^{-a}
  
  logint = 0.0;
  
  if (c >= b && b > 0) {
    logint = lgamma(c) - lgamma(b) - lgamma(c - b);
  }
  
  if ( z == 1.0) {
    if (c > b + a)
      logint =  lgamma(c) + lgamma(c - a - b) - lgamma(c - b)
      - lgamma(c - a);
    else logint = log(0.0);
  }
  else{
    /*  Laplace approximation in terms of exp(tau) = g  */
    //
    //  Integrate[g^(b-1) (1 + g)^(a - c) (1 + (1 - z) g)^(-a) dg
    //  Integrate[e^tau b (1 + e^g)^(a - c) *( 1 + (1 -z)e^g)^(-a) dtau
    //  Jacobian e^tau
    
    A = (b - c)*(1. - z);
    B = 2.*b - c + (a-b)*z;
    C = b;
    D = B*B - 4.*A*C;
    
    if (D < 0 )    Rprintf("ERROR: In Laplace approximation to 2F1");
    
    /*  root1 = (-B - sqrt(D))/(2.*A);*/
    /*  root2 = (-B + sqrt(D))/(2.*A);*/
    
    ghat1 = (2.*b)/(c + b*(z - 2.0) - a*z + sqrt(c*c + 4*a*b*z -  2.0*(a + b)*c *z + (a -b )*(a-b)*z*z));
    
    /*    ghat2 =  -(c +b*(-2. + z) - a*z + sqrt(c*c + 4*a*b*z -  2.0*(a + b)*c *z + (a -b )*(a-b)*z*z))/(2.*(b - c)*(z - 1.));*/
    
    // Rprintf("+root = %lf -root = %lf ghat1  = %lf hgat2 = %lf\n", root2, root1, ghat1, ghat2);
    ghat = ghat1;
    if ( ghat1 < 0) Rprintf("ERROR: In Laplace approximation to 2F1");
    
    sigmahat =1.0/((-a + c)*((ghat/(1 + ghat))*(1 - ghat/(1 + ghat))) +
      a*((ghat*(1.-z))/(1.+ghat*(1.-z))*
      (1. - (ghat*(1.-z))/(1.+ghat*(1.-z)))));
    
    if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to in 2F1 sigmhat = %f, ghat =  %f  z = %f \n", sigmahat, ghat, z);
    logint = logint
      + .5*( log(2.0*PI) +  log(sigmahat)) +
        b*log(ghat) + (a - c)*log(1.0 + ghat)
      -  a*log(1.0 + ghat*(1. - z));
  }
  return(logint);
}



/*log binomial distribution*/
double lbinom(double x, double size, double prob)
{
  double returner;
  returner=x*log(prob)+(size-x)*log(1-prob) + lgamma(size+1)-lgamma(x+1)-lgamma(size-x+1);
  return(returner);
}
/*calculates prior (branching process) for our knots*/

double calculate_prior1(arma::vec resolutions,double a_pi,double b_pi,int dim,int method)
{
  double returner=0;
  /*This is for beta binomial*/
  if(method==1)
  {
    arma::vec tabler(max(resolutions)+1);
    tabler.fill(0);
    for(int i=0;i<(int)max(resolutions);++i)
    {
      for(int j=0;j<(int)resolutions.n_elem;++j)
      {
        if(resolutions(j)==i+1)
        {
          tabler(i)=tabler(i)+1;
        }
      }
    }
    double a_new=a_pi+accu(tabler)-tabler(0);
    double b_new=b_pi;
    for(int i=1;i<(int)tabler.n_elem;++i)
    {
      b_new=b_new+tabler(i-1)*pow(2,dim)-tabler(i);
    }
    returner=R::lbeta(a_new,b_new)-R::lbeta(a_pi,b_pi);
  }
  /*this is for fixed pi*/
  /*note that a_pi and b_pi are the numerator and denominator of the fraction if method==2*/
  if(method==2)
  {
    double pii;
    pii=a_pi/b_pi;
    arma::vec tabler(max(resolutions)+1);
    tabler.fill(0);
    for(int i=0;i<(int)max(resolutions);++i)
    {
      for(int j=0;j<(int)resolutions.n_elem;++j)
      {
        if(resolutions(j)==i+1)
        {
          tabler(i)=tabler(i)+1;
        }
      }
    }
    double log_prior_prob=0;
    for(int i=1;i<(int)tabler.n_elem;++i)
    {
      log_prior_prob=log_prior_prob+lbinom(tabler(i),tabler(i-1)*pow(2,dim),pii);
    }
    returner=log_prior_prob;
  }
  
  
  return(returner);
}


struct plus_marginal_holder
{
  arma::mat plus_marginals;
  arma::mat mus;
};
/*helper function for below*/

double hhh(double g, double n, double p, double r2)
{
  double returner;
  returner=((n-p-1)/2)*log(1+g)+(-(n-1)/2)*log(1+(1-r2)*g)+(-1.5)*log(g)-n/(2*g);
  return(returner);
}
double hhhd(double g, double n, double p, double r2,double D)
{
  double returner;
  returner=((n-p-1)/2)*log(1+g)+(-(n-1)/2)*log(1+(1-r2)*g)+(-1.5)*log(g)-n/(2*D*g);
  return(returner);
}
double hhhg(double g, double n, double p, double r2,double D,double a, double b)
{
  double returner;
  returner=((n-p-1+2*b)/2)*log(1+g)+(-(n-1)/2)*log(1+(1-r2)*g)+(a-1.5)*log(g)-D/(2*g);
  return(returner);
}

/*laplace approximation for ZS marginal*/

double zs_marg(arma::uword nn, arma::uword pp, double r2)
{
  double rooter;
  double sigmahat;
  double lapapprox;
  double loginter;
  
  arma::vec coefs(4);
  coefs(3)=-(1-r2)*(pp+3);
  coefs(2)=(nn-pp-4-2*(1-r2));
  coefs(1)=(nn*(2-r2)-3);
  coefs(0)=nn;
  rooter=cube_rootfind(coefs);
  sigmahat=pow(-.5*(
    .5*(
    ((nn-1)*pow(1-r2,2))/(pow((1+rooter*(1-r2)),2))-
    (nn-pp-1)/(pow(1+rooter,2))+
    (3)/(pow(rooter,2))-
    2*nn/(pow(rooter,3)))),-.5);
  loginter=hhh(rooter,nn,pp,r2);
  
  lapapprox=log(pow(2*arma::datum::pi,.5))+log(sigmahat)+(loginter);

  return(lapapprox);
}

double zs_margsig(arma::uword nn, arma::uword pp, double r2)
{
  double rooter;
  double sigmahat;
  double lapapprox;
  double loginter;
  
  arma::vec coefs(4);
  coefs(3)=-(1-r2)*(pp+3);
  coefs(2)=(nn-pp-4-2*(1-r2));
  coefs(1)=(nn*(2-r2)-3);
  coefs(0)=nn;
  rooter=cube_rootfind(coefs);
  sigmahat=pow(-.5*(
                (
  ((nn-1.)*pow(1.-r2,2))/(pow((1.+rooter*(1.-r2)),2))-
  (nn-pp-1)/(pow(1.+rooter,2))+
  (3.)/(pow(rooter,2))-
  2.*nn/(pow(rooter,3)))),-.5);
  loginter=hhh(rooter,nn,pp,r2);
  
  Rcpp::Rcout <<  log(pow(2.*arma::datum::pi,.5))<< std::endl;
  Rcpp::Rcout <<  log(sigmahat)<< std::endl;
  Rcpp::Rcout <<  log(loginter)<< std::endl;
  
  
  lapapprox=log(pow(2.*arma::datum::pi,.5))+log(sigmahat)+(loginter);
  
  return(lapapprox);
}
/*laplace approximation for ZS marginal with linearly penalized n in mixing dist*/

double zsd_marg(arma::uword nn, arma::uword pp, double r2,double D)
{
  double rooter;
  double sigmahat;
  double lapapprox;
  double loginter;
  
  arma::vec coefs(4);
  coefs(3)=-(1-r2)*(pp+3);
  coefs(2)=(nn-pp-4+(nn*(1-D)/D-2)*(1-r2));
  coefs(1)=((nn/D)*(2-r2)-3);
  coefs(0)=nn/D;
  rooter=cube_rootfind(coefs);
/*  Rcpp::Rcout << rooter << std::endl;*/
  sigmahat=pow(-.5*(
                .5*(
  ((nn-1)*pow(1-r2,2))/(pow((1+rooter*(1-r2)),2))-
  (nn-pp-1)/(pow(1+rooter,2))+
  (3)/(pow(rooter,2))-
  2*nn/(D*pow(rooter,3)))),-.5);
  loginter=hhhd(rooter,nn,pp,r2,D);
/*  Rcpp::Rcout << sigmahat << std::endl;*/
 /* Rcpp::Rcout << loginter << std::endl;*/
  lapapprox=log(pow(2*arma::datum::pi,.5))+log(sigmahat)+(loginter);
  
  return(lapapprox);
}

/*laplace approximation for ZS modified marginal with D/2 as the beta parameter for the mixing dist*/
/*aa and bb are the a and b from the paper, indicate which integral we're solving, a=b=0 is the marginal likelihood*/
 /*D=n gets us zellner siow */

double zsg_marg(arma::uword nn, arma::uword pp, double r2,double D,double aa, double bb)
{
  double rooter;
  double sigmahat;
  double lapapprox;
  double loginter;
  
  arma::vec coefs(4);
  coefs(3)=-(1-r2)*(pp+3-2*(aa+bb));
  coefs(2)=(nn-pp-4+2*aa+2*bb+(1-r2)*(D+2*aa-nn-2));
  coefs(1)=(D*(1-r2)+D+2*aa-3);
  coefs(0)=D;
  rooter=cube_rootfind(coefs);
  /* Rcpp::Rcout << rooter << std::endl;*/
  sigmahat=pow(-.5*(
                .5*(
  ((nn-1)*pow(1-r2,2))/(pow((1+rooter*(1-r2)),2))-
  (nn-pp-1)/(pow(1+rooter,2))+
  (3)/(pow(rooter,2))-
  2*D/(pow(rooter,3)))),-.5);
  loginter=hhhg(rooter,nn,pp,r2,D,aa,bb);
  /* Rcpp::Rcout << sigmahat << std::endl;*/
   /*Rcpp::Rcout << loginter << std::endl;*/
  lapapprox=log(pow(2*arma::datum::pi,.5))+log(sigmahat)+(loginter);
  
  return(lapapprox);
}




double bic_marg(arma::uword nn, arma::uword pp, double r2)
{
  double mle_g;
  double negBIC;
  
  mle_g=(nn*r2-r2-pp)/(pp*(1-r2));
  negBIC=-log(nn)*pp+2*(
    ((nn-pp-1)/2)*log(1+mle_g)-
    ((nn-1)/2)*log(1+(1-r2)*mle_g)
    );
  negBIC=negBIC/2;
  

  
  return(negBIC);
}


double bruno_marg(arma::uword nn, arma::uword pp, double r2, double yssq, double ybar)
{
  double rss;
  double del;
  double returner;
  rss=yssq*(1.0-r2);
  del=(r2*yssq+(double)nn*pow(ybar,2))/(double)nn;
  returner=lgamma(((double)pp+1)/2)+lgamma(((double)nn-(double)pp)/2.0)-(((double)pp+1.0)/2.0)*log((del+1.0)/rss)-(((double)nn-1.0)/2.0)*log(rss);
  return(returner);
}



double bic_wangli_marg(arma::uword nn, arma::uword pp, double r2)
{
  double mle_g;
  double negBIC;
  
  mle_g=(nn*r2-r2-pp)/(pp*(1-r2));
  negBIC=-sqrt(nn)*pp+2*(
    ((nn-pp-1)/2)*log(1+mle_g)-
    ((nn-1)/2)*log(1+(1-r2)*mle_g)
  );
  negBIC=negBIC/2;
  
  
  
  return(negBIC);
}

/*Now we write the plus_marginal_fun function which updates a single plus marginal*/
/*can no long rcpp export it since I'm using a custom structure*/
plus_marginal_holder plus_marginal_fun1(arma::mat knot,
                                        arma::mat current_knots,
                                        int resolution,
                                        arma::mat locations,
                                        arma::vec max_dist,
                                        double nu,
                                        arma::SpMat<double> XX_old,
                                        arma::SpMat<double> sigmastar_old,
                                        arma::SpMat<double> mu,
                                        arma::mat yy,
                                        double y_ssq,
                                        arma::vec res_dist,
                                        double a_g,
                                        arma::vec previous_knot_res,
                                        double a_pi,
                                        double b_pi,
                                        int pi_method, int g_method, double r2_old, double m0_size, int Kernel_type_c)
{
  arma::mat pot_knots(2*knot.n_cols,knot.n_cols);
  pot_knots.fill(0);
  /*create potential knots from previous knot*/
  pot_knots=create_knots1(knot,res_dist(resolution));

  /*remove from consideration knots that are already in our model by forcing their log likelihoods to be NaN*/
  
  arma::vec find_elem(pot_knots.n_rows);
  find_elem.fill(0);
  for(int j=0;j<(int)pot_knots.n_rows;++j)
  {
    for(int k=0;k<(int)current_knots.n_rows;++k)
    {
      find_elem(j)=0;
      for(int l=0;l<(int)pot_knots.n_cols;++l)
      {
        find_elem(j)=find_elem(j)+abs(current_knots(k,l)-pot_knots(j,l));
      }
      /*this is pretty terrible, think of better maybe, can't use equal because of floating point stuff*/
      if(find_elem(j)<=.00000000001)
      {
        find_elem(j)=0;
        break;
      }
    }
  }
  /*if find_elem(j) is zero then we don't want that knot included,NaN*/
  
  arma::mat returner(pot_knots.n_rows,knot.n_cols+1);
  /*need to initialize mus, should be as many rows as there are potential knots*/
  /*needs max_knot number of columns, which is 100*(res1 knot count) so must calculate that*/
  int counter=0;
  for(int i=0;i<(int)previous_knot_res.n_elem;++i)
  {
    if(previous_knot_res(i)==1) 
    {
      counter=counter+1;
    }
  }
  arma::mat mus(pot_knots.n_rows,100*counter);
  
  returner.fill(arma::datum::nan);
  mus.fill(0);
  returner.submat(0,0,pot_knots.n_rows-1,knot.n_cols-1)=pot_knots;
  plus_marginal_holder return_fun;
  return_fun.plus_marginals=returner;
  return_fun.mus=mus;
  for(int i=0;i<(int)pot_knots.n_rows;++i)
  {
    /*Here's where they get coded to nan*/
    double log_marginal_lik;
    if(find_elem(i)==0)
    {
      log_marginal_lik=arma::datum::nan;
    }
    else
    {
      /*create new column*/
      arma::sp_mat XX_new_holder(locations.n_rows,1);
      arma::sp_mat XX_new(yy.n_rows,1);
      XX_new_holder=create_column1(pot_knots.row(i),resolution+1,locations,max_dist,nu,Kernel_type_c);
      XX_new.submat(0,0,locations.n_rows-1,0)=XX_new_holder;
      /*deal with empty columns*/
      if(XX_new.n_nonzero<5)
      {
        returner(i,knot.n_cols)=arma::datum::nan;
        continue;
      }
      /*ends up as NaN for log likelihood for empty xx_new*/
      log_marginal_lik=0;

      double qstar11=accu(square(XX_new));
      arma::mat qstar1;
      /*following add 1 regression update, using notation and work from peter mueller wavelet notes*/
      qstar1=XX_old.t()*(XX_new);
      
      arma::mat hh=sigmastar_old*qstar1;
      arma::mat new_temp(1+sigmastar_old.n_rows,1+sigmastar_old.n_cols);
      new_temp.fill(0);
      new_temp(0,0)=1;
      new_temp.submat(1,1,new_temp.n_cols-1,new_temp.n_cols-1)=hh*hh.t();
      new_temp.submat(0,1,0,new_temp.n_cols-1)=-hh.t();
      new_temp.submat(1,0,new_temp.n_rows-1,0)=-hh;
      arma::mat sigmastar_new(sigmastar_old.n_rows+1,sigmastar_old.n_cols+1);
      sigmastar_new.fill(0);
      sigmastar_new.submat(1,1,sigmastar_new.n_cols-1,sigmastar_new.n_cols-1)=sigmastar_old;
      sigmastar_new=sigmastar_new+sigmastar_new+new_temp * (1/as_scalar((qstar11-(qstar1.t())*hh)));
      
      arma::SpMat<double> mu_possible(mu.n_rows+1,1);
      /* mu_possible.fill(0);*/
      mu_possible.submat(1,0,mu.n_rows,0)=mu;
      arma::SpMat<double> new_design=join_rows(arma::sp_mat(XX_new),XX_old);
      
      
      mu_possible=mu_possible+(1/as_scalar((qstar11-(qstar1.t())*hh)))*new_temp.col(0)*(new_temp.row(0)*new_design.t()*yy);
      mus.submat(i,0,i,mu.n_rows)=mu_possible.t();
      double r2_new;
      /*no avoiding this costly update but sparse matrix multiplication isn't as costly as regular matrix mult*/
      r2_new=1-accu(pow(yy-new_design*mu_possible,2))/y_ssq;
      
      /*for now automatically hyper g*/    /*for now automatically hyper g*/    /*for now automatically hyper g*/
      /*Bayarri and berger*/
      if(g_method==1)
      {
        log_marginal_lik=log(2*a_g)-log((double)mu_possible.n_rows-m0_size+2*a_g)+log_laplace_2F1(((double)yy.n_rows-m0_size)/2,1,a_g+1+((double)mu_possible.n_rows-m0_size)/2, 1-(1-r2_new)/(1-r2_old));
      }
      /*FOR HYPER G PRIOR*/
      /*Blatantly stolen from BAS within bayesreg.c, quite obfuscated to get hyper g prior marginal*/
      /*FOR ZELLNER SIOW PRIOR, uses laplace approximation, only a function of n, r2, and p*/
      if(g_method==2)
      {
        log_marginal_lik=(log(a_g-2)-log(mu_possible.n_rows+a_g-2))+log_laplace_2F1(((double)yy.n_rows-1)/2,1,((double)mu_possible.n_rows+a_g)/2, r2_new);
      }
      /*DK edit of zellner siow, just adds a parameter to linearly penalize n*/
      if(g_method==3)
      {
        log_marginal_lik=zsg_marg(yy.n_elem,mu_possible.n_rows,r2_new,yy.n_elem/a_g,0,0);
      }
      /*wangli bic*/
      if(g_method==4)
      {
/*        log_marginal_lik=bic_wangli_marg(yy.n_elem, mu_possible.n_rows, r2_new);*/
        log_marginal_lik=bruno_marg(yy.n_elem, mu_possible.n_rows, r2_new, y_ssq, accu(yy)/yy.n_elem);
      }
      arma::vec temp(1);
      temp=resolution+1;
      arma::mat holder(previous_knot_res.n_elem+1,1);
      holder=join_cols(previous_knot_res,temp);
      /*add in prior*/
      log_marginal_lik=log_marginal_lik+calculate_prior1(holder,a_pi,b_pi,locations.n_cols,pi_method);
    }
    returner(i,knot.n_cols)=log_marginal_lik;
    
    return_fun.plus_marginals=returner;
    return_fun.mus=mus;
  }
  return(return_fun);
}
/*this code applies the above function in parallel*/

plus_marginal_holder plus_marginal_fun_parallel1(arma::mat knots,
                                                 arma::mat locations,
                                                 arma::vec max_dist,
                                                 double nu,
                                                 arma::SpMat<double> XX_old,
                                                 arma::SpMat<double> sigmastar_old,
                                                 arma::SpMat<double> mu,
                                                 arma::mat yy,
                                                 double y_ssq,
                                                 arma::vec res_dist,
                                                 double a_g,
                                                 arma::vec previous_knot_res,
                                                 double a_pi,
                                                 double b_pi,
                                                 int numthreads,
                                                 int pi_method, int g_method, double r2_old, double m0_size, int Kernel_type_c)
{
  /*need to initialize mus, should be as many rows as there are potential knots*/
  /*needs max_knot number of columns, which is 100*(res1 knot count) so must calculate that*/
  int counter=0;
  for(int i=0;i<(int)previous_knot_res.n_elem;++i)
  {
    if(previous_knot_res(i)==1) 
    {
      counter=counter+1;
    }
  }
  
  
  arma::mat returner(2*knots.n_cols*knots.n_rows,knots.n_cols+1);
  arma::mat mu_returner(2*knots.n_cols*knots.n_rows,100*counter);
  returner.fill(0);
  mu_returner.fill(0);
  /*talk to someone if there's a better pragma call here*/
  omp_set_num_threads(numthreads);
#pragma omp parallel for shared(returner,mu_returner)
  for(int i=0;i<(int)knots.n_rows;++i)
  {
    plus_marginal_holder return_fun=plus_marginal_fun1(knots.row(i),
                                                       knots,
                                                       previous_knot_res(i),
                                                       locations,
                                                       max_dist,
                                                       nu,
                                                       XX_old,
                                                       sigmastar_old,
                                                       mu,
                                                       yy,
                                                       y_ssq,
                                                       res_dist,
                                                       a_g,
                                                       previous_knot_res,
                                                       a_pi,
                                                       b_pi,
                                                       pi_method, g_method,r2_old, m0_size, Kernel_type_c);




    
    returner.submat(i*2*knots.n_cols,0,(i+1)*knots.n_cols*2-1,knots.n_cols)=return_fun.plus_marginals;
    mu_returner.submat(i*2*knots.n_cols,0,(i+1)*knots.n_cols*2-1,100*counter-1)=return_fun.mus;
    
  }
  plus_marginal_holder return_holder;
  return_holder.mus=mu_returner;
  return_holder.plus_marginals=returner;
  return return_holder;
}


struct minus_marginal_holder
{
  double log_l;
  arma::mat mus;
};

/*now for minus_marginals*/

minus_marginal_holder minus_marginal_fun1(   /*this is actually xx_old but we're going to modify it in place*/
arma::SpMat<double> XX_new,
int colnum_delete,
arma::SpMat<double> sigmastar_old,
/*this is actually old mu but we're going to modify it in place*/
arma::SpMat<double> mu_new,
arma::mat yy,
/*need this explicitly because we're not passing the locations to this function*/
int dim,
double y_ssq,
double a_g,
arma::vec previous_knot_res,
double a_pi,
double b_pi,int pi_method, int g_method, double r2_old, double m0_size)
{
  int counter=0;
  for(int i=0;i<(int)previous_knot_res.n_elem;++i)
  {
    if(previous_knot_res(i)==1) 
    {
      counter=counter+1;
    }
  }
  arma::mat mus(1,100*counter);
  
  XX_new.shed_col(colnum_delete);
  
  
  double sigma_ll=sigmastar_old(colnum_delete,colnum_delete);
  arma::SpMat<double> sigma_l=  sigmastar_old.col(colnum_delete);
  
  sigma_l.shed_row(colnum_delete);
  
  
  double old_muu=mu_new(colnum_delete,0);
  mu_new.shed_row(colnum_delete);
  mu_new=mu_new-sigma_l*old_muu/sigma_ll;
  mus.submat(0,0,0,mu_new.n_elem-1)=mu_new.t();
  double r2_new=1-accu(pow(yy-XX_new*mu_new,2))/y_ssq;
  
  /*Blatantly stolen from BAS within bayesreg.c, quite obfuscated to get hyper g prior marginal*/
  double log_marginal_lik=0;
  
  /*Blatantly stolen from BAS within bayesreg.c, quite obfuscated to get hyper g prior marginal*/
  /*Bayarri and berger*/
  if(g_method==1)
  {
    /*log_marginal_lik=(log(a_g-2)-log(mu_new.n_rows+a_g-2))+log_laplace_2F1(((double)yy.n_rows-1)/2,1,((double)mu_new.n_rows+a_g)/2, r2_new);*/
    log_marginal_lik=log(2*a_g)-log((double)mu_new.n_rows-m0_size+2*a_g)+log_laplace_2F1(((double)yy.n_rows-m0_size)/2,1,a_g+1+((double)mu_new.n_rows-m0_size)/2, 1-(1-r2_new)/(1-r2_old));
  }
  
  /*Hyper G*/
  if(g_method==2)
  {
    log_marginal_lik=(log(a_g-2)-log(mu_new.n_rows+a_g-2))+log_laplace_2F1(((double)yy.n_rows-1)/2,1,((double)mu_new.n_rows+a_g)/2, r2_new);
  }
  /*DK edit of zellner siow, just adds a paremeter to linearly penalize n*/
  if(g_method==3)
  {
    log_marginal_lik=zsg_marg(yy.n_elem,mu_new.n_rows,r2_new,yy.n_elem/a_g,0,0);
  }
  /*increase df of t*/
  if(g_method==4)
  {
    /*log_marginal_lik=bic_wangli_marg(yy.n_elem, mu_new.n_rows, r2_new);*/
    log_marginal_lik=bruno_marg(yy.n_elem, mu_new.n_rows, r2_new, y_ssq, accu(yy)/yy.n_elem);
  }
  
  
  /*previous_knot_res will be updated in the same order that XX_old is in, be very careful*/
  previous_knot_res.shed_row(colnum_delete);
  /*add in prior*/
  log_marginal_lik=log_marginal_lik+calculate_prior1(previous_knot_res,a_pi,b_pi,dim,pi_method);
  minus_marginal_holder returner;
  returner.log_l=log_marginal_lik;
  returner.mus=mus;
  return returner;
}

/*this code applies the above function in parallel*/

plus_marginal_holder minus_marginal_fun_parallel1(arma::SpMat<double> XX_old,
                                       arma::SpMat<double> sigmastar_old,
                                       arma::SpMat<double> mu,
                                       arma::mat yy,
                                       double y_ssq,
                                       double a_g,
                                       arma::vec previous_knot_res,
                                       double a_pi,
                                       double b_pi,
                                       arma::uword numthreads,
                                       arma::uvec childless,
                                       int dim,
                                       int pi_method, int g_method, double r2_old, double m0_size)
{
  int counter=0;
  for(int i=0;i<(int)previous_knot_res.n_elem;++i)
  {
    if(previous_knot_res(i)==1) 
    {
      counter=counter+1;
    }
  }
  /*create vector childless for knots that can be deleted*/
  
  arma::mat returner(childless.n_elem,2);
  arma::mat mus(childless.n_elem,100*counter);
  mus.fill(0);
  int supernum=min(numthreads,childless.n_elem);
  /*talk to someone if there's a better pragma call here*/
  omp_set_num_threads(supernum);
#pragma omp parallel for shared(returner)
  for(int i=0;i<(int)childless.n_elem;++i)
  {
    minus_marginal_holder temp=minus_marginal_fun1(XX_old,
                                                   childless(i),
                                                   sigmastar_old,
                                                   mu,
                                                   yy,
                                                   dim,
                                                   y_ssq,
                                                   a_g,
                                                   previous_knot_res,
                                                   a_pi,
                                                   b_pi,
                                                   pi_method, g_method,r2_old,m0_size);
    returner(i,1)=temp.log_l;
    returner(i,0)=childless(i);
    mus.row(i)=temp.mus;

  }
  plus_marginal_holder return_temp;
  return_temp.mus=mus;
  return_temp.plus_marginals=returner;
  return(return_temp);
}


IntegerVector dan_sample(IntegerVector population, arma::vec probablities)
{
  probablities=probablities/accu(probablities);
  probablities=cumsum(probablities);
  IntegerVector returner(1);
  NumericVector randu=Rcpp::runif(1,0,1);
  for(int i=0;i<(int)probablities.n_elem;++i)
  {
    if(randu(0)<probablities(i))
    {
      returner(0)=population(i);
      break;
    }
  }
  return returner;
}

arma::uvec my_setdiff(arma::uvec& x, const arma::uvec& y){
  
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y[j]));
    x.shed_row(q1);
  }
  return x;
}
// this is the only function that must be exported
// [[Rcpp::export(.full_looper)]]
List full_looper(arma::mat locations,
                 arma::mat yy,
                 arma::mat knots,
                 arma::vec previous_knot_res,
                 arma::uword maxiters,
                 arma::SpMat<double> XX_old,
                 arma::SpMat<double> sigmastar_old,
                 arma::SpMat<double> mu,
                 double a_pi,
                 double b_pi,
                 double a_g,
                 double nu,
                 arma::vec max_dist,
                 arma::vec res_dist,
                 double init_log_marginal_lik,
                 double y_ssq,
                 arma::uword numthreads,
                 int pi_method,
                 int g_method,
                 double r2_old,
                 double m0_size,
                 int Kernel_type_c)
{
  
  /*storage for top models by likelihood and their knots*/
  arma::vec top_models_log_likelihood(300);
  
  arma::vec top_model_storage(maxiters);
  top_model_storage.fill(0);
  arma::vec current_model_r2(maxiters);
  current_model_r2.fill(0);
  arma::vec top_model_size(maxiters);
  top_model_size.fill(0);
  arma::mat top_model_knot_res(300,100*knots.n_rows);
  arma::mat top_model_mus(300,100*knots.n_rows);
  arma::mat top_model_knot_res_iters(100*knots.n_rows,maxiters);
  top_model_knot_res.fill(0);
  top_model_mus.fill(0);
  int iterator;
  /*this will break if the model tries to return 10 times as many knots as it started with, should almost never happen*/
  arma::cube top_models_knots(100*knots.n_rows,knots.n_cols,300);

  top_models_knots.fill(arma::datum::nan);
  /*for debugging */
  arma::mat plus_marginals;
  arma::mat plus_mus;
  arma::mat minus_marginals;
  arma::uvec childless_knot_columns;
  /*this is storage for the current set of knots' parent columns in XX*/
  /*R1 knots will be set as their own parent*/
  arma::uvec parent_knot_columns(knots.n_rows);
  
  for(int i=0;i<(int)knots.n_rows;++i)
  {
    parent_knot_columns(i)=i;
  }
  arma::uvec plus_top100_indicies;
  arma::uvec minus_top100_indicies;
  
  IntegerVector add_knot(1);
  add_knot(0)=1; 
  for(int i=0;i<(int)maxiters;++i)
  {

    add_knot(0)=1;
   
    /*initialize storage for stuff*/

    IntegerVector minus_row;
    IntegerVector plus_row;
    plus_marginal_holder temp;

    temp=plus_marginal_fun_parallel1(knots,locations,max_dist,nu,XX_old,sigmastar_old,mu,
                                     yy,y_ssq,res_dist,a_g,previous_knot_res,a_pi,b_pi,numthreads,pi_method, g_method,r2_old,m0_size, Kernel_type_c);


    plus_marginals=temp.plus_marginals;
    plus_mus=temp.mus;
    arma::vec parent_col(plus_marginals.n_rows);
    
    for(int b=0;b<(int)plus_marginals.n_rows;++b)
    {
      parent_col(b)=floor(b/(2*locations.n_cols));
    }
    

    
    /*exclude NAN to deal with empty columns / useless knots*/
    for(int k=plus_marginals.n_rows;k>(int)0;--k)
    {
      /*recent fix, knots.n_cols is the spatial dimension, that's the column that has the log likelihood*/
      if(std::isnan(plus_marginals(k-1,knots.n_cols)))
      {
        plus_marginals.shed_row(k-1);
        plus_mus.shed_row(k-1);
        parent_col.shed_row(k-1);
      }
    }

    
    
    /*Now sample a new plus one model a la SSS*/
    IntegerVector plusrownums=seq_len(plus_marginals.n_rows);
    arma::vec plusprobabilities=(exp(plus_marginals.col(knots.n_cols)-max(plus_marginals.col(knots.n_cols))));
    plus_row=dan_sample(plusrownums,plusprobabilities);
    plus_row(0)=plus_row(0)-1;

    if((i==0))
    {
      top_models_log_likelihood.fill((min(plus_marginals.col(knots.n_cols))-1));
    }
    /*now we grab the top 100 plus models, or fewer if there are fewer to choose from*/
    int plus_topn=min(100,(int)plus_marginals.n_rows);
    
    arma::vec plus_log_likelihoods=plus_marginals.col(knots.n_cols);
    /*   arma::uvec plus_top100_indicies;*/
    plus_top100_indicies=sort_index(plus_log_likelihoods,"descend");
    /*this now contains the indicies of the top 100 new plus models (or fewer)*/
    
    plus_top100_indicies=plus_top100_indicies.head(plus_topn);
    
    
    
    /*now do the same for the knots of these models, put them in the slices*/
    /*first preinitilize a matrix that has the knots, but is padded with NA, including an empty first element*/
    arma::mat knotsplus(top_models_knots.n_rows,top_models_knots.n_cols);
    /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/
    
    arma::vec resplus(top_model_knot_res.n_cols);
    resplus.fill(0);
    /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/
    resplus.subvec(1,previous_knot_res.n_elem)=previous_knot_res;
    knotsplus.fill(arma::datum::nan);
    knotsplus.rows(1,knots.n_rows)=knots;
 
    /*now loop through and assign the slice, reassigning the first element to the knot*/
    for(int l=0;l<(int)plus_topn;++l)
    {
      top_models_knots.slice(100+l)=knotsplus;
      top_models_knots(0,0,100+l)=plus_marginals(plus_top100_indicies(l),0);
      if(knots.n_cols==2)
      {
        top_models_knots(0,1,100+l)=plus_marginals(plus_top100_indicies(l),1);
      }
      /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/ /*NEW*/
      resplus(0)=previous_knot_res(parent_col(plus_top100_indicies(l)))+1;
      top_model_knot_res.row(100+l)=resplus.t();
      top_model_mus.row(100+l)=plus_mus.row(plus_top100_indicies(l));
      top_models_log_likelihood(100+l)=plus_marginals(plus_top100_indicies(l),knots.n_cols);
    }
    arma::uvec all_parent_cols=unique(parent_knot_columns);
    arma::uvec all_knot_cols(previous_knot_res.n_elem);
    for(int l=0;l<(int)previous_knot_res.n_elem;++l)
    {
      all_knot_cols(l)=l;
    }
    /*now we need set difference between all_knot_cols and all_parent_cols*/

    childless_knot_columns=my_setdiff(all_knot_cols,all_parent_cols);
    
    /*only does this calculation if there are possible knots to remove*/
    if(childless_knot_columns.n_elem>0)
    {
      arma::vec plusminusprobs(2);
      arma::vec plusminusrow;
      /*plus marginal holder just hodls 2 matricies*/

      plus_marginal_holder temp=minus_marginal_fun_parallel1(XX_old,sigmastar_old,mu,yy,y_ssq,a_g,previous_knot_res,a_pi,
                                                              b_pi,numthreads,childless_knot_columns,locations.n_cols,pi_method, g_method,r2_old,m0_size);

      minus_marginals=temp.plus_marginals;
        
      
      IntegerVector minusrownums=seq_len(minus_marginals.n_rows);
      arma::vec minusprobabilities=(exp(minus_marginals.col(1)-max(minus_marginals.col(1))));
      minus_row=dan_sample(minusrownums,minusprobabilities);
      minus_row(0)=minus_row(0)-1;
      
      plusminusprobs(0)=minus_marginals(minus_row(0),1);
      plusminusprobs(1)=plus_marginals(plus_row(0),knots.n_cols);
      plusminusprobs=(exp(plusminusprobs-max(plusminusprobs)));
      
      IntegerVector plusminusrownums=seq_len(2);
      add_knot=dan_sample(plusminusrownums ,plusminusprobs);
      add_knot(0)=add_knot(0)-1;
      
      /*now grab the top 100 minus models (usually won't be 100 of them)*/
      arma::vec minus_log_likelihoods=minus_marginals.col(1);
      int minus_topn=min(100,(int)minus_marginals.n_rows);
      /*   arma::uvec minus_top100_indicies;*/
      minus_top100_indicies=sort_index(minus_log_likelihoods,"descend");
      /*this now contains the indicies of the top 100 new plus models (or fewer)*/
      minus_top100_indicies=minus_top100_indicies.head(minus_topn);
      
      
      
      
      
      
      
      /*now loop through and assign the slice to be current knots with correct one removed*/
      
      
      
      
      for(int l=0;l<(int)minus_topn;++l)
      {
        arma::mat knotsminus(top_models_knots.n_rows+1,top_models_knots.n_cols);
        arma::vec resminus(top_model_knot_res.n_cols+1);
        resminus.fill(0);
        resminus.subvec(0,previous_knot_res.n_elem-1)=previous_knot_res;
        knotsminus.fill(arma::datum::nan);
        knotsminus.rows(0,knots.n_rows-1)=knots;
        
        
        knotsminus.shed_row(minus_marginals(minus_top100_indicies(l),0));
        resminus.shed_row(minus_marginals(minus_top100_indicies(l),0));
        top_models_knots.slice(200+l)=knotsminus;
        top_models_log_likelihood(200+l)=minus_log_likelihoods(minus_top100_indicies(l));
        /*GIANT SUBTLE CHANGE TO FIX EVERYTHING*/
        top_model_mus.row(200+l)=temp.mus.row(minus_top100_indicies(l));
        top_model_knot_res.row(200+l)=resminus.t();
      }
      
      
      
    }
    
    /*now for the updating, add_knot==1 if we're adding a knot, add_knot==0 if we're subtracting*/
    
    if(add_knot(0)==1)
    {
      
      /*update knot matrix*/
      arma::mat new_knots(knots.n_rows+1,knots.n_cols);
      new_knots.tail_rows(knots.n_rows)=knots;
      
      new_knots(0,0)=plus_marginals(plus_row(0),0);
      if(knots.n_cols==2)
      {
        new_knots(0,1)=plus_marginals(plus_row(0),1);
      }

      knots=new_knots;
      arma::sp_mat XX_new(yy.n_rows,1);
      /*create new column*/
      arma::sp_mat XX_new_holder(locations.n_rows,1);

      /*How to extract the resolution of the parent knot*/
      int new_resolution=previous_knot_res(parent_col(plus_row(0)))+1;
      
      XX_new_holder=create_column1(new_knots.row(0),new_resolution,locations,max_dist,nu,Kernel_type_c);
      
      XX_new.submat(0,0,locations.n_rows-1,0)=XX_new_holder;
      
     /* Rcpp::Rcout << XX_new.max()<<std::endl;*/
      /*update mathematical quantities*/
      double qstar11=accu(square(XX_new));
      arma::mat qstar1;
      /*following add 1 regression update, using notation and work from peter mueller wavelet notes*/
      qstar1=XX_old.t()*(XX_new);
      
 
      
      
      
      arma::mat hh=sigmastar_old*qstar1;
      arma::mat new_temp(1+sigmastar_old.n_rows,1+sigmastar_old.n_cols);
      new_temp.fill(0);
      new_temp(0,0)=1;

      new_temp.submat(1,1,new_temp.n_cols-1,new_temp.n_cols-1)=hh*hh.t();
      new_temp.submat(0,1,0,new_temp.n_cols-1)=-hh.t();
   
      
      new_temp.submat(1,0,new_temp.n_rows-1,0)=-hh;

      arma::mat sigmastar_new(sigmastar_old.n_rows+1,sigmastar_old.n_cols+1);
      sigmastar_new.fill(0);
 
      sigmastar_new.submat(1,1,sigmastar_new.n_cols-1,sigmastar_new.n_cols-1)=sigmastar_old;

      sigmastar_new=sigmastar_new+new_temp * (1/as_scalar((qstar11-(qstar1.t())*hh)));
      sigmastar_old=sigmastar_new;
      arma::SpMat<double> mu_possible(mu.n_rows+1,1);
      
      

      

      mu_possible.submat(1,0,mu.n_rows,0)=mu;

      XX_old=join_rows(arma::sp_mat(XX_new),XX_old);
      
      /*update resolutions */
      
      mu_possible=mu_possible+(1/as_scalar((qstar11-(qstar1.t())*hh)))*new_temp.col(0)*(new_temp.row(0)*XX_old.t()*yy);
      /*finish updating*/
      arma::vec temp(1);
      
      temp=new_resolution;
      previous_knot_res=join_cols(temp,previous_knot_res);
      mu=mu_possible;
   /*update parent_knot_columns*/
      for(int l=0;l<(int)parent_knot_columns.n_elem;++l)
      {
        parent_knot_columns(l)=parent_knot_columns(l)+1;
      }
      /*insert a new zero row into this*/
      parent_knot_columns.insert_rows(0,1);
      /*add the parent knot column to this spot*/
      parent_knot_columns(0)=parent_col(plus_row(0))+1;
      
    }
    else
    {
     XX_old.shed_col(childless_knot_columns(minus_row(0)));
      knots.shed_row(childless_knot_columns(minus_row(0)));
      previous_knot_res.shed_row(childless_knot_columns(minus_row(0)));
      /*update mathematical quantities*/
      double sigma_ll=sigmastar_old(childless_knot_columns(minus_row(0)),childless_knot_columns(minus_row(0)));
      arma::SpMat<double> sigma_l=sigmastar_old.col(childless_knot_columns((minus_row(0))));
      
      sigma_l.shed_row(childless_knot_columns(minus_row(0)));
      

      double old_muu=mu(childless_knot_columns(minus_row(0)),0);
      mu.shed_row(childless_knot_columns(minus_row(0)));
      mu=mu-sigma_l*old_muu/sigma_ll;
      sigmastar_old.shed_col(childless_knot_columns(minus_row(0)));
      sigmastar_old.shed_row(childless_knot_columns(minus_row(0)));
      sigmastar_old=sigmastar_old-(sigma_l)*sigma_l.t()/sigma_ll;
      /*update childless_knot_columns*/
      parent_knot_columns.shed_row(childless_knot_columns(minus_row(0)));
      /*update parent_knot_columns*/
      for(int l=0;l<(int)parent_knot_columns.n_elem;++l)
      {
        if(parent_knot_columns(l)>childless_knot_columns(minus_row(0)))
        {
          parent_knot_columns(l)=parent_knot_columns(l)-1;
        }
      }
    }
    /*now sort the log likelihoods largest to smallest, dedupe, and sort the cube the same way*/
    /*first remove duplicates*/
    /*Double precision makes this comparison too sensitive and will include two log likelihoods that are equal*/
    /*converting it to a float fixes things, maybe this is bad practice but works great*/
    arma::Col<float> temp_top= arma::conv_to<arma::fcolvec>::from(top_models_log_likelihood);
    arma::uvec unique_indicies=find_unique(temp_top);
    arma::vec unique_top_models_log_likelihood=top_models_log_likelihood(unique_indicies);
    arma::cube unique_top_model_knots=top_models_knots.slices(unique_indicies);
    arma::mat unique_top_model_res=top_model_knot_res.rows(unique_indicies);
    arma::mat unique_top_model_mus=top_model_mus.rows(unique_indicies);
    /*   now find indicies of top 100*/
    arma::uvec all_top100_indicies=sort_index(unique_top_models_log_likelihood,"descend");

    /*account for under 100 possible plus models*/
    int all_topn=min(100,(int)all_top100_indicies.n_elem);
    all_top100_indicies=all_top100_indicies.head(all_topn);
    
    /*now fill with a low value so they won't ever get picked*/
    top_models_log_likelihood.fill((min(plus_log_likelihoods)-1));
    top_models_knots.fill(arma::datum::nan);
    top_model_knot_res.fill(0);
    top_model_mus.fill(0);
    /*now fill with top 100 now that they'be been sorted*/
    top_models_log_likelihood.head(all_topn)=unique_top_models_log_likelihood(all_top100_indicies);
    
    top_models_knots.slices(0,all_topn-1)=unique_top_model_knots.slices(all_top100_indicies);
    top_model_knot_res.rows(0,all_topn-1)=unique_top_model_res.rows(all_top100_indicies);
    top_model_mus.rows(0,all_topn-1)=unique_top_model_mus.rows(all_top100_indicies);
    top_model_storage(i)=top_models_log_likelihood(0);
    top_model_size(i)=mu.n_rows;

    top_model_knot_res_iters.submat(0,i,previous_knot_res.n_elem-1,i)=previous_knot_res;
    double r2_new=1-accu(pow(yy-XX_old*mu,2))/y_ssq;
    current_model_r2(i)=r2_new;
    /*if the top model doesn't change for 10 iterations, then exit */
    if(i>=20)
    {
      if(top_model_storage(i)==top_model_storage(i-19))
      {
        Rcpp::Rcout << i << " iterations total" << std::endl;
        iterator=i+1;
        break;
      }
    }
    

    if((i+1)%10==0)
    {
      Rcpp::Rcout << 100*(double)(i+1)/ (double)maxiters << " percent complete " << std::endl;
      Rcpp::Rcout << i<< " iterations complete " << std::endl;
      Rcpp::Rcout << top_model_storage(i)<< " top log likelihood " << std::endl;
      Rcpp::Rcout << r2_new<< " current r2 " << std::endl;
      Rcpp::Rcout << XX_old.n_cols<< " current columns"<< std::endl;
    }


  }
  return Rcpp::List::create(Rcpp::Named("top_models_log_likelihood")=top_models_log_likelihood,
                            Rcpp::Named("top_models_knots")=top_models_knots,
                            Rcpp::Named("top_model_knot_res")=top_model_knot_res,
                            Rcpp::Named("top_model_mus")=top_model_mus,
							              Rcpp::Named("top_model_storage_iters")=top_model_storage,
							              Rcpp::Named("current_model_r2")=current_model_r2,
                            Rcpp::Named("top_model_size")= top_model_size,
                            Rcpp::Named("top_model_knot_res_iters")=top_model_knot_res,
                            Rcpp::Named("total_iters")=iterator
  )  ;
}


