#' Create Markov Random Field Precision matrix
#' 
#' @param knots_r1 1d vector or 2d matrix of knot locations.
#' @param precision parameter, positive constant.
#' @return nrow(knots_r1) square matrix with intrinsic MRF matrix.
#' @examples
#' Intrinsic_MRF_Prec(c(1,2,3),2)
Intrinsic_MRF_Prec<-function(knots_r1,tau)
{
  #create GMRF inverse covariance for resolution 1
  #nearest neighbor is the neighborhood definition for GMRF that we're using
  aaa=as.matrix(dist(knots_r1))
  diag(aaa)=max(aaa)
  whichmin=which(aaa==min(aaa))
  WW=matrix(0,dim(aaa)[1],dim(aaa)[2])
  WW[whichmin]=-1
  diag(WW)=-colSums(WW)+.01 #diagonal add noise
  return(WW*tau)
}
#' Create R1 knots on a grid a specificed number of knots wide
#' 
#' @param locations 1d vector or 2d matrix of locations of the data.
#' @param spatial_dimension 1 or 2
#' @param n_x how wide the grid should be in number of knots
#' @param buffer how much wider should the grid be than then data range
#' @return vector or n by 2 matrix of knot locations
#' @examples
#' Intrinsic_MRF_Prec(c(1,2,3),1,10,2)
r1_create<-function(locations,spatial_dimension,n_x=10,buffer=0)
{
  if(spatial_dimension==1)
  {
    knots_r1=seq(floor(min(locations)-buffer),ceiling(max(locations)+buffer),length=(n_x))
    return(knots_r1)
  }
  else if(spatial_dimension==2)
  {
    knots_r1_x=seq(floor(min(locations[,1])-buffer),ceiling(max(locations[,1]+buffer)),length=(n_x))
    domain_width_x=knots_r1_x[length(knots_r1_x)]-knots_r1_x[1]
    partition_width_x=domain_width_x/(length(knots_r1_x)-1)
    knots_r1_y=seq(floor(min(locations[,2])-buffer),
                   ceiling(max(locations[,2])+buffer+partition_width_x),by=partition_width_x)
    centerer=(max(knots_r1_y)-max(locations[,2]))/2
    knots_r1_y=knots_r1_y-centerer
    return(expand.grid(knots_r1_x,knots_r1_y))
  }
  else{
    warning("spatial dimension must be 1 or 2")
  }
}

#' Create R1 knots on a grid a specificed number of knots wide
#' 
#' @param locations locations of spatial observations.
#' @param yy observation at each location
#' @param knots_r1 locations of R1 knots, equally spaced grid, output of r1_create
#' @param spatial_dimension dimension of the locations, should be 1 or 2
#' @param maxiters number of iterations, default is 100
#' @param cores number of cores for parallel processing, ~20 is optimal
#' @param design_mat fixed effects, at least an intercept is recommended
#' @param stopping_rule 'maxiter' is default, could be 'log likelihood' if you want it to exit early when finding no new models
#' @param g_method 1 for hyper g, 1 i sdefault
#' @param a_pi parameter 1 for beta binomial sparsity if pi_method=1, or numerator if pi_method=2, 1 is default
#' @param b_pi parameter 2 for beta binomial sparsity if pi_method=1, or denominator if pi_method=2, 5 is default
#' @param a_g   hyper g prior parameter, 3 is default, 2-4 recommended
#' @param kernel_width bezier kernel width parameter, should be 1.5 or greater to be sensible, 1.5 is default
#' @param nu bezier smoothness parameter, default is 1
#' @param pi_method beta binomial is 1 (default), binomial (fixed pi) is 2
#' @param R1_prior covariance matrix for R1 knots, default is NULL
#' 
#' @return large list, exists to be processed by the function msss_predict forsummaries and predictions
#' @examples
#' blank
msss_fit<-function(locations,
                                           yy, #observation at each location
                                           knots_r1, #locations of R1 knots
                                           spatial_dimension, #dimension of the locations, should be 1 or 2
                                           maxiters=100, #number of iterations
                                           cores=1, #number of cores for parallel processing, ~20 is optimal
                                           design_mat=NULL, #fixed effects, at least
                                           #an intercept is recommended
                                           stopping_rule="maxiter", #this is 'maxiter' if model should always run maxiter times
                                           #could be 'log likelihood' if you want the algorithm to stop once it's
                                           #not exploring new models (stops after 20 iterations with no new best model)
                                           g_method=1,              #1 for hyper g
                                           a_pi=1,                  #parameter 1 for beta binomial sparsity if pi_method=1, or numerator if pi_method=2
                                           b_pi=5,                  #parameter 2 for beta binomial sparsity if pi_method=1, or denominator if pi_method=2
                                           a_g=3,                   #if hyper g, we have a parameter here, means different in Bayarri
                                           kernel_width=1.5,        #bezier kernel width parameter, should be 1.5 or greater to be sensible
                                           nu=1,                    #bezier smoothness parameter
                                           pi_method=1,              #beta binomial is 1
                                           R1_prior=NULL            #covariance matrix for R1 knots
)      
{
  start=proc.time()
  #bezier kernel for initial knots_r1
  niu=kernel_width
  
  #number of possible kernels
  p=2^spatial_dimension
  #maximum distance for compact kernels at each resolution
  max_dist=niu*min(dist(knots_r1))*((.5)^((0):99))
  #how far the kernels are from each other
  resolution_dist=min(dist(knots_r1))*((.5)^((0):99))
  
  #form design matrix for R1
  K_dists_r1=rdist(locations,knots_r1)
  K_kernel_r1=(1-(K_dists_r1/max_dist[1])^2)^nu
  K_kernel_r1[which(K_dists_r1>max_dist[1])]<-0
  #if there is no R1 prior then we need to exclude empty columns
  #we also exclude near empty columns for stability
  if(is.null(R1_prior))
  {
    #bad knots_r1
    #bad knots are where there are fewer than 5 datapoints in the kernel
    bad_knot_index=which(apply(K_kernel_r1!=0,2,sum)<5)
    if(spatial_dimension==2){
      bad_knots=knots_r1[bad_knot_index,]
    }
    if(spatial_dimension==1){
      bad_knots=knots_r1[bad_knot_index]
    }
    if(length(bad_knot_index)>0&spatial_dimension==2)
    {
      knots_r1=knots_r1[-bad_knot_index,]
      XX=K_kernel_r1[,-bad_knot_index]
    }
    if(length(bad_knot_index)>0&spatial_dimension==1)
    {
      knots_r1=knots_r1[-bad_knot_index]
      XX=K_kernel_r1[,-bad_knot_index]
    }
    if(length(bad_knot_index)==0)
    {
      XX=K_kernel_r1
    }
    
    XX_spam=as.spam(XX)
  }else{
    if(dim(R1_prior)[1]!=dim(K_kernel_r1)[2])
    {stop("Prior is not same dimension as J(1)")}
    XX=K_kernel_r1
    XX_spam=as.spam(XX)
    choler=chol(R1_prior)
    yy=c(yy,rep(0,dim(R1_prior)[1]))
    XX=rbind(XX,choler)
    XX_spam=rbind.spam(XX_spam,choler)
    if(is.null(design_mat)==F)
    {
      design_mat_filler=matrix(0,dim(R1_prior)[1],dim(design_mat)[2])
      design_mat=rbind(design_mat,design_mat_filler)
    }
  }

  
  if(is.null(design_mat)==F)
  {
    XX=cbind(XX,design_mat)
    XX_spam=cbind.spam(XX_spam,as.spam(design_mat))
  }
  
  #start with just r1 knots, fit linear model, extract design matrix and coefficients/cov mat
  lm_temp=lm(as.numeric(as.matrix(yy))~XX+0)
  sigmastar_old=as.spam(vcov(lm_temp)/(summary(lm_temp)$sigma^2))
  mu=as.spam(coef(lm_temp))
  XX_old=XX_spam
  y_ssq=sum((yy-mean(as.numeric(as.matrix(yy))))^2)
  r2_old=1-sum((yy-XX_old%*%mu)^2)/y_ssq
  
  
  n=dim(XX_spam)[1]
  
  if(spatial_dimension==2){
    knotID=apply(as.matrix(cbind(1,knots_r1$Var1,knots_r1$Var2)),1,paste0,collapse="_")
  }
  if(spatial_dimension==1){
    knotID=apply(as.matrix(cbind(1,knots_r1)),1,paste0,collapse="_")
  }
  if(spatial_dimension==2)
  {
    resolution=rep(1,nrow(knots_r1))
    current_knots=data.frame(knotID=knotID,
                             resolution=rep(1,nrow(knots_r1)),
                             xloc=knots_r1$Var1,
                             yloc=knots_r1$Var2,
                             child1ID=rep(NA,nrow(knots_r1)),
                             child2ID=rep(NA,nrow(knots_r1)),
                             child3ID=rep(NA,nrow(knots_r1)),
                             child4ID=rep(NA,nrow(knots_r1)),
                             XX_column=1:nrow(knots_r1),
                             stringsAsFactors=FALSE)
  }
  if(spatial_dimension==1)
  {
    resolution=rep(1,length(knots_r1))
    current_knots=data.frame(knotID=knotID,
                             resolution=rep(1,length(knots_r1)),
                             xloc=knots_r1,
                             child1ID=rep(NA,length(knots_r1)),
                             child2ID=rep(NA,length(knots_r1)),
                             child3ID=rep(NA,length(knots_r1)),
                             child4ID=rep(NA,length(knots_r1)),
                             XX_column=1:length(knots_r1),
                             stringsAsFactors=FALSE)
  }
  
  class(current_knots$child1ID)<-class(current_knots$child2ID)<-class(current_knots$child3ID)<-class(current_knots$child4ID)<-"character"
  
  #cpp to r convert pi method
  if(pi_method==1){pi_method_r="beta binomial pi"}
  if(pi_method==2){pi_method_r="fixed pi"}
  if(pi_method==3){pi_method_r="giantdata"}
  #calculate marginal for this initialmodel
  if(pi_method_r=="fixed pi")
  {
    params=list(resolution,a_pi/b_pi,spatial_dimension)
  }
  if(pi_method_r=="beta binomial pi")
  {
    params=list(resolution,a_pi,b_pi,spatial_dimension)
  }
  if(pi_method_r=="giantdata")
  {
    params=list(resolution,a_pi,spatial_dimension)
  }
  if(g_method==1)
  {
    log_marginal_lik=0 #hyper g with Bayarri and Berger correction, base model gets 1 (since its marginal is part of the marginal in the rest)
  }
  #this will be general bayarri and berger but still is a zero here
  if(g_method==2)
  {
    log_marginal_lik=(log(a_g-2)-log(dim(XX)[2]+a_g-2))+hypergeometric2F1((n-1)/2,1,(dim(XX)[2]+a_g)/2, r2_old, method = "Laplace", log = TRUE)+calculate_prior(params,pi_method_r)
  }
  
  params=list(locations=locations,yy=yy,knots_r1=knots_r1,maxiters=maxiters,spatial_dimension=spatial_dimension,
              stopping_rule=stopping_rule,pi_method=pi_method,g_method=g_method,a_pi=a_pi,b_pi=b_pi,
              a_g=a_g,kernel_width=kernel_width,nu=nu,design_mat=design_mat,R1_prior=R1_prior)
  #here we will do the data augmentation for the proper prior on the first resolution
  
  #now comes the iterative part
  
  #note that I've added r2_old to the c++ call, this is the base model's r squared value, which is needed
  #in the computation of correct Bayes factors for models with more than an intercept
  if(spatial_dimension==2){
    cpp_run=.full_looper(matrix(locations,dim(locations)[1],dim(locations)[2]),
                        as.matrix(yy),
                        unname(as.matrix(knots_r1)),
                        resolution,
                        maxiters,
                        as.dgRMatrix.spam(XX_old),
                        as.dgRMatrix.spam(sigmastar_old),
                        as.dgRMatrix.spam(mu),
                        a_pi,
                        b_pi,
                        a_g,
                        nu,
                        max_dist,
                        resolution_dist,
                        log_marginal_lik,
                        y_ssq,
                        cores,
                        pi_method,
                        g_method,
                        r2_old,
                        m0_size=dim(XX_old)[2]
    )
  }
  
  if(spatial_dimension==1)
  {
    cpp_run=.full_looper(matrix(locations),
                        as.matrix(yy),
                        unname(as.matrix(knots_r1)),
                        resolution,
                        maxiters,
                        as.dgRMatrix.spam(XX_old),
                        as.dgRMatrix.spam(sigmastar_old),
                        as.dgRMatrix.spam(mu),
                        a_pi,
                        b_pi,
                        a_g,
                        nu,
                        max_dist,
                        resolution_dist,
                        log_marginal_lik,
                        y_ssq,
                        cores,
                        pi_method,
                        g_method,
                        r2_old,
                        m0_size=dim(XX_old)[2]
    )
  }
  runtime=proc.time()-start
  returner=list(cpp_run=cpp_run,params=params,runtime=runtime)
  return(returner)
}