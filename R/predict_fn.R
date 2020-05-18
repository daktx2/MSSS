.create_column<-function (knot, resolution, locations, max_dist, nu, Kernel_type, spatial_dimension) 
{
  temp = fields::rdist(locations, knot)
  if(Kernel_type=="bezier")
  {
    K_kernel_temp = (1 - (temp/max_dist[resolution])^2)^nu
    K_kernel_temp[is.nan(K_kernel_temp)] <- 0
    K_kernel_temp[which(temp > max_dist[resolution])] <- 0
    K_kernel_temp = spam::as.spam(K_kernel_temp)
  }
  if(Kernel_type=="wendland")
  {
    ll=floor(spatial_dimension/2)+2
    K_kernel_temp=(pmax(1-(temp/max_dist[resolution]),0))^{ll+1}
    K_kernel_temp=K_kernel_temp*(1+(ll+1)*(temp/max_dist[resolution]))
  }
  return(K_kernel_temp)
}
#' Predictions and summaries from msss_fit
#' 
#' @param locations 1d vector or 2d matrix of locations of the locations you'd like to predict at
#' @param results output from msss_fit
#' @param design_mat if fixed effects are part of the model, design matrix for the locations
#' @param level confidence level, not always needed, .95 is defaykt
#' @param model_used vector of how many top models to use for Bayesian model averaging, 1:100 is default, if slow try 1:10
#' @param type 'pred' is default for prediction interval, 'mean' is interval for the mean, 'noint' is just a prediction without lower and upper bounds and is much faster, 'rescount' counts how many resoltuions are active at each point
#' @return list with either a row matrix of predictions, a row matrix of resolutions active, or 3 row matricies with upper bounds, predicted values, and lower bound
#' @examples
#' FILLER
msss_pred<-function(locations,results,design_mat=NULL,level=.95,model_used=1:100,type="pred")
{
  Kernel_type=results$params$Kernel_type
  spatial_dimension=results$params$spatial_dimension
  nu=results$params$nu
  niu=results$params$kernel_width
  spatial_dimension=results$params$spatial_dimension
  pred_val=matrix(0,length(model_used),dim(locations)[1])
  upper_bound=matrix(0,length(model_used),dim(locations)[1])
  lower_bound=matrix(0,length(model_used),dim(locations)[1])
  for(modelno in model_used)
  {
    
    #from 1 model set up knots
    #first set up data frame with resolution and location of each knot
    #extract knots and the resolutions for modelno
    knots_dataframe=data.frame(results$cpp_run$top_models_knots[,,modelno])
    knots_dataframe=knots_dataframe[complete.cases(knots_dataframe),]
    knots_dataframe=as.matrix(knots_dataframe)
    knot_resolutions=results$cpp_run$top_model_knot_res[modelno,]
    knot_resolutions=knot_resolutions[which(knot_resolutions!=0)]
    knot_resolutions=(knot_resolutions)
    #create column fn
    
    r1_knot_mindist=min(dist(knots_dataframe[which(knot_resolutions==1),]))
    maxdist=niu*r1_knot_mindist*((.5)^((0):99))
    design_matrix=.create_column(knots_dataframe[1,1:results$params$spatial_dimension,drop=F],knot_resolutions[1],locations,maxdist,nu, Kernel_type, spatial_dimension)
    for(i in 2:dim(knots_dataframe)[1])
    {
      design_matrix=spam::cbind.spam(design_matrix,.create_column(knots_dataframe[i,1:results$params$spatial_dimension,drop=F],knot_resolutions[i],locations,maxdist,nu, Kernel_type, spatial_dimension))
    }
    #for resolution counting
    if(type=="rescount")
    {
      useful=design_matrix[]>0
    }

    if(is.null(design_mat)==F)
    {
      design_matrix=spam::cbind.spam(design_matrix,spam::as.spam(design_mat))
    }
    muu=results$cpp_run$top_model_mus[modelno,]
    muu=muu[1:(dim(design_matrix)[2])]
    pred_val[modelno,]=as.matrix(design_matrix%*%muu)
    if(!(type %in%(c("rescount","noint"))))
    {
      #must reconstruct sigmahat for intervals
      design_matrix_data=.create_column(knots_dataframe[1,1:results$params$spatial_dimension,drop=F],knot_resolutions[1],results$params$locations,maxdist,nu, Kernel_type, spatial_dimension)
      for(i in 2:dim(knots_dataframe)[1])
      {
        design_matrix_data=spam::cbind.spam(design_matrix_data,.create_column(knots_dataframe[i,1:results$params$spatial_dimension,drop=F],knot_resolutions[i],results$params$locations,maxdist,nu, Kernel_type, spatial_dimension))
      }
      if(is.null(design_mat)==F)
      {
        if(is.null(results$params$R1_prior))
        {
          design_matrix_data=spam::cbind.spam(design_matrix_data,spam::as.spam(results$params$design_mat))
        }
        else{
          design_mat_r1_prior=matrix(0,dim(results$params$R1_prior)[1],ncol(design_matrix_data))
          design_mat_r1_prior[,which(knot_resolutions==1)]=chol(results$params$R1_prior)
          design_matrix_data=rbind(design_matrix_data,design_mat_r1_prior)
          design_matrix_data=spam::cbind.spam(design_matrix_data,spam::as.spam(results$params$design_mat))
        }
        
      }
      pred_val_data=design_matrix_data%*%muu
      choler=chol(spam::crossprod.spam(design_matrix_data))
      sigma_est=apply(spam::forwardsolve.spam(choler,t(design_matrix))^2,2,sum)
      random.error=sum((results$params$yy-pred_val_data)^2)/length(results$params$yy)
    }
    
    if(type=="pred") {SE_Preds=sqrt((sigma_est+1)*random.error)}
    if(type=="mean") {SE_Preds=sqrt((sigma_est)*random.error)}
    if(type=="rescount")
    {
      useful=useful*t(replicate(nrow(useful),knot_resolutions))
      pred_val[modelno,]=apply(useful,1,max)
    }
    if(!(type%in%c("rescount","noint")))
    {
      lower_bound[modelno,]=pred_val[modelno,]+qt(((1-level)/2),dim(design_matrix_data)[1]-dim(design_matrix_data)[2])*SE_Preds
      upper_bound[modelno,]=pred_val[modelno,]+qt(1-(1-level)/2,dim(design_matrix_data)[1]-dim(design_matrix_data)[2])*SE_Preds
    }
  }
  logweights=results$cpp_run$top_models_log_likelihood[model_used]-max(results$cpp_run$top_models_log_likelihood[model_used])
  weightss=exp(logweights)/sum(exp(logweights))
  summary_pred=t(weightss)%*%pred_val
  if(type=="rescount")
  {
    return(list(preds=summary_pred))
  }
  if(type=="noint")
  {
    return(list(preds=summary_pred))
  }
  summary_lower=t(weightss)%*%lower_bound
  summary_upper=t(weightss)%*%upper_bound
  return(list(preds=summary_pred,lower=summary_lower,upper=summary_upper))
}