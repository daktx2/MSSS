mvsss_interval_c_bma<-function(locations,results,design_mat=NULL,level,model_used=1:10,type="mean")
{
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
    design_matrix=create_column(knots_dataframe[1,1:results$params$spatial_dimension,drop=F],knot_resolutions[1],locations,maxdist,nu)
    for(i in 2:dim(knots_dataframe)[1])
    {
      design_matrix=cbind.spam(design_matrix,create_column(knots_dataframe[i,1:results$params$spatial_dimension,drop=F],knot_resolutions[i],locations,maxdist,nu))
    }
    
    if(is.null(design_mat)==F)
    {
      design_matrix=cbind.spam(design_matrix,as.spam(design_mat))
    }
    muu=results$cpp_run$top_model_mus[modelno,]
    muu=muu[1:(dim(design_matrix)[2])]
    pred_val[modelno,]=as.matrix(design_matrix%*%muu)
    
    #must reconstruct sigmahat, will get removed once i edit the C++ code
    design_matrix_data=create_column(knots_dataframe[1,1:results$params$spatial_dimension,drop=F],knot_resolutions[1],results$params$locations,maxdist,nu)
    for(i in 2:dim(knots_dataframe)[1])
    {
      design_matrix_data=cbind.spam(design_matrix_data,create_column(knots_dataframe[i,1:results$params$spatial_dimension,drop=F],knot_resolutions[i],results$params$locations,maxdist,nu))
    }
    if(is.null(design_mat)==F)
    {
      if(is.null(results$params$R1_prior))
      {
        design_matrix_data=cbind.spam(design_matrix_data,as.spam(results$params$design_mat))
      }
      else{
        design_mat_r1_prior=matrix(0,dim(results$params$R1_prior)[1],ncol(design_matrix_data))
        design_mat_r1_prior[,which(knot_resolutions==1)]=chol(results$params$R1_prior)
        design_matrix_data=rbind(design_matrix_data,design_mat_r1_prior)
        design_matrix_data=cbind.spam(design_matrix_data,as.spam(results$params$design_mat))
      }
      
    }
    pred_val_data=design_matrix_data%*%muu
    choler=chol(crossprod.spam(design_matrix_data))
    sigma_est=apply(forwardsolve.spam(choler,t(design_matrix))^2,2,sum)
    random.error=sum((results$params$yy_spam-pred_val_data)^2)/length(results$params$yy_spam)
    if(type=="pred") {SE_Preds=sqrt((sigma_est+1)*random.error)}
    if(type=="mean") {SE_Preds=sqrt((sigma_est)*random.error)}
    lower_bound[modelno,]=pred_val[modelno,]+qt(((1-level)/2),dim(design_matrix_data)[1]-dim(design_matrix_data)[2])*SE_Preds
    upper_bound[modelno,]=pred_val[modelno,]+qt(1-(1-level)/2,dim(design_matrix_data)[1]-dim(design_matrix_data)[2])*SE_Preds
  }
  logweights=results$cpp_run$top_models_log_likelihood[model_used]-max(results$cpp_run$top_models_log_likelihood[model_used])
  weightss=exp(logweights)/sum(exp(logweights))
  summary_pred=t(weightss)%*%pred_val
  summary_lower=t(weightss)%*%lower_bound
  summary_upper=t(weightss)%*%upper_bound
  return(list(preds=summary_pred,lower=summary_lower,upper=summary_upper))
}