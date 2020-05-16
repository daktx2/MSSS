#' Predictions and summaries from msss_fit
#' 
#' @param locations 1d vector or 2d matrix of locations of the locations you'd like to predict at
#' @param results output from mr_optim_fit
#' @param design_mat if fixed effects are part of the model, design matrix for the locations
#' @param model_used which of the models to use, less than or equal to results$total_iters
#' @param type 'pred' is default for prediction interval, 'mean' is interval for the mean, 'noint' is just a prediction without lower and upper bounds and is much faster, 'rescount' counts how many resoltuions are active at each point
#' @return list with either a row matrix of predictions or a row matrix of resolutions active
#' FILLER
mkl_fn_seq_em_tau_pred<-function(locations,results,design_mat=NULL,model_used=1,type="pred")
{
  nu=results$params$smoothness
  niu=results$params$kernel_width
  spatial_dimension=results$params$spatial_dimension
  if(spatial_dimension==2)
  {
    pred_val=matrix(0,1,dim(locations)[1])
  }
  if(spatial_dimension==1)
  {
    pred_val=matrix(0,1,length(locations))
  }
  counter=0
  for(modelno in model_used)
  {
    counter=counter+1
    #from 1 model set up knots
    #first set up data frame with resolution and location of each knot
    #extract knots and the resolutions for modelno
    knot_resolutions=results$results[[modelno]]$current_knots[,2]
    if(spatial_dimension==2)
    {
      knots_dataframe=results$results[[modelno]]$current_knots[,3:4]
      r1_knot_mindist=min(dist(knots_dataframe[which(knot_resolutions==1),]))
    }
    if(spatial_dimension==1)
    {
      knots_dataframe=results$results[[modelno]]$current_knots[,3]
      r1_knot_mindist=min(dist(knots_dataframe[which(knot_resolutions==1)]))
    }
    
    #create column fn
    
    maxdist=niu*r1_knot_mindist*((.5)^((0):99))
    if(spatial_dimension==2)
    {
      design_matrix=create_column(knots_dataframe[1,1:results$params$spatial_dimension,drop=F],knot_resolutions[1],locations,maxdist,nu)
      for(i in 2:dim(knots_dataframe)[1])
      {
        design_matrix=cbind.spam(design_matrix,create_column(knots_dataframe[i,1:results$params$spatial_dimension,drop=F],knot_resolutions[i],locations,maxdist,nu))
      }
    }
    if(spatial_dimension==1)
    {
      design_matrix=create_column(knots_dataframe[1],knot_resolutions[1],locations,maxdist,nu)
      for(i in 2:length(knots_dataframe))
      {
        design_matrix=cbind.spam(design_matrix,create_column(knots_dataframe[i],knot_resolutions[i],locations,maxdist,nu))
      }
    }
    
    if(is.null(design_mat)==F)
    {
      design_matrix=cbind.spam(design_matrix,as.spam(design_mat))
    }
    if(type=="rescount")
    {
      useful=design_matrix[]>0
      useful=useful*t(replicate(nrow(design_matrix),knot_resolutions))
      pred_val[counter,]=apply(useful,1,max)
    }
    else{
      muu=results$results[[modelno]]$fit[[1]]
      pred_val[counter,]=as.matrix(design_matrix%*%muu)
    }
    
  }
  return(pred_val)
}