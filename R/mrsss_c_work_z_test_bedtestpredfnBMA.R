msss_pred<-function(locations,results,design_mat=NULL,model_used=1:100,type="standard")
{
  nu=results$params$nu
  niu=results$params$kernel_width
  spatial_dimension=results$params$spatial_dimension
  pred_val=matrix(0,length(model_used),dim(locations)[1])

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


    if(type=="standard")
    { 
      if(is.null(design_mat)==F)
      {
        design_matrix=cbind.spam(design_matrix,as.spam(design_mat))
      }
      muu=results$cpp_run$top_model_mus[modelno,]
      muu=muu[1:(dim(design_matrix)[2])]
      pred_val[modelno,]=as.matrix(design_matrix%*%muu)
    }
    if(type=="rescount")
    {
      useful=design_matrix[]>0
      useful=useful*t(replicate(nrow(design_matrix),knot_resolutions))
      pred_val[modelno,]=apply(useful,1,max)
    }
  }
  logweights=results$cpp_run$top_models_log_likelihood[model_used]-max(results$cpp_run$top_models_log_likelihood[model_used])
  weightss=exp(logweights)/sum(exp(logweights))
  summary_pred=t(weightss)%*%pred_val
  return(summary_pred)
}