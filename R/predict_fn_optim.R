#' Predictions and summaries from mr_optim_fit
#' 
#' @param locations 1d vector or 2d matrix of locations of the locations you'd like to predict at
#' @param results output from mr_optim_fit
#' @param design_mat if fixed effects are part of the model, design matrix for the locations
#' @param marg_type howt to calculate the marginal, "bic_CAP" is default and strongly recommended but "hyperg_refit" is an option
#' @param int_type "pred" for prediction, "recount" to count resolutions, "noint" for a prediction only
#' @param level confidence level for the interval
#' @param a_beta a parameter for beta prior on model space, 1 is default
#' @param b_beta b parameter for beta prior on model space, 1 is default, but if data is large, use 10^(3*n/10^4)
#' @param a_g if "hyperg_refit" is used
#' @return list with either a row matrix of predictions or a row matrix of resolutions active
#' FILLER
mr_optim_pred<-function(locations, results, design_mat = NULL, marg_type="bic_CAP", int_type = "pred",level=.95,a_beta=1,b_beta=1,a_g=3)
{
  #locate MAP estimates
  lambda_tracker=rep(0,length(results$results))
  for(i in 1:length(lambda_tracker))
  {
    lambda_tracker[i]=results$results[[i]]$lambda1
  }
  
  
  
  
  #find converged indices
  converged=c(which(diff(lambda_tracker)!=0),length(results$results))
  
  
  
  nu = results$params$kernel_width
  spatial_dimension = results$params$spatial_dimension
  
  #first compute BMA weights
  
  log_BMA_weights=rep(0,length(converged))
  log_prior=rep(0,length(converged))
  log_marginal=rep(0,length(converged))
  #extract some useful quantities
  yy=as.numeric(as.matrix(results$params$yy))
  total_ssq=sum((yy - mean(as.numeric(as.matrix(yy))))^2)
  p0=sum(results$results[[1]]$fit[[1]]!=0)
  r2_m0= 1 - sum((yy - results$results[[1]]$current_design %*% (results$results[[1]]$fit[[1]]))^2)/total_ssq
  counter=0
  
  store_muu=list()
  store_design=list()
  
  for(i in converged)
  {
    counter=counter+1
    resolutions=results$results[[i]]$current_knots[which(results$results[[i]]$fit[[1]]!=0),2]
    params=list(resolutions,a_beta,b_beta,results$params$spatial_dimension)
    log_prior[counter]=calculate_prior_beta_pi(params)
    
    #extract beta
    muu=results$results[[i]]$fit[[1]][which(results$results[[i]]$fit[[1]]!=0)]
    #extract design mat
    design_mat_extract=results$results[[i]]$current_design[,which(results$results[[i]]$fit[[1]]!=0)]
    store_design[[counter]]=design_mat_extract
    if(marg_type=="bic_CAP")
    {
      log_marginal[counter]=-(log(length(yy))*length(muu)+.5*length(yy)*log(sum((yy - results$results[[i]]$current_design %*% (results$results[[i]]$fit[[1]]))^2)/length(yy)))
    }
    if(marg_type=="hyperg_refit")
    {
      #need to reconstruct GMRF for conjugate gradient descent
      aaa=as.matrix(dist(results$params$knots_r1))
      diag(aaa)=max(aaa)
      whichmin=which(aaa==min(aaa))
      WW=matrix(0,dim(aaa)[1],dim(aaa)[2])
      WW[whichmin]=-1
      diag(WW)=-colSums(WW)+.01 #diagonal add noise
      GMRF=WW*results$results[[i]]$tau
      prior=bdiag(GMRF,diag(rep(0,length(muu)-ncol(GMRF))))
      
      refit_muu=multires.conjGrad(as.matrix(prior)+as.matrix(crossprod(design_mat_extract)),c(as.matrix(t(design_mat_extract)%*%yy)),muu)
      store_muu[[counter]]=refit_muu
      pp=length(refit_muu)
      r2_new=1-sum((yy-design_mat_extract%*%refit_muu)^2)/total_ssq;
      log_marginal[counter]=log(2*a_g)-log(pp-p0+2*a_g)+BAS::hypergeometric2F1((length(yy)-p0)/2,1,(a_g+1+pp-p0)/2,1-(1-r2_new)/(1-r2_m0),method="Laplace")
    }
    log_BMA_weights[counter]=log_marginal[counter]+log_prior[counter]
  }
  
  
  returner=list()
  
  
  
  counter = 0
  
  for (modelno in converged) 
  {
    iterator=list()
    counter = counter + 1
    muu = results$results[[modelno]]$fit[[1]]
    good_knot_index=which(muu!=0)
    muu = muu[good_knot_index]
    
    knot_resolutions = results$results[[modelno]]$current_knots[,2]
    if (spatial_dimension == 2) 
    {
      knots_dataframe = results$results[[modelno]]$current_knots[good_knot_index+length(design_mat[,1]),3:4]
      r1_knot_mindist = min(dist(knots_dataframe[which(knot_resolutions ==  1), ]))
    }
    if (spatial_dimension == 1) 
    {
      knots_dataframe = results$results[[modelno]]$current_knots[,3]
      r1_knot_mindist = min(dist(knots_dataframe[which(knot_resolutions == 1)]))
    }
    maxdist = nu * r1_knot_mindist * ((0.5)^((0):99))
    if (spatial_dimension == 2) 
    {
      design_matrix = create_column(knots_dataframe[1,1:results$params$spatial_dimension, drop = F],knot_resolutions[1], locations, maxdist, nu)
      for(i in 2:nrow(knots_dataframe))
      {
        design_matrix = cbind.spam(design_matrix, create_column(knots_dataframe[i,1:results$params$spatial_dimension, drop = F],knot_resolutions[i], locations, maxdist, nu))
      }
    }
    if (spatial_dimension == 1) 
    {
      design_matrix = create_column(knots_dataframe[1],knot_resolutions[1], locations, maxdist, nu)
      for (i in 2:length(knots_dataframe)) 
      {
        design_matrix = cbind.spam(design_matrix, create_column(knots_dataframe[i],knot_resolutions[i], locations, maxdist, nu))
      }
    }
    if (is.null(design_mat) == F) {
      design_matrix = cbind.spam(design_matrix, as.spam(design_mat))
    }
    
    
    if (int_type == "rescount") {
      useful = design_matrix[] > 0
      useful = useful * t(replicate(nrow(design_matrix), 
                                    knot_resolutions))
      preds = apply(useful, 1, max)
      iterator$preds=preds
      
    }
    else 
    {
      if(marg_type=='bic_CAP')
      {
        preds = as.matrix(design_matrix %*%muu)
        iterator$preds=preds
      }
      if(marg_type=='hyperg_refit')
      {
        preds = as.matrix(design_matrix %*%store_muu[[counter]])
        iterator$preds=preds
      }
      if(int_type!="noint")
      {
        
        if(marg_type=="bic_CAP")
        {
          pred_val_data=store_design[[counter]]%*%muu
          choler=chol(spam::crossprod.spam(store_design[[counter]]))
          sigma_est=apply(forwardsolve(choler,t(design_matrix))^2,2,sum)
          random.error=sum((results$params$yy-pred_val_data)^2)/length(results$params$yy)
        }
        
        if(marg_type=="hyperg_refit")
        {
          pred_val_data=store_design[[counter]]%*%store_muu[[counter]]
          choler=chol(spam::crossprod.spam(store_design[[counter]]))
          sigma_est=apply(forwardsolve(choler,t(design_matrix))^2,2,sum)
          random.error=sum((results$params$yy-pred_val_data)^2)/length(results$params$yy)
        }
        
        
        
        SE_Preds=sqrt((sigma_est+1)*random.error)
        lower=preds+qt(((1-level)/2),dim(store_design[[counter]])[1]-dim(store_design[[counter]])[2])*SE_Preds
        upper=preds+qt(1-(1-level)/2,dim(store_design[[counter]])[1]-dim(store_design[[counter]])[2])*SE_Preds
        iterator$upper=upper
        iterator$lower=lower
      }
    }
    returner[[counter]]=iterator 
  }
  
  return(list(returner=returner,log_BMA_weights=log_BMA_weights,log_marginal=log_marginal,log_prior=log_prior))
}