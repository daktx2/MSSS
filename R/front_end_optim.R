#' Runs multiresolution optimization routine
#' 
#' @param yy response values
#' @param locationss locations of observations, either 1d vector or matrix with 2 columns
#' @param spatial_dimension spatial dimension of observations, either 1 or 2
#' @param kernel_width width of the kernel
#' @param smoothness smoothness parameter for kernel
#' @param shrinkage vector of parameter d, the shrinkage that increases in resolution, long vector
#' @param maxiter maximum number of iterations to run
#' @param lambdatree_seq decreasing sequence of lambdas, the parameter that controls sparsity
#' @param tau_init initial tau value
#' @param tau_a a parameter for gamma prior on tau
#' @param tau_b b parameter for gamma prior on tau
#' @param em_tol tolerance for when tau has converged
#' @param m_tol default is .00001, tolerance for when the coefficients have converged
#' @param design_mat optional design matrix for fixed effects
#' @return list to be processed by prediction_function
mr_optim_fit<-function(yy,locationss,knots_r1,spatial_dimension,kernel_width,smoothness,shrinkage,maxiter,lambdatree_seq,tau_init=tau_init,
                            tau_a=tau_a,
                            tau_b=tau_b, em_tol=.001,m_tol=1e-5, design_mat=NULL)
{
  params=list(yy=yy,
              locationss=locationss,
              knots_r1=knots_r1,
              spatial_dimension=spatial_dimension,
              kernel_width=kernel_width,
              smoothness=smoothness,
              shrinkage=shrinkage,
              maxiter=maxiter,
              lambdatree_seq=lambdatree_seq,
              tau_init=tau_init,
              tau_a=tau_a,
              tau_b=tau_b,
              design_mat=design_mat)
  if(is.null(design_mat))
  {
    design_cols=0
  }
  else{
    design_cols=ncol(design_mat)
  }
  init_calc="GMRF"
  tree_it=1
  times=0
  tau=tau_init
  lambda1=lambdatree_seq[tree_it]
  max_dist=min(dist(knots_r1))*(kernel_width)*(.5^(0:100))
  nu=smoothness
  #this will contain the knot locations, parent index, and XX index
  #initialize first resolution
  if(spatial_dimension==1)
  {
    current_knots=data.frame(knot_row=1:length(knots_r1),
                             resolution=rep(1,length(knots_r1)),
                             xloc=knots_r1,
                             yloc=NA,
                             child1row=NA,
                             child2row=NA,
                             child3row=NA,
                             child4row=NA,
                             XX_column=1:length(knots_r1),
                             stringsAsFactors=FALSE)
  }
  if(spatial_dimension==2)
  {
    current_knots=data.frame(knot_row=1:nrow(knots_r1),
                             resolution=rep(1,nrow(knots_r1)),
                             xloc=knots_r1[,1],
                             yloc=knots_r1[,2],
                             child1row=NA,
                             child2row=NA,
                             child3row=NA,
                             child4row=NA,
                             XX_column=1:nrow(knots_r1),
                             stringsAsFactors=FALSE)
  }
  
  current_knots=as.matrix(current_knots)
  #now create and fill in 2nd res
  prevknotct=nrow(current_knots)
  for(i in 1:prevknotct)
  {
    current_loc=which(current_knots[,1]==i)
    current_knots[current_loc,5]=nrow(current_knots)+1
    current_knots[current_loc,6]=nrow(current_knots)+2
    #EDIT FOR 2D
    if(spatial_dimension==1)
    {
      aa=.create_knots_1d(current_knots[current_loc,3],current_knots[current_loc,2],2,min(dist(knots_r1))*.5^(0:10))
      newrow=c(NA,aa[1,2],aa[1,1],NA,NA,NA,NA,NA,NA)
      current_knots=rbind(current_knots,newrow)
      newrow=c(NA,aa[2,2],aa[2,1],NA,NA,NA,NA,NA,NA)
      current_knots=rbind(current_knots,newrow)
    }
    if(spatial_dimension==2)
    {
      current_knots[current_loc,7]=nrow(current_knots)+3
      current_knots[current_loc,8]=nrow(current_knots)+4
      aa=.create_knots_2d(matrix(current_knots[current_loc,3:4],1,2),current_knots[current_loc,2],2*spatial_dimension,min(dist(knots_r1))*.5^(0:10))
      newrow=c(NA,aa[1,3],aa[1,1],aa[1,2],NA,NA,NA,NA,NA)
      current_knots=rbind(current_knots,newrow)
      newrow=c(NA,aa[2,3],aa[2,1],aa[2,2],NA,NA,NA,NA,NA)
      current_knots=rbind(current_knots,newrow)
      newrow=c(NA,aa[3,3],aa[3,1],aa[3,2],NA,NA,NA,NA,NA)
      current_knots=rbind(current_knots,newrow)
      newrow=c(NA,aa[4,3],aa[4,1],aa[4,2],NA,NA,NA,NA,NA)
      current_knots=rbind(current_knots,newrow)
    }
    
  }
  #now fix the row assignments
  current_knots[,1]<-current_knots[,9]<-1:nrow(current_knots)
  
  
  #now we create the design matrix
  #kernel parameters
  
  #might want to rewrite this into C at some point or at least parallelize
  lister=list()
  for(i in 1:nrow(current_knots))
  {
    if(spatial_dimension==1)
    {
      knott=current_knots[i,3]
    }
    if(spatial_dimension==2)
    {
      knott=matrix(current_knots[i,3:4],1,2)
    }
    ress=current_knots[i,2]
    #this can be improved
    lister[[i]]=as(spam::as.dgCMatrix.spam(.create_column(knott,ress,locationss,max_dist,nu,'bezier')),"sparseVector")
    #print(i)
  }
  #function sourced from stackexchange for cbinding a bunch of sparse vectors
  sv.cbind <- function (...) {
    input <- lapply( ..., as, "dsparseVector" )
    thelength <- unique(sapply(input,length))
    stopifnot( length(thelength)==1 )
    return( Matrix::sparseMatrix( 
      x=unlist(lapply(input,slot,"x")), 
      i=unlist(lapply(input,slot,"i")), 
      p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
      dims=c(thelength,length(input))
    ) )
  }
  KK_mat=sv.cbind(lister)
  
  resolutions_non1=current_knots[which(current_knots[,2]!=1),2]
  n_groups=1+sum(current_knots[,2]!=1)
  n_r1=sum(current_knots[,2]==1)
  eta_g=c(10^-100,shrinkage[resolutions_non1])
  
  groups=matrix(0,n_groups,n_groups)
  tracker=1
  for(i in 1:nrow(current_knots))
  {
    if(current_knots[i,2] > 1)
    {
      tracker=tracker+1
    }
    children=current_knots[i,5:8]
    if(sum(is.na(children))==4)
    {
      next
    }
    children=children[!is.na(children)]
    groups[children-n_r1+1,tracker]=1
  }
  
  groups=Matrix::Matrix(groups)
  groups=groups>0
  #remember a minus 1 is needed since C is 0 indexed
  own_variables=as.integer(c(1,current_knots[which(current_knots[,2]>1),9]+design_cols)-1)
  if(spatial_dimension==1)
  {
    N_own_variables=as.integer(c(length(knots_r1)+design_cols,rep(1,n_groups-1)))
    num_r1=length(knots_r1)
  }
  if(spatial_dimension==2)
  {
    N_own_variables=as.integer(c(nrow(knots_r1)+design_cols,rep(1,n_groups-1)))
    num_r1=nrow(knots_r1)
  }
  
  tree=list(eta_g=eta_g,
            groups=groups,
            own_variables=own_variables,
            N_own_variables=N_own_variables)
  if(init_calc=="LMR1")
  {
    if(spatial_dimension==1)
    {
      lmm=MatrixModels:::lm.fit.sparse(KK_mat[,1:length(knots_r1)],c(yy))
    }
    if(spatial_dimension==2)
    {
      lmm=MatrixModels:::lm.fit.sparse(KK_mat[,1:(dim(knots_r1)[1])],c(yy))
    }
    init_val=matrix(c(lmm,rep(0,ncol(KK_mat)-length(lmm))))
  }
  if(init_calc=="ridge")
  {
    init_val=as.matrix(solve(Diagonal(dim(KK_mat)[2])+Matrix::t(KK_mat)%*%KK_mat)%*%Matrix::t(KK_mat)%*%c(yy))
  }
  if(init_calc=="lasso")
  {
    init_val=as.matrix((glmnet(KK_mat,yy,family="gaussian", alpha=1,lambda=lambda1,intercept = F,penalty.factor = c(rep(0,num_r1),rep(1,(dim(KK_mat)[2]-num_r1)))))$beta)
  }
  #init_value_for_sigma
  #fitt=spams.fistaTree((matrix(c(yy))),KK_mat,init_val,tree,TRUE,numThreads = -1,verbose = FALSE,
  #                     lambda1 =lambda1, it0 = 10, max_it = 100,
  #                     L0 = 0.1, tol = 1e-5, intercept = FALSE,
  #                     pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-l2')
  #fitt needs to be updated to fit the GMRF in resolution 1
  #update cholesky/eigenval/eigenvec to have block form with zeros for non R1
  
  #create GMRF inverse covariance for resolution 1
  #nearest neighbor is the neighborhood definition for GMRF that we're using
  aaa=as.matrix(dist(knots_r1))
  diag(aaa)=max(aaa)
  whichmin=which(aaa==min(aaa))
  WW=matrix(0,dim(aaa)[1],dim(aaa)[2])
  WW[whichmin]=-1
  diag(WW)=-colSums(WW)+.01 #diagonal add noise
  #choler = choler,eigenvecinv=(aaa$vectors),eigenvalinv = diaginveigen
  #now compute cholesky
  choler=chol(WW)*sqrt(tau)
  bbb=eigen(WW)
  eigval=bbb$values
  eigval[which(eigval<0)]=0
  diaginveigen=tau*diag(bbb$values)
  
  #now we create the new design matrix that includes the fixed design matrix
  KK_mat_new=cbind(design_mat,KK_mat)
  
  filler=matrix(0,design_cols,design_cols)
  supercholer=as.matrix(Matrix::bdiag(filler,choler))
  if(init_calc=="GMRF")
  {
    init_val=as.matrix(rbind(solve(Matrix::t(supercholer)%*%supercholer+Matrix::t(KK_mat_new[,1:(num_r1+design_cols)])%*%KK_mat_new[,1:(num_r1+design_cols)])%*%Matrix::t(KK_mat_new[,1:(num_r1+design_cols)])%*%c(yy),as.matrix(rep(0,(dim(KK_mat_new)[2]-(num_r1+design_cols))))))
  }
  
  #adjust the GMRF to the right size with zeros for the non R1
  cholerbig=as.matrix(Matrix::bdiag(supercholer,diag(rep(0,length(init_val)-dim(supercholer)[1]))))
  
  
  
  eigenvecinvbig=as.matrix(Matrix::bdiag(lava::revdiag(rep(1,design_cols)),as.matrix(Matrix::bdiag(bbb$vectors,lava::revdiag(rep(1,length(init_val)-dim(choler)[1]-design_cols))))))
  diaginveigenbig=as.matrix(Matrix::bdiag(lava::revdiag(rep(0,design_cols)),diaginveigen,lava::revdiag(rep(0,length(init_val)-(dim(choler)[1]+design_cols)))))
  
  
  # init_val=as.matrix(rbind(solve(t(choler)%*%choler+t(KK_mat[,1:num_r1])%*%KK_mat[,1:num_r1])%*%t(KK_mat[,1:num_r1])%*%c(yy),as.matrix(rep(0,(dim(KK_mat)[2]-num_r1)))))
  fitt=multires.fistaTree((matrix(c(yy))),KK_mat_new,init_val,tree,TRUE,numThreads = -1,verbose = F,
                          lambda1 =lambda1, it0 = 10, max_it = 100000,
                          L0 = 0.1, tol = m_tol, intercept = FALSE,
                          pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-l2-tikh',
                          choler = cholerbig,eigenvecinv=eigenvecinvbig,eigenvalinv = diaginveigenbig)
  
  #create sigma^2 initial value
  sigma2=(sum((c(yy)-(KK_mat_new)%*%Matrix::Matrix(fitt[[1]]))^2)+tau*Matrix::t(fitt[[1]])%*%Matrix::t(cholerbig)%*%cholerbig%*%fitt[[1]])/(length(yy)+2)
  
  results=list()
  results[[1]]=list(current_design=KK_mat,fit=fitt,tree=tree,current_knots=current_knots,lambda1=lambda1,tau=tau,sigma2=sigma2)
  #maximum iters
  start=proc.time()
  for(iter in 1:(maxiter-1))
  {
    #update knots
    #first find nonzero non res1 knots with no children yet
    #we can change this to all knots if we want maybe to test
    #nonzero indicies
    nonzero_coef=which(results[[iter]]$fit[[1]]!=0)
    #extract knots
    if(design_cols==0)
    {
      nonzero_nondesign_coef=nonzero_coef
    }
    else{
      #need to adjust it for design matrix
      nonzero_nondesign_coef=nonzero_coef[-(1:design_cols)]-design_cols
    }
    
    
    nonzero_knots=results[[iter]]$current_knots[nonzero_nondesign_coef,]
    #find the ones without children, and get the row numbers
    nonzero_nochild_rows=nonzero_knots[which(apply(nonzero_knots[,c(5:8)],1,sum,na.rm=T)==0),1]
    #now we update the knot matrix
    prevknotct=nrow(current_knots)
    for(i in nonzero_nochild_rows)
    {
      #update children
      current_loc=which(current_knots[,1]==i)
      current_knots[which(current_knots[,5]>current_loc),5]=current_knots[which(current_knots[,5]>current_loc),5]+2^spatial_dimension
      current_knots[which(current_knots[,6]>current_loc),6]=current_knots[which(current_knots[,6]>current_loc),6]+2^spatial_dimension
      current_knots[which(current_knots[,7]>current_loc),7]=current_knots[which(current_knots[,7]>current_loc),7]+2^spatial_dimension
      current_knots[which(current_knots[,8]>current_loc),8]=current_knots[which(current_knots[,8]>current_loc),8]+2^spatial_dimension
      #handles rbind for last row
      if(i==prevknotct){lastrow=T}
      else{lastrow=F}
      
      #EDIT FOR 2D
      current_knots[current_loc,5]=current_loc+1
      current_knots[current_loc,6]=current_loc+2
      if(spatial_dimension==1)
      {
        aa=.create_knots_1d(current_knots[current_loc,3],current_knots[current_loc,2],2,min(dist(knots_r1))*.5^(0:100))
        
        if(lastrow)
        {
          newrow=c(NA,aa[1,2],aa[1,1],NA,NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots,newrow)
          newrow=c(NA,aa[2,2],aa[2,1],NA,NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots,newrow)
        }
        else
        {
          newrow=c(NA,aa[1,2],aa[1,1],NA,NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots[1:current_loc,],newrow,
                              current_knots[((current_loc)+1):nrow(current_knots),])
          newrow=c(NA,aa[2,2],aa[2,1],NA,NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots[1:current_loc,],newrow,
                              current_knots[((current_loc)+1):nrow(current_knots),])
          
        }
      }
      if(spatial_dimension==2)
      {
        current_knots[current_loc,7]=current_loc+3
        current_knots[current_loc,8]=current_loc+4
        aa=.create_knots_2d(matrix(current_knots[current_loc,3:4],1,2),current_knots[current_loc,2],2*spatial_dimension,min(dist(knots_r1))*.5^(0:10))
        if(lastrow)
        {
          newrow=c(NA,aa[1,3],aa[1,1],aa[1,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots,newrow)
          newrow=c(NA,aa[2,3],aa[2,1],aa[2,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots,newrow)
          newrow=c(NA,aa[3,3],aa[3,1],aa[3,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots,newrow)
          newrow=c(NA,aa[4,3],aa[4,1],aa[4,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots,newrow)
        }
        else
        {
          newrow=c(NA,aa[1,3],aa[1,1],aa[1,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots[1:current_loc,],newrow,
                              current_knots[((current_loc)+1):nrow(current_knots),])
          newrow=c(NA,aa[2,3],aa[2,1],aa[2,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots[1:current_loc,],newrow,
                              current_knots[((current_loc)+1):nrow(current_knots),])
          newrow=c(NA,aa[3,3],aa[3,1],aa[3,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots[1:current_loc,],newrow,
                              current_knots[((current_loc)+1):nrow(current_knots),])
          newrow=c(NA,aa[4,3],aa[4,1],aa[4,2],NA,NA,NA,NA,NA)
          current_knots=rbind(current_knots[1:current_loc,],newrow,
                              current_knots[((current_loc)+1):nrow(current_knots),])       
        }
      }
      
      
    }
    #fix row assignments
    current_knots[,1]<-1:nrow(current_knots)
    #create design matrix
    lister=list()
    for(i in 1:dim(current_knots)[1])
    {
      XX_col=current_knots[i,9]
      if(is.na(XX_col))
      {
        if(spatial_dimension==1)
        {
          knott=current_knots[i,3]
        }
        else
        {
          knott=matrix(current_knots[i,3:4],1,2)
        }
        ress=current_knots[i,2]
        lister[[i]]=as(spam::as.dgCMatrix.spam(.create_column(knott,ress,locationss,max_dist,nu,"bezier")),"sparseVector")
      }
      else
      {
        lister[[i]]=as(results[[iter]]$current_design[,XX_col],"sparseVector")
      }
    }
    KK_mat=sv.cbind(lister)
    
    #now we update the group matrix
    
    resolutions_non1=current_knots[which(current_knots[,2]!=1),2]
    n_groups=1+sum(current_knots[,2]!=1)
    n_r1=sum(current_knots[,2]==1)
    eta_g=c(10^-100,shrinkage[resolutions_non1])
    
    groups=matrix(0,n_groups,n_groups)
    tracker=1
    for(i in 1:nrow(current_knots))
    {
      if(current_knots[i,2] > 1)
      {
        tracker=tracker+1
      }
      children=current_knots[i,5:8]
      if(sum(is.na(children))==4)
      {
        next
      }
      children=children[!is.na(children)]
      groups[children-n_r1+1,tracker]=1
    }
    
    groups=Matrix::Matrix(groups)
    groups=groups>0
    
    
    
    
    #remember a minus 1 is needed since C is 0 indexed
    own_variables=as.integer(c(1,current_knots[which(current_knots[,2]>1),1]+design_cols)-1)
    if(spatial_dimension==1)
    {
      N_own_variables=as.integer(c(length(knots_r1)+design_cols,rep(1,n_groups-1)))
    }
    if(spatial_dimension==2)
    {
      N_own_variables=as.integer(c(nrow(knots_r1)+design_cols,rep(1,n_groups-1)))
    }
    tree=list(eta_g=eta_g,
              groups=groups,
              own_variables=own_variables,
              N_own_variables=N_own_variables)
    #get initial value using previous iteration, will need to pad zeros since knots aren't in same positions
    old_fit=results[[iter]]$fit[[1]]
    
    KK_mat_new=cbind(design_mat,KK_mat)
    
    init_val=matrix(0,dim(KK_mat_new)[2],1)
    init_val[1:(num_r1+design_cols)]=old_fit[1:(num_r1+design_cols)]
    for(i in (num_r1+design_cols+1):dim(init_val)[1])
    {
      if(!is.na(current_knots[i-design_cols,9]))
      {
        init_val[i,]=old_fit[current_knots[i-design_cols,9]]
      }
    }
    
    current_knots[,9]=1:dim(current_knots)[1]
    
    
    
    
    
    #adjust the GMRF to the right size with zeros for the non R1
    
    
    
    
    diaginveigenbig=as.matrix(Matrix::bdiag(lava::revdiag(rep(0,design_cols)),diaginveigen,lava::revdiag(rep(0,length(init_val)-(dim(choler)[1]+design_cols)))))
    
    
    
    
    cholerbig=as.matrix(Matrix::bdiag(supercholer,diag(rep(0,length(init_val)-(dim(choler)[1]+design_cols)))))
    eigenvecinvbig=as.matrix(Matrix::bdiag(lava::revdiag(rep(1,design_cols)),bbb$vectors,lava::revdiag(rep(1,length(init_val)-(dim(choler)[1]+design_cols)))))
    diaginveigenbig=as.matrix(Matrix::bdiag(lava::revdiag(rep(0,design_cols)),diaginveigen,lava::revdiag(rep(0,length(init_val)-(dim(choler)[1]+design_cols)))))
    init_val1=init_val
    #correction to speed up initial if no knots have been added
    if(length(init_val)==length(results[[1]]$fit[[1]]))
    {
      init_val1=as.matrix(rbind(solve(Matrix::t(supercholer)%*%supercholer+Matrix::t(KK_mat_new[,1:(num_r1+design_cols)])%*%KK_mat_new[,1:(num_r1+design_cols)])%*%Matrix::t(KK_mat_new[,1:(num_r1+design_cols)])%*%c(yy),as.matrix(rep(0,(dim(KK_mat_new)[2]-num_r1-design_cols)))))
    }
    fitt=multires.fistaTree((matrix(c(yy))),KK_mat_new,init_val1,tree,TRUE,numThreads = -1,verbose = F,
                            lambda1 = lambda1, it0 = 10, max_it = 10000,
                            L0 = 0.1, tol = m_tol, intercept = FALSE,
                            pos = FALSE,compute_gram = TRUE, loss = 'square',regul = 'tree-l2-tikh',
                            choler = cholerbig,eigenvecinv=eigenvecinvbig,eigenvalinv = diaginveigenbig)
    #create sigma^2 initial value
    sigma2=(sum((c(yy)-KK_mat_new%*%fitt[[1]])^2)+tau*Matrix::t(fitt[[1]])%*%Matrix::t(cholerbig)%*%cholerbig%*%fitt[[1]])/(length(yy)+2)
    
    
    results[[iter+1]]=list(current_design=KK_mat,fit=fitt,tree=tree,
                           current_knots=current_knots,lambda1=lambda1,tau=tau,runtime=(proc.time()-start),sigma2=sigma2)
    print(iter)
    print(lambda1)
    print(tau)
    print(sigma2)
    print(proc.time()-start)
    #now we must decide if we need to update tau (when we've maximized beta for a particular tau), or if tau has converged and we
    #need a new lambda
    if(length(results[[iter+1]]$fit[[1]])==length(results[[iter]]$fit[[1]]))
    {
      times=times+1
      if(times==1)
      {
        next
      }
      #if we get to the next else, then beta has converged
      else
      {
        times=0
        #now we decide if tau is still changing
        #we do this by making a new tau and seeing if it is the same as the old tau
        if(em_tol!=-1)
        {
          #update R1 tau
          r1_index=current_knots[which(current_knots[,2]==1),9]
          r1coef=fitt[[1]][r1_index+design_cols]
          tau_newa=tau_a+length(r1coef)/2
          tau_newb=tau_b+t(r1coef)%*%WW%*%(r1coef)/sigma2
          tau=c(tau_newa/tau_newb)
          
          #now update cholesky
          choler=chol(WW)*sqrt(tau)
          bbb=eigen(WW)
          eigval=bbb$values
          eigval[which(eigval<0)]=0
          diaginveigen=tau*diag(bbb$values)
          #if new and old tau are the same, then we need a new lambda or to break
          
          if(abs(tau-results[[iter+1]]$tau)<em_tol)
          {
            if(tree_it==length(lambdatree_seq))
            {
              break
            }
            else
            {
              tree_it=tree_it+1
              lambda1=lambdatree_seq[tree_it]
            }
            
          }
        }
        else
        {
          if(tree_it==length(lambdatree_seq))
          {
            break
          }
          else
          {
            tree_it=tree_it+1
            lambda1=lambdatree_seq[tree_it]
          }
        }
        
        
      }
    }
  }
  return(list(results=results,params=params,total_iter=(iter+1),runtime=proc.time()-start))
}

.create_knots_1d<-function (knots, resolutions, p, resolution_dist) 
{
  repknot = knots[rep(1:length(knots), each = p)]
  temp = 0.5 * resolution_dist[c(resolutions + 1)]
  tempx = rep(temp, each = p) * rep(c(-1, 1), length(knots))
  return(cbind(repknot + tempx, rep(resolutions + 1, each = p)))
}

.create_knots_2d<-function (knots, resolutions, p, resolution_dist) 
{
  repknot = knots[rep(1:nrow(knots), each = p), ]
  temp = 0.5 * resolution_dist[c(resolutions + 1)]
  tempx = rep(temp, each = p) * rep(c(-1, -1, 1, 1), nrow(knots))
  tempy = rep(temp, each = p) * rep(c(-1, 1, -1, 1), nrow(knots))
  return(cbind(repknot + cbind(tempx, tempy), rep(resolutions + 
                                                    1, each = p)))
}