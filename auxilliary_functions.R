#Inputs needed
#Y , X, nr, ntot, nsave, B0inv, b0
bvartvpm <- function(nr,Y,X,nburn,nsave,B0inv,b0){
  require(bayesm)
  require(stochvol)
  #Get data/model dimensions (full TVP-VAR)
  T <- nrow(Y)
  K <- ncol(X)
  M <- ncol(Y)
  Rcpp::sourceCpp('carter_kohn_cpp.cpp')
  #Prior stuff for SV
  svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,T))
  hv <- svdraw$latent
  para <- list(mu=-10,phi=.9,sigma=.2)
  H <- matrix(-10,T,M)
  #prs on state eq for alpha
  Qprmean <- 0.001
  Q_prvar <- 40
  #Create full-data matrices
  if (nr==1) slct <- NULL else slct <- 1:(nr-1)
  Y__ <- Y[,nr]
  X__ <- cbind(Y[,slct],X)
  K_ <- ncol(X__)
  M_ <- M-length(slct)
  
  #storage matrices
  H_store <- matrix(0,T,nsave)
  ALPHA_store <- array(0,c(T,K_,nsave))
  u_ <- matrix(0,T,1)
  ntot <- nburn+nsave
  for (irep in 1:ntot){
    #------------------ normalize data-------------------#
    #rescaled data, only needed in the constant VAR case
    X_ <- X__*as.numeric(exp(-hv/2))
    Y_ <- Y__*as.numeric(exp(-hv/2))
    
    if (irep==1) Qdraw <- B0inv[(M_):(K+M-1),(M_):(K+M-1)]/1000
    ALPHA <- KF(t(as.matrix(Y__)),X__,as.matrix(exp(hv)),Qdraw,K_,1,T,b0[(M_):(K+M-1),nr,drop=FALSE],B0inv[(M_):(K+M-1),(M_):(K+M-1)])
    
    # Take the SSE in the state equation of ALPHA_t
    Atemp  <-  t(ALPHA[,2:T]) - t(ALPHA[,1:(T-1)])
    SSE  <-  matrix(0,K_,K_)
    for (i in 1:(T-1)){
      SSE  <-  SSE + t(Atemp[i,,drop=FALSE])%*%Atemp[i,,drop=FALSE]
    }
    # ...and subsequently draw Q, the covariance matrix of B(t)
    Qinv  <-  solve(SSE + diag(K_)*Qprmean)
    Qinvdraw  <-  rwishart(T+Q_prvar,Qinv)$W
    Qdraw  <-  solve(Qinvdraw) # this is a draw from Q
    
    for (i in 1:T){
      u_[i,1] <- Y__[i]-X__[i,]%*%ALPHA[,i]
    }
    
    #-----------------sample log volas-------------------#
    svdraw <- .svsample(u_,startpara=para(svdraw),startlatent=latent(svdraw),priorphi=c(25,1.5))
    hv <- latent(svdraw)
    if (irep>nburn){
      H_store[,irep-nburn] <- hv
      ALPHA_store[,,irep-nburn] <- t(ALPHA)
    }
    print(irep)
    
  }
  return(list(A=ALPHA_store,H=H_store,slct=slct))
}
mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag) 
}

dmean <- function(x){
  
  xnew <- (x-mean(x))/sd(x)
  return(xnew)
}


