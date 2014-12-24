require(snowfall)
require(stochvol)
#-------------a-Gibbs-related prelims-------------------------
p <- 2#number of lags of the dependent variable
nsave <- 5000 # final number of draws to save
nburn <- 5000 #Burn-INs
ntot <- nsave+nburn # total draws
#------------------------------------------------------------
#Standardize data
#factors_dat <- lapply(fac_list,standardize)
#now retrieve group specific factors and create a large X matrix

#Xraw <-  read.table("~/Dropbox/UNI/PhD/Fiscal Multiplier in Financial Crises/R_project/dataset.csv", sep=";", quote="\"")#read.table("C:/Users/Flo/Google Drive/UNI/PhD/Fiscal Multiplier in Financial Crises/R_project/dataset.csv", sep=";", quote="\"")
Xraw <- read.csv("~/Dropbox/efficient_tvp_sv/mudata.csv", header=FALSE)
Xraw <- as.matrix(Xraw[,1:4])

class(Xraw) <- "numeric"
Y <- Xraw
#Create lagged Y matrix
Xlag <- mlag(Y,p)

Y <- Y[(p+1):nrow(Y),]
X <- cbind(1,Xlag[(p+1):nrow(Xlag),])

M <- ncol(Y);K <- ncol(X);T <- nrow(X) #get dimensions of X and Y
A0 <- diag(M)
A1 <- matrix(0,K,M)
#Prior setup
b0 <- matrix(0,K,M)
b0[1:M,1:M] <- diag(M)
b0 <- rbind(matrix(0,M-1,M),b0)

B0 <- diag(K+(M-1))*1/100
B0inv <- solve(B0)


#Storage matrices for the full system
H_store <- array(0,c(T,M,nsave))
ALPHA_store <- array(0,c(T,K,M,nsave))
A0_store <- array(0,c(T,M,M,nsave))
Ht <- kronecker(matrix(1,T,1),diag(M)*0.1)

post_draws <- list()
sfInit(parallel=TRUE,cpus=4)
sfExport(list=list("Y","X","nburn","nsave","B0inv","b0","bvartvpm","KF"))
post_draws <- sfLapply(1:ncol(Y),function(i) bvartvpm(i,Y,X,nburn,nsave,B0inv,b0))
sfStop()

for (ii in 1:nsave){
  for (jj in 1:M){
    slct <- post_draws[[jj]]$slct
    #split and create structural matrix
    A0_store[,jj,slct,ii] <- (-1*post_draws[[jj]]$A[,slct,ii])
    if (jj==1){
      ALPHA_store[,,jj,ii] <- (post_draws[[jj]]$A[,,ii])
    }else{
      ALPHA_store[,,jj,ii] <- (post_draws[[jj]]$A[,-slct,ii])
    }
  }
}
A_mean <- apply(ALPHA_store,c(1,2,3),mean)
A0_mean <- apply(A0_store,c(1,2,3),mean)
H_mean <- apply(H_store,c(1,2),mean)
fit <- NULL
S_post <- array(0,c(T,M,M))


for (jj in 1:T){
  A0 <- A0_mean[jj,,]
  diag(A0) <- 1
  Atilda <- t(solve(A0)%*%t(A_mean[jj,,]))
  fit <- rbind(fit,X[jj,]%*%Atilda)
  S_post[jj,,] <- solve(A0)%*%diag(exp(H_mean[jj,]))%*%t(solve(A0))
}
plot(Y[,3])
lines(fit[,3])




