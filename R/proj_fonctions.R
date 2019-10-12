####-------------- If not installed, we need the here package
#install.packages("here")

####------------- Packages
library(glmnet)
library(ncvreg)
library(spatstat)
library(readxl)

########################################################################################################################################
####------------ Function for estimating the parameters of the models (Strauss, Geyer and Area) using regularization techniques ---#####
###------------------ (Lasso, Ridge, Elastic Net, Adaptive Lasso, Adaptive Elastic Net, SCAD and MCP) or not ----------------------#####
########################################################################################################################################
gpp_reg <- function(pp, Qimage, int, method='nopen', tuning='BIC', f.dummy=2){
  
  # Range and area of the observation window
  int.x <- pp$window$xrange
  int.y <- pp$window$yrange
  D <- (int.x[2]-int.x[1])*(int.y[2]-int.y[1])
  
  # Covariates list and names
  cov.list <- NULL
  lstnames <- 'x1'  
  
  # Standardize the covariates
  Qmeans <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, mean)
  Qsd <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, sd)
  cov.list[[1]] <- as.im((Qimage[,,1]-Qmeans[1])/Qsd[1], W=as.owin(c(int.x,int.y)))
  for(k in 2:length(Qimage[1,1,])){
    lstnames <- c(lstnames, paste('x', k, sep = ''))
    cov.list[[k]] <- as.im((Qimage[,,k]-Qmeans[k])/Qsd[k], W=as.owin(c(int.x,int.y)))
  }
  names(cov.list) <- lstnames
  
  # Fitting Gibbs 
  n <- npoints(pp)
  n.dummy <- round(f.dummy*sqrt(n))
  temp <- ppm(pp ~ ., interaction=int, data=cov.list, nd=n.dummy)
  temp_glm <- temp$internal$glmfit$data[,1:(length(temp$coef)+1)]
  yy <- temp_glm[,2]
  N <- length(yy)
  wts <- temp_glm[,1]
  Q <- as.matrix(temp_glm[,c(-1,-2)])
  par.ini <- coef(temp)
  ini.theta <- par.ini[-1]
  n1 <- length(ini.theta)
  tot.cov <- dim(Qimage)[3]
  
  # Variable selection procedures
  if(method == "nopen"){
    return(par.ini)}
  
  if(method == "lasso"){
    fit = glmnet(Q, yy, family="poisson", alpha=1, weights=wts, penalty.factor=c(rep(D,tot.cov),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "ridge"){
    fit = glmnet(Q, yy, family="poisson", alpha=0, weights=wts, penalty.factor=c(rep(D,tot.cov),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "enet"){
    fit = glmnet(Q, yy, family="poisson", alpha=0.5, weights=wts, penalty.factor=c(rep(D,tot.cov),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "al"){
    fit = glmnet(Q, yy, family="poisson", alpha=1, weights=wts, penalty.factor=c(D/abs(ini.theta[-n1]),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "aenet"){
    fit = glmnet(Q, yy, family="poisson", alpha=.5, weights=wts, penalty.factor=c(D/abs(ini.theta[-n1]),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "scad"){
    fit = ncvreg(Q, yy, family="poisson", penalty="SCAD", wts=wts, penalty.factor=c(rep(D,tot.cov),0))
    coef.theta = as.matrix(fit$beta)
    lambda = as.matrix(fit$lambda)
  }
  if(method == "mcp"){
    fit = ncvreg(Q, yy, family="poisson", penalty="MCP", wts=wts, penalty.factor=c(rep(D,tot.cov),0))
    coef.theta = as.matrix(fit$beta)
    lambda = as.matrix(fit$lambda)
  }
  
  # Tuning parameter selection
  optim.inf.criteria = 1e30
  index.optim.theta = NULL
  Q.new <- cbind(1,Q)
  H.theta <- vcov(temp, what="fisher", hessian=TRUE, fine=TRUE)
  vcov.essai <- try(V.theta <- vcov(temp, fine=TRUE), silent=TRUE)
  if (class(vcov.essai)[1] == "try-error") V.theta <- Emp_vcov(Temp=temp,cov.nb=tot.cov,int.m=int,cov.list=cov.list,nb.it=100,f.dum=f.dummy)
  for (i in 1:length(lambda))
  {
    theta.cur = coef.theta[,i]
    ind.zero = which(theta.cur==0)
    logpseudolike = sum(wts*(yy*(Q.new%*%theta.cur)-exp(Q.new%*%theta.cur)))
    V <- V.theta
    H <- H.theta
    if(length(ind.zero)!=0){
      V = V[-ind.zero,-ind.zero]
      H = H[-ind.zero,-ind.zero]
    }
    if (tuning == "BIC")     {inf.criteria = -2*logpseudolike + (sum(diag(H%*%V)))*log(n)}
    if (tuning == "ERIC")    {inf.criteria = -2*logpseudolike + (sum(diag(H%*%V)))*log((n/(N*lambda[i])))}
    
    if(inf.criteria < optim.inf.criteria){
      index.optim.theta = i
      optim.inf.criteria = inf.criteria
    }
    
  }
  theta.cur = coef.theta[,index.optim.theta]
  lambda.opt = lambda[index.optim.theta]
  
  return(list(theta.reg=coef.theta, theta=theta.cur, lambda=lambda, lambda.opt=lambda.opt))
}



####------------ Computation of vcov when the number of quadrature points is very large i.e. vcov of spatstat fails.
Emp_vcov <- function(Temp,cov.nb,int.m,cov.list,nb.it,f.dum){
  var.cov <- matrix(NA, nrow = nb.it, ncol = cov.nb+2)
  for(i in 1:nb.it){
    X.sim <- rmh(Temp, verbose=FALSE)
    n.dum.sim <- round(f.dum*sqrt(npoints(X.sim)))
    var.cov[i,] <- ppm(X.sim ~ ., interaction=int.m, data=cov.list, nd=n.dum.sim)$coef
  }
  return(var(var.cov))
}


####------------------ Function for simulation in Scenario I and Scenario II
Sim.GPP <- function(Model="strauss",number.iterations=10,beta0,par.interact=.2,R=9.25,wdow=square(1),tr,ns=1000,
                           f.dum=2,satur=NULL,Qim,int.mod){
  if(Model=="strauss") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R), w=wdow, trend=tr)
  if(Model=="geyer") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R, sat=satur), w=wdow, trend=tr)
  if(Model=="areaint") mod <- list(cif=Model, par=list(beta=exp(beta0), eta=par.interact, r=R), w=wdow, trend=tr)
  Theta_Lasso_BIC <-  Theta_Lasso_ERIC <-  Theta_Ridge_BIC <- Theta_Ridge_ERIC <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)
  Theta_Enet_BIC <-  Theta_Enet_ERIC <-  Theta_ALasso_BIC <- Theta_ALasso_ERIC <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)  
  Theta_AEnet_BIC <-  Theta_AEnet_ERIC <- Theta_Scad_BIC <- Theta_Scad_ERIC  <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)
  Theta_Mcp_BIC <-  Theta_Mcp_ERIC  <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)
  nbre.points <- rep(NA,number.iterations)
  for(i in 1:number.iterations){
    X <- rmh(mod, start=list(n.start=ns), control=list(nrep=1e6), verbose=FALSE)
    nbre.points[i] <- npoints(X)
    Theta_Lasso_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="lasso", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_Lasso_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="lasso", tuning="ERIC", f.dummy=f.dum)[[2]]
    Theta_Ridge_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="ridge", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_Ridge_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="ridge", tuning="ERIC", f.dummy=f.dum)[[2]]
    Theta_Enet_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="enet", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_Enet_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="enet", tuning="ERIC", f.dummy=f.dum)[[2]]
    Theta_ALasso_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="al", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_ALasso_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="al", tuning="ERIC", f.dummy=f.dum)[[2]]
    Theta_AEnet_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="aenet", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_AEnet_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="aenet", tuning="ERIC", f.dummy=f.dum)[[2]]
    Theta_Scad_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="scad", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_Scad_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="scad", tuning="ERIC", f.dummy=f.dum)[[2]]
    Theta_Mcp_BIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="mcp", tuning="BIC", f.dummy=f.dum)[[2]]
    Theta_Mcp_ERIC[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="mcp", tuning="ERIC", f.dummy=f.dum)[[2]]
  }
  m1 <- mean(nbre.points)
    return(list(l1=Theta_Lasso_BIC,l2=Theta_Lasso_ERIC,l3=Theta_Ridge_BIC,l4=Theta_Ridge_ERIC,l5=Theta_Enet_BIC,
                l6=Theta_Enet_ERIC,l7=Theta_ALasso_BIC,l8=Theta_ALasso_ERIC,l9=Theta_AEnet_BIC,l10=Theta_AEnet_ERIC,
                l11=Theta_Scad_BIC,l12=Theta_Scad_ERIC,l13=Theta_Mcp_BIC,l14=Theta_Mcp_ERIC,l15=m1))
}

####------------------ Function for simulation in Scenario 0
Sim.GPP.Nopen <- function(Model="strauss",number.iterations=10,beta0,par.interact=.2,R=9.25,wdow=square(1),tr,ns=1000,
                    f.dum=2,satur=NULL,Qim,int.mod){
  if(Model=="strauss") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R), w=wdow, trend=tr)
  if(Model=="geyer") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R, sat=satur), w=wdow, trend=tr)
  if(Model=="areaint") mod <- list(cif=Model, par=list(beta=exp(beta0), eta=par.interact, r=R), w=wdow, trend=tr)
  Theta_Nopen  <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)
  nbre.points <- rep(NA,number.iterations)
  for(i in 1:number.iterations){
    X <- rmh(mod, start=list(n.start=ns), control=list(nrep=1e6), verbose=FALSE)
    nbre.points[i] <- npoints(X)
    Theta_Nopen[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method="nopen", f.dummy=f.dum)
  }
  return(list(l1=Theta_Nopen,l2=mean(nbre.points)))
}



#############################################################################################

## ---- Standardize the covariates ----
Standardize.cov <- function(Qimage, Wdw){
  cov.list <- NULL
  Qmeans <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, mean)
  Qsd <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, sd)
  for(k in 1:length(Qimage[1,1,])){
    cov.list[[k]] <- as.im((Qimage[,,k]-Qmeans[k])/Qsd[k], W=Wdw)
  }
  return(cov.list)
} 


## ---- Variables simulées ----
sum_vec_array <- function(Vect,Arr){
  s <- 0
  for(i in 1:length(Vect)) s <- s + Vect[i]*Arr[,,i]
  return(s)
}
built_scenario1 <- function(W, Nbre.points, fact.col=0.7){
  p <- floor(3*area(W)^.25)
  Qimage <- array(0, dim=c(101,201,p))
  Qimage[,,1] <- bei.extra$elev$v
  Qimage[,,2] <- bei.extra$grad$v
  for(i in 3:p) { Qimage[,,i] <- rnoise(rnorm,dimyx=c(101,201),w=W,mean=0,sd=1)$v }
  Sigma <- matrix(NA, nrow=p, ncol=p)
  for(i in 1:p) Sigma[i,] <- fact.col^(abs(i-(1:p)))
  Sigma[1,2] <- Sigma[2,1] <- 0
  #Sigma.svd <- svd(Sigma)
  #V <- t(Sigma.svd$u%*%sqrt(diag(Sigma.svd$d)))
  V <- chol(Sigma)
  if(all.equal(Sigma,t(V)%*%V)){
    Qim <- array(0, dim=c(101,201,p))
    for(i in 1:p) { Qim[,,i] <- sum_vec_array(V[,i],Qimage) }
    Qim.cr <- Standardize.cov(Qim,W)
    b <- round(log((Nbre.points/integral(exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]]),W))),4)
    trend.function <- exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]])
  }
  else return("erreur de decomposition")
  
  return(list(l0=p,l1=b,l2=Qim,l3=trend.function))
}

## ---- Variables réelles ----
soil_nut <- read_xls(here::here("data","bci.block20.data.xls"),sheet=2)
soil_nut_data <- as.data.frame(soil_nut)
soil_nut_data_im <- soil_nut_data[,-c(1,2)]
## conversion of the 13 pixel images of soil nutrients
soil_nut_im <- list()
for (i in 1:length(soil_nut_data_im)) {
  M <- matrix(soil_nut_data_im[,i],nrow=25,ncol=50)
  z <- im(M,xrange=c(-2.5,1002.5),yrange=c(-2.5,502.5),unitname="metres")
  soil_nut_im[[i]] <- as.im(z,W=commonGrid(bei.extra$elev,z))
}
built_scenario2 <- function(W, Nbre.points){
  p <- floor(3*area(W)^.25)
  Qimage <- array(0, dim=c(101,201,p))
  Qimage[,,1] <- bei.extra$elev$v
  Qimage[,,2] <- bei.extra$grad$v
  for(i in 3:15) Qimage[,,i] <- soil_nut_im[[(i-2)]]$v
  k <- 16; j <- 1; q <- p - 15
  while(j <= q){
    l <- 1
    while(l <= length(soil_nut_im)){
      i <- l + 1
      while( (i <= length(soil_nut_im)) && (k <= p) ){
        intermed <- soil_nut_im[[l]]*soil_nut_im[[i]]
        Qimage[,,k] <- intermed$v
        i <- i + 1
        k <- k + 1
      }
      if( (i > length(soil_nut_im)) ) l <- l + 1
      else l <- length(soil_nut_im) + 1
    }
    if( k > p) j <- q + 1
    else j <- j + 1
  }
  Qimage.cr <- Standardize.cov(Qimage,W)
  b <- round(log((Nbre.points/integral(exp(2*Qimage.cr[[1]]+0.75*Qimage.cr[[2]]),W))),4)
  trend.function <- exp(2*Qimage.cr[[1]]+0.75*Qimage.cr[[2]])
  return(list(l0=p,l1=b,l2=Qimage,l3=trend.function))
}
  
####################################################################################################################

## ---- Selection Performance ----
# Computation of TPR, True positive rate
TPR <- function(Vect){
  tpr <- sum(Vect[c(2,3)]!=0)/2
  return(tpr)
}

# Computation of FPR, False positive rate
FPR <- function(Vect){
  p <- length(Vect)
  p1 <- p - 4
  fpr <- sum(Vect[4:(p-1)]!=0)/p1
  return(fpr)
}

# Computation of PPV, Positive predictive value
PPV <- function(Vect){
  p <- length(Vect)
  nb.select.true.cov <- sum(Vect[c(2,3)]!=0)
  tot.nb.select.cov <- sum(Vect[-c(1,p)]!=0)
  ppv <-  nb.select.true.cov/tot.nb.select.cov
  return(ppv)
}

########################################################################################################################


## ---- Properties of the estimate ----
## Properties of the estimate: Bias, Variance, MSE, TPR and FPR
Estimate.Properties <- function(Est.theta, init.theta, method="nopen") {
  Bias.est <- round(sqrt(sum((apply(Est.theta,2,mean)[-1] - init.theta[-1])^2)),2)
  SD.est <- round(sqrt(sum(apply(Est.theta,2,var)[-1])),2)
  RMSE.est <- round(sqrt(Bias.est^2 + SD.est^2),2)
  if(method=="nopen") {return(list(Bias=Bias.est, SD=SD.est, RMSE=RMSE.est))}
  else{
    FPR.est <- round(100*mean(apply(Est.theta, 1, FPR)))
    TPR.est <- round(100*mean(apply(Est.theta, 1, TPR)))
    #PPV.est <- round(100*mean(apply(Est.theta, 1, PPV)))
    return(list(Bias=Bias.est, SD=SD.est, RMSE=RMSE.est, FPR=FPR.est, TPR=TPR.est))
  }
}






