library(glmnet)
library(spatstat)

gpp_reg <- function(pp, Qimage, int, method='nopen', tuning='BIC', f.dummy=2){
  
  # Range of the observation window
  int.x <- pp$window$xrange
  int.y <- pp$window$yrange
  
  # Covariates list and names
  cov.list <- NULL
  rhs <- 'x1'
  lstnames <- 'x1'  
  
  # Standardize the covariates
  Qmeans <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, mean)
  Qsd <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, sd)
  cov.list[[1]] <- as.im((Qimage[,,1]-Qmeans[1])/Qsd[1], W=as.owin(c(int.x,int.y)))
  for(k in 2:length(Qimage[1,1,])){
    rhs <- paste(rhs, '+', 'x', k, sep = '')
    lstnames <- c(lstnames, paste('x', k, sep = ''))
    cov.list[[k]] <- as.im((Qimage[,,k]-Qmeans[k])/Qsd[k], W=as.owin(c(int.x,int.y)))
  }
  names(cov.list) <- lstnames
  
  # Fitting Gibbs 
  n <- npoints(pp)
  n.dummy <- round(f.dummy*sqrt(n))
  temp <- ppm(pp, trend=as.formula(paste('~',rhs)), interaction=int, covariates=cov.list, nd=n.dummy) 
  temp_glm <- temp$internal$glmfit$data[,1:(length(temp$coef)+1)]
  yy <- temp_glm[,2]
  N <- length(yy)
  wts <- temp_glm[,1]
  Q <- as.matrix(temp_glm[,c(-1,-2)])
  par.ini <- coef(temp)
  theta0.hat <- par.ini[1]
  ini.theta <- par.ini[-1]
  n1 <- length(ini.theta)
  tot.cov <- dim(Qimage)[3]
  
  # Variable selection procedures
  if(method == "nopen"){
    return(list(theta=par.ini))}
  
  if(method == "lasso"){
    fit = glmnet(Q, yy, family="poisson", alpha=1, weights=wts, penalty.factor=c(rep(1,tot.cov),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "ridge"){
    fit = glmnet(Q, yy, family="poisson", alpha=0, weights=wts, penalty.factor=c(rep(1,tot.cov),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "enet"){
    fit = glmnet(Q, yy, family="poisson", alpha=0.5, weights=wts, penalty.factor=c(rep(1,tot.cov),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "al"){
    fit = glmnet(Q, yy, family="poisson", alpha=1, weights=wts, penalty.factor=c(1/abs(ini.theta[-n1]),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  if(method == "aenet"){
    fit = glmnet(Q, yy, family="poisson", alpha=.5, weights=wts, penalty.factor=c(1/abs(ini.theta[-n1]),0))
    coef.theta = as.matrix(rbind(fit$a0, fit$beta))
    lambda = as.matrix(fit$lambda)
  }
  
  # Tuning parameter selection
  lambda.min <- min(lambda)
  lambda.max <- max(lambda)
  optim.inf.criteria = 1e30
  index.optim.theta = NULL
  Q.new <- cbind(1,Q)
  H.theta <- vcov(temp, hessian=TRUE)
  if(n > 2500) V.theta <- Emp_vcov(Temp=temp,cov.nb=tot.cov,trend.f=as.formula(paste('~',rhs)),int.m=int,    
                                   cov.list=cov.list,nb.it=50,f.dum=f.dummy)
  else{V.theta <- vcov(temp)}
  V.theta <- solve(V.theta)
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
    if (tuning == "BIC")     {inf.criteria = -2*logpseudolike + (sum(diag(V%*%H)))*log(n)}
    if (tuning == "ERIC")    {inf.criteria = -2*logpseudolike + (sum(diag(V%*%H)))*log((n/N)*lambda[i])}
    
    if(inf.criteria < optim.inf.criteria){
      index.optim.theta = i
      optim.inf.criteria = inf.criteria
    }
    
  }
  theta.cur = coef.theta[,index.optim.theta]
  lambda.opt = lambda[index.optim.theta]
  
  return(list(theta.reg=coef.theta, theta=theta.cur, lambda=lambda, lambda.opt=lambda.opt))
}

## Main function 
Simulation.GPP <- function(Model="strauss",number.iterations=50,beta0,par.interact=.2,R=12,wdow=square(1),tr,ns=1000,
                           f.dum=2,satur=NULL,Qim,int.mod,met="lasso",inf.criteria="BIC"){
  if(Model=="strauss") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R), w=wdow, trend=tr)
  if(Model=="geyer") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R, sat=satur), w=wdow, trend=tr)
  if(Model=="areaint") mod <- list(cif=Model, par=list(beta=exp(beta0), eta=par.interact, r=R), w=wdow, trend=tr)
  Theta <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)
  for(i in 1:number.iterations){
    X <- rmh(mod, start=list(n.start=ns), control=list(nrep=1e6), verbose=FALSE)
    Theta[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method=met, tuning=inf.criteria, f.dummy=f.dum)[[2]]
  }
  return(Theta)
}

Emp_vcov <- function(Temp,cov.nb,trend.f,int.m,cov.list,nb.it,f.dum){
  var.cov <- matrix(NA, nrow = nb.it, ncol = cov.nb+2)
  for(i in 1:nb.it){
    X.sim <- simulate(Temp, drop=TRUE)
    n.dum.sim <- round(f.dum*sqrt(npoints(X.sim)))
    var.cov[i,] <- ppm(X.sim, trend=trend.f, interaction=int.m, covariates=cov.list, nd=n.dum.sim)$coef
    }
  return(var(var.cov))
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


## ---- Scenario I: Variables simulées ----
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
  Sigma.svd <- svd(Sigma)
  V <- t(Sigma.svd$u%*%sqrt(diag(Sigma.svd$d)))
  if(all.equal(Sigma,t(V)%*%V)){
    Qim <- array(0, dim=c(101,201,p))
    for(i in 1:p) { Qim[,,i] <- sum_vec_array(V[,i],Qimage) }
    Qim.cr <- Standardize.cov(Qim,W)
    b <- round(log((Nbre.points/integral(exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]]),W))),4)
    trend.function <- exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]])
  }
  else return("erreur decomposition")
  
  return(list(l0=p,l1=b,l2=Qim,l3=trend.function))
}

####################################################################################################################

## ---- Selection Performance ----
# Computation of TPR, True positive rate
TPR.Scenario12 <- function(Vect){
  return(sum(Vect[c(2,3)]!=0)/2)
}

# Computation of FPR, False positive rate
FPR.Scenario12 <- function(Vect){
  return(sum(Vect[4:51]!=0)/48)
}

# Computation of FPR, False positive rate for scenario 3
FPR.Scenario3 <- function(Vect){
  return(sum(Vect[4:16]!=0)/13)
}

# Computation of PPV, Positive predictive value
PPV.Scenario12 <- function(Vect){
  n <- length(Vect)
  number.selected.true.covariates <- sum(Vect[c(2,3)]!=0)
  total.number.selected.cov.mod <- sum(Vect[-c(1,n)]!=0)
  return(number.selected.true.covariates/total.number.selected.cov.mod)
}

########################################################################################################################


## ---- Properties of the estimate ----
## Properties of the estimate: Bias, Variance, MSE, TPR and FPR
Estimate.Properties <- function(Est.theta, init.theta, Scenario) {
  Bias.est <- round(sqrt(sum((apply(Est.theta,2,mean)[-1] - init.theta[-1])^2)),2)
  Var.est <- round(sum(apply(Est.theta,2,var)[-1]),2)
  MSE.est <- round(Bias.est^2 + Var.est,2)
  if(Scenario=="1" | Scenario=="2"){ FPR.est <- round(100*mean(apply(Est.theta, 1, FPR.Scenario12)),0)}
  else FPR.est <- round(100*mean(apply(Est.theta, 1, FPR.Scenario3)),0)
  TPR.est <- round(100*mean(apply(Est.theta, 1, TPR.Scenario12)),0)
  return(list(Bias=Bias.est, Var=Var.est, MSE=MSE.est, FPR=FPR.est, TPR=TPR.est))
}






