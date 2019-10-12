library(spatstat)
W <- owin(c(0,1000),c(0,500))
Qi <- array(0, dim=c(101,201,2))
Qi[,,1] <- bei.extra$elev$v
Qi[,,2] <- bei.extra$grad$v
Qi.cr <- Standardize.cov(Qi,W)

cov.list <- NULL
rhs <- 'x1'
lstnames <- 'x1'
Qmeans <- apply(matrix(Qi, prod(dim(Qi[,,1])), length(Qi[1,1,])), 2, mean)
Qsd <- apply(matrix(Qi, prod(dim(Qi[,,1])), length(Qi[1,1,])), 2, sd)
cov.list[[1]] <- as.im((Qi[,,1]-Qmeans[1])/Qsd[1], W=W)
for(k in 2:length(Qi[1,1,])){
  rhs <- paste(rhs, '+', 'x', k, sep = '')
  lstnames <- c(lstnames, paste('x', k, sep = ''))
  cov.list[[k]] <- as.im((Qi[,,k]-Qmeans[k])/Qsd[k], W=W)
}
names(cov.list) <- lstnames

# Fitting Gibbs 
n <- npoints(bei)
n.dummy <- round(8*sqrt(n))
temp <- ppm(bei, trend=as.formula(paste('~',rhs)), covariates=cov.list, nd=n.dummy) 
H <- vcov(temp, what="fisher",hessian=TRUE, fine=TRUE)
VC <- vcov(temp, fine=TRUE)
sum(diag(H%*%VC))

## Lennard Jones pair potential function
len.Jones <- function(r,a,b){
  return(a*r^(-12)+b*r^(-6))
}



Simulation.GPP <- function(Model="strauss",number.iterations=10,beta0,par.interact=.2,R=12,wdow=square(1),tr,ns=1000,
                           f.dum=2,satur=NULL,Qim,int.mod,met="nopen",inf.criteria="BIC"){
  if(Model=="strauss") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R), w=wdow, trend=tr)
  if(Model=="geyer") mod <- list(cif=Model, par=list(beta=exp(beta0), gamma=par.interact, r=R, sat=satur), w=wdow, trend=tr)
  if(Model=="areaint") mod <- list(cif=Model, par=list(beta=exp(beta0), eta=par.interact, r=R), w=wdow, trend=tr)
  Theta <- matrix(NA, nrow=number.iterations,ncol=dim(Qim)[3]+2)
  for(i in 1:number.iterations){
    X <- rmh(mod, start=list(n.start=ns), control=list(nrep=1e6), verbose=FALSE)
    if(met=="nopen") {Theta[i,] <- unlist(gpp_reg(pp=X, Qimage=Qim, int=int.mod, method=met, tuning=inf.criteria, f.dummy=f.dum))}
    else {Theta[i,] <- gpp_reg(pp=X, Qimage=Qim, int=int.mod, method=met, tuning=inf.criteria, f.dummy=f.dum)[[2]]}
  }
  return(Theta)
}


######---------- gamma=0.5 pour strauss model

## Lasso method
Str.BIC.Lasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                     Qim=ScenI1$l2,int.mod=Strauss(12),ns=500,f.dum=2,met="lasso")
Str.ERIC.Lasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                      Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=500,met="lasso")

## Ridge method
Str.BIC.Ridge.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                     Qim=ScenI1$l2,int.mod=Strauss(12),ns=500,f.dum=2,met="ridge")
Str.ERIC.Ridge.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                      Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=500,met="ridge")

## Elastic Net method
Str.BIC.Enet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                    Qim=ScenI1$l2,int.mod=Strauss(12),ns=500,f.dum=2,met="enet")
Str.ERIC.Enet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                     Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=500,met="enet")

## Adaptive Lasso method
Str.BIC.ALasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                      Qim=ScenI1$l2,int.mod=Strauss(12),ns=500,f.dum=2,met="al")
Str.ERIC.ALasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                       Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=500,met="al")

## Adaptive Elastic Net method 
Str.BIC.AEnet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                     Qim=ScenI1$l2,int.mod=Strauss(12),ns=500,f.dum=2,met="aenet")
Str.ERIC.AEnet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                      Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=500,met="aenet")

save(Str.BIC.Lasso.ScI1,Str.ERIC.Lasso.ScI1,Str.BIC.Ridge.ScI1,Str.ERIC.Ridge.ScI1,Str.BIC.Enet.ScI1,
     Str.ERIC.Enet.ScI1,Str.BIC.ALasso.ScI1,Str.ERIC.ALasso.ScI1,Str.BIC.AEnet.ScI1,
     Str.ERIC.AEnet.ScI1, file="data/Strauss.ScenarioI1.RData")

######------------ gamma=1.5 for geyer model

## Lasso method
Gey.BIC.Lasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                     tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=500,f.dum=2,met="lasso")
Gey.ERIC.Lasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=500,met="lasso")

## Ridge method
Gey.BIC.Ridge.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                     tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=500,f.dum=2,met="ridge")

Gey.ERIC.Ridge.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=500,met="ridge")

## Elastic Net method
Gey.BIC.Enet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                    tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=500,f.dum=2,met="enet")
Gey.ERIC.Enet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                     tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                     f.dum=2,ns=500,met="enet")

## Adaptive Lasso method
Gey.BIC.ALasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=500,f.dum=2,met="al")
Gey.ERIC.ALasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                       tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=2,ns=500,met="al")

## Adaptive Elastic Net 
Gey.BIC.AEnet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                     tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=500,f.dum=2,met="aenet")

Gey.ERIC.AEnet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=500,met="aenet")

save(Gey.BIC.Lasso.ScI1, Gey.ERIC.Lasso.ScI1, Gey.BIC.Ridge.ScI1, Gey.ERIC.Ridge.ScI1, 
     Gey.BIC.Enet.ScI1, Gey.ERIC.Enet.ScI1, Gey.BIC.ALasso.ScI1, Gey.ERIC.ALasso.ScI1, 
     Gey.BIC.AEnet.ScI1, Gey.ERIC.AEnet.ScI1, file="data/Geyer.ScenarioI1.RData")


######---------- gamma=0.5 pour strauss model

## Lasso method
Str.BIC.Lasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),ns=2000,f.dum=2,met="lasso")
Str.ERIC.Lasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=2000,met="lasso")


## Ridge method
Str.BIC.Ridge.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),ns=2000,f.dum=2,met="ridge")
Str.ERIC.Ridge.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=2000,met="ridge")

## Elastic Net method
Str.BIC.Enet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                    Qim=ScenI2$l2,int.mod=Strauss(12),ns=2000,f.dum=2,met="enet")
Str.ERIC.Enet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=2000,met="enet")

## Adaptive Lasso method
Str.BIC.ALasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),ns=2000,f.dum=2,met="al")
Str.ERIC.ALasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                       Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=2000,met="al")

## Adaptive Elastic Net method 
Str.BIC.AEnet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),ns=2000,f.dum=2,met="aenet")
Str.ERIC.AEnet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=2000,met="aenet")

save(Str.BIC.Lasso.ScI2,Str.ERIC.Lasso.ScI2,Str.BIC.Ridge.ScI2,Str.ERIC.Ridge.ScI2,Str.BIC.Enet.ScI2,
     Str.ERIC.Enet.ScI2,Str.BIC.ALasso.ScI2,Str.ERIC.ALasso.ScI2,Str.BIC.AEnet.ScI2,
     Str.ERIC.AEnet.ScI2, file="data/Strauss.ScenarioI2.RData")

######------------ gamma=1.5 for geyer model

## Lasso method
Gey.BIC.Lasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=2000,f.dum=2,met="lasso")
Gey.ERIC.Lasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=2000,met="lasso")

## Ridge method
Gey.BIC.Ridge.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=2000,f.dum=2,met="ridge")
Gey.ERIC.Ridge.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=2000,met="ridge")

## Elastic Net method
Gey.BIC.Enet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                    tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=2000,f.dum=2,met="enet")
Gey.ERIC.Enet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                     f.dum=2,ns=2000,met="enet")

## Adaptive Lasso method
Gey.BIC.ALasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=2000,f.dum=2,met="al")
Gey.ERIC.ALasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                       tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=2,ns=2000,met="al")

## Adaptive Elastic Net 
Gey.BIC.AEnet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=2000,f.dum=2,met="aenet")
Gey.ERIC.AEnet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=2000,met="aenet")

save(Gey.BIC.Lasso.ScI2, Gey.ERIC.Lasso.ScI2, Gey.BIC.Ridge.ScI2, Gey.ERIC.Ridge.ScI2, 
     Gey.BIC.Enet.ScI2, Gey.ERIC.Enet.ScI2, Gey.BIC.ALasso.ScI2, Gey.ERIC.ALasso.ScI2, 
     Gey.BIC.AEnet.ScI2, Gey.ERIC.AEnet.ScI2, file="data/Geyer.ScenarioI2.RData")

######---------- gamma=0.5 pour strauss model

## Lasso method
Str.BIC.Lasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=2,met="lasso")
Str.ERIC.Lasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=4000,met="lasso")

## Ridge method
Str.BIC.Ridge.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=2,met="ridge")
Str.ERIC.Ridge.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=4000,met="ridge")

## Elastic Net method
Str.BIC.Enet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                    Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=2,met="enet")
Str.ERIC.Enet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=4000,met="enet")

## Adaptive Lasso method
Str.BIC.ALasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=2,met="al")
Str.ERIC.ALasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                       Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=4000,met="al")

## Adaptive Elastic Net method 
Str.BIC.AEnet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=2,met="aenet")
Str.ERIC.AEnet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=2,ns=4000,met="aenet")

save(Str.BIC.Lasso.ScI3,Str.ERIC.Lasso.ScI3,Str.BIC.Ridge.ScI3,Str.ERIC.Ridge.ScI3,Str.BIC.Enet.ScI3,
     Str.ERIC.Enet.ScI3,Str.BIC.ALasso.ScI3,Str.ERIC.ALasso.ScI3,Str.BIC.AEnet.ScI3,
     Str.ERIC.AEnet.ScI3, file="data/Strauss.ScenarioI3.RData")

######------------ gamma=1.5 for geyer model

## Lasso method
Gey.BIC.Lasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=2,met="lasso")
Gey.ERIC.Lasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=4000,met="lasso")

## Ridge method
Gey.BIC.Ridge.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=2,met="ridge")
Gey.ERIC.Ridge.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=4000,met="ridge")

## Elastic Net method
Gey.BIC.Enet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                    tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=2,met="enet")
Gey.ERIC.Enet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                     f.dum=2,ns=4000,met="enet")

## Adaptive Lasso method
Gey.BIC.ALasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=2,met="al")
Gey.ERIC.ALasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                       tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=2,ns=4000,met="al")

## Adaptive Elastic Net 
Gey.BIC.AEnet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=2,met="aenet")
Gey.ERIC.AEnet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=2,ns=4000,met="aenet")

save(Gey.BIC.Lasso.ScI3, Gey.ERIC.Lasso.ScI3, Gey.BIC.Ridge.ScI3, Gey.ERIC.Ridge.ScI3, 
     Gey.BIC.Enet.ScI3, Gey.ERIC.Enet.ScI3, Gey.BIC.ALasso.ScI3, Gey.ERIC.ALasso.ScI3, 
     Gey.BIC.AEnet.ScI3, Gey.ERIC.AEnet.ScI3, file="data/Geyer.ScenarioI3.RData")