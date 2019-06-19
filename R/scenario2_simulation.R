# Scenario 2: Building multicolinearity
## Runner scenario 1 avant
## ---- Scenario2 ----
library(spatstat)
Sigma <- matrix(NA, nrow=50, ncol=50) 
for(i in 1:50) Sigma[i,] <- 0.7^(abs(i-(1:50)))
Sigma[1,2] <- Sigma[2,1] <- 0
V <- chol(Sigma)  
sum_vec_array <- function(Vect,Arr){
  s <- 0
  for(i in 1:length(Vect)) s <- s + Vect[i]*Arr[,,i]
  return(s)
}
Qim2 <- array(0, dim=c(101,201,50))
for(i in 1:50) { Qim2[,,i] <- sum_vec_array(V[,i],Qim1)}
Qim2.cr <- Standardize.cov(Qim2,W)
beta0.sc2 <- round ( log((4000/ integral( exp(2*Qim2.cr[[1]]+0.75*Qim2.cr[[2]]),W) ) ),4)
trend.function2 <- exp(2*Qim2.cr[[1]]+0.75*Qim2.cr[[2]])

## ---- Theta2 ----
Theta.Init.Strauss2 <- c(beta0.sc2,2,0.75,rep(0,48),log(0.5))
Theta.Init.Poisson2 <- c(beta0.sc2,2,0.75,rep(0,48),log(1))
Theta.Init.Geyer2 <- c(beta0.sc2,2,0.75,rep(0,48),log(1.5))

########################-------------------- Simulation Results for Scenario 2

######-------------------------- gamma=0.5 for strauss model

## Lasso method
Strauss.BIC.Lasso.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                        Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,Pr=.6)
Strauss.ERIC.Lasso.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                         Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,Pr=.6)

## Ridge method
Strauss.BIC.Ridge.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                        Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge",Pr=.6)
Strauss.ERIC.Ridge.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                         Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge",Pr=.6)

## Elastic Net method
Strauss.BIC.Enet.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                       Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet",Pr=.6)
Strauss.ERIC.Enet.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                        Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet",Pr=.6)

## Adaptive Lasso method
Strauss.BIC.ALasso.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                         Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="al",Pr=.6)
Strauss.ERIC.ALasso.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                          Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al",Pr=.6)

## Adaptive Elastic Net method 
Strauss.BIC.AEnet.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                        Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet",Pr=.6)
Strauss.ERIC.AEnet.Sc2 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc2,par.interact=.5,wdow=W,tr=trend.function2,
                                         Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet",Pr=.6)

save(Strauss.BIC.Lasso.Sc2,Strauss.ERIC.Lasso.Sc2,Strauss.BIC.Ridge.Sc2,Strauss.ERIC.Ridge.Sc2,Strauss.BIC.Enet.Sc2,
     Strauss.ERIC.Enet.Sc2,Strauss.BIC.ALasso.Sc2,Strauss.ERIC.ALasso.Sc2,Strauss.BIC.AEnet.Sc2,
     Strauss.ERIC.AEnet.Sc2, file="data/Strauss.Scenario2.RData")


######-------------------------- gamma=1 for Poisson model

## Lasso method
Pois.BIC.Lasso.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                     Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,Pr=.4)
Pois.ERIC.Lasso.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                      Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,Pr=.4)

## Ridge method
Pois.BIC.Ridge.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                     Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge",Pr=.4)
Pois.ERIC.Ridge.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                      Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge",Pr=.8)

## Elastic Net method
Pois.BIC.Enet.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                    Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet",Pr=.4)
Pois.ERIC.Enet.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                     Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet",Pr=.4)

## Adaptive Lasso method
Pois.BIC.ALasso.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                      Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="al",Pr=.4)
Pois.ERIC.ALasso.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                       Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al",Pr=.4)

## Adaptive Elastic Net 
Pois.BIC.AEnet.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                     Qim=Qim2,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet",Pr=.4)
Pois.ERIC.AEnet.Sc2 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc2,par.interact=1,wdow=W,tr=trend.function2,
                                      Qim=Qim2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet",Pr=.4)

save(Pois.BIC.Lasso.Sc2, Pois.ERIC.Lasso.Sc2, Pois.BIC.Ridge.Sc2, Pois.ERIC.Ridge.Sc2, 
     Pois.BIC.Enet.Sc2, Pois.ERIC.Enet.Sc2, Pois.BIC.ALasso.Sc2, Pois.ERIC.ALasso.Sc2, 
     Pois.BIC.AEnet.Sc2, Pois.ERIC.AEnet.Sc2, file="data/Poisson.Scenario2.RData")

######-------------------------- gamma=1.5 for geyer model

## Lasso method
Geyer.BIC.Lasso.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                      tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,Pr=.25)
Geyer.ERIC.Lasso.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                       tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,Pr=.25)

## Ridge method
Geyer.BIC.Ridge.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                      tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="ridge",Pr=.25)
Geyer.ERIC.Ridge.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                       tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="ridge",Pr=.25)

## Elastic Net method
Geyer.BIC.Enet.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                     tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="enet",Pr=.25)
Geyer.ERIC.Enet.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                      tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="enet",Pr=.25)

## Adaptive Lasso method
Geyer.BIC.ALasso.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                       tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="al",Pr=.25)
Geyer.ERIC.ALasso.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                        tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                        f.dum=8,ns=4000,met="al",Pr=.25)

## Adaptive Elastic Net 
Geyer.BIC.AEnet.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                      tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="aenet",Pr=.25)
Geyer.ERIC.AEnet.Sc2 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc2,par.interact=1.5,wdow=W,
                                       tr=trend.function2,Qim=Qim2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="aenet",Pr=.25)

save(Geyer.BIC.Lasso.Sc2, Geyer.ERIC.Lasso.Sc2, Geyer.BIC.Ridge.Sc2, Geyer.ERIC.Ridge.Sc2, 
     Geyer.BIC.Enet.Sc2, Geyer.ERIC.Enet.Sc2, Geyer.BIC.ALasso.Sc2, Geyer.ERIC.ALasso.Sc2, 
     Geyer.BIC.AEnet.Sc2, Geyer.ERIC.AEnet.Sc2, file="data/Geyer.Scenario2.RData")


