## ---- Scenario1 ----
library(spatstat)
W <- owin(c(0,1000),c(0,500))
Qim1 <- array(0, dim=c(101,201,50))
Qim1[,,1] <- bei.extra$elev$v
Qim1[,,2] <- bei.extra$grad$v
for(i in 3:50) { Qim1[,,i] <- rnoise(rnorm,dimyx=c(101,201),w=W,mean=0,sd=1)$v }
Qim1.cr <- Standardize.cov(Qim1,W)
beta0.sc1 <- round ( log((4000/ integral( exp(2*Qim1.cr[[1]]+0.75*Qim1.cr[[2]]),W) ) ),4)
trend.function1 <- exp(2*Qim1.cr[[1]]+0.75*Qim1.cr[[2]])

## ---- Theta1 ----
Theta.Init.Strauss1 <- c(beta0.sc1,2,0.75,rep(0,48),log(0.5))
Theta.Init.Poisson1 <- c(beta0.sc1,2,0.75,rep(0,48),log(1))
Theta.Init.Geyer1 <- c(beta0.sc1,2,0.75,rep(0,48),log(1.5))

## ---- Simulation Results for Scenario 1 ----

######-------------------------- gamma=0.5 for strauss model

## Lasso method
Strauss.BIC.Lasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                        Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8)
Strauss.ERIC.Lasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                         Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000)

## Ridge method
Strauss.BIC.Ridge.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                        Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge")
Strauss.ERIC.Ridge.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                         Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Strauss.BIC.Enet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                       Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet")
Strauss.ERIC.Enet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                        Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Strauss.BIC.ALasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                         Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="al")
Strauss.ERIC.ALasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                          Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net method 
Strauss.BIC.AEnet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                        Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet")
Strauss.ERIC.AEnet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=.5,wdow=W,tr=trend.function1,
                                         Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet")

save(Strauss.BIC.Lasso.Sc1,Strauss.ERIC.Lasso.Sc1,Strauss.BIC.Ridge.Sc1,Strauss.ERIC.Ridge.Sc1,Strauss.BIC.Enet.Sc1,
     Strauss.ERIC.Enet.Sc1,Strauss.BIC.ALasso.Sc1,Strauss.ERIC.ALasso.Sc1,Strauss.BIC.AEnet.Sc1,
     Strauss.ERIC.AEnet.Sc1, file="data/Strauss.Scenario1.RData")

######-------------------------- gamma=1 for Poisson model

## Lasso method
Pois.BIC.Lasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                     Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8)
Pois.ERIC.Lasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                      Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000)

## Ridge method
Pois.BIC.Ridge.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                     Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge")
Pois.ERIC.Ridge.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                      Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Pois.BIC.Enet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                    Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet")
Pois.ERIC.Enet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                     Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Pois.BIC.ALasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                      Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="al")
Pois.ERIC.ALasso.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                       Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net 
Pois.BIC.AEnet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                     Qim=Qim1,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet")
Pois.ERIC.AEnet.Sc1 <- Simulation.GPP(number.iterations=500,beta0=beta0.sc1,par.interact=1,wdow=W,tr=trend.function1,
                                      Qim=Qim1,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet")

save(Pois.BIC.Lasso.Sc1, Pois.ERIC.Lasso.Sc1, Pois.BIC.Ridge.Sc1, Pois.ERIC.Ridge.Sc1, 
     Pois.BIC.Enet.Sc1, Pois.ERIC.Enet.Sc1,Pois.BIC.ALasso.Sc1, Pois.ERIC.ALasso.Sc1, 
     Pois.BIC.AEnet.Sc1, Pois.ERIC.AEnet.Sc1, file="data/Poisson.Scenario1.RData")

######-------------------------- gamma=1.5 for geyer model

## Lasso method
Geyer.BIC.Lasso.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                      tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8)
Geyer.ERIC.Lasso.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                       tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000)

## Ridge method
Geyer.BIC.Ridge.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                      tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="ridge",Pr=.25)
Geyer.ERIC.Ridge.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                       tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Geyer.BIC.Enet.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                     tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="enet")
Geyer.ERIC.Enet.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                      tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Geyer.BIC.ALasso.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                       tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="al")
Geyer.ERIC.ALasso.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                        tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                        f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net 
Geyer.BIC.AEnet.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                      tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="aenet")
Geyer.ERIC.AEnet.Sc1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=beta0.sc1,par.interact=1.5,wdow=W,
                                       tr=trend.function1,Qim=Qim1,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="aenet")

save(Geyer.BIC.Lasso.Sc1, Geyer.ERIC.Lasso.Sc1, Geyer.BIC.Ridge.Sc1, Geyer.ERIC.Ridge.Sc1, 
     Geyer.BIC.Enet.Sc1, Geyer.ERIC.Enet.Sc1,Geyer.BIC.ALasso.Sc1, Geyer.ERIC.ALasso.Sc1, 
     Geyer.BIC.AEnet.Sc1, Geyer.ERIC.AEnet.Sc1, file="data/Geyer.Scenario1.RData")

