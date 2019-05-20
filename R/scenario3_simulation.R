## ---- Scenario3 ----
library(spatstat)
library(readxl)
## If not installed, we need the here package
soil_nut <- read_xls(here::here("data","bci.block20.data.xls"),sheet=2)
soil_nut_data <- as.data.frame(soil_nut)
soil_nut_im <- as.im(soil_nut_data)
W <- owin(c(0,1000),c(0,500))
Qim3 <- array(0, dim=c(101,201,15)) 
Qim3[,,1] <- bei.extra$elev$v
Qim3[,,2] <- bei.extra$grad$v
for(i in 3:15) Qim3[,,i][1:25,1:50] <- soil_nut_im[[(i-2)]]$v   
Qim3.cr <- Standardize.cov(Qim3,W)
beta0.sc3 <- round ( log((4000/ integral( exp(2*Qim3.cr[[1]]+0.75*Qim3.cr[[2]]),W) ) ),4)
trend.function3 <- exp(2*Qim3.cr[[1]]+0.75*Qim3.cr[[2]])

## ---- Theta3 ----
Theta.Init.Strauss3 <- c(beta0.sc3,2,0.75,rep(0,13),log(0.5))
Theta.Init.Poisson3 <- c(beta0.sc3,2,0.75,rep(0,13),log(1))
Theta.Init.Geyer3 <- c(beta0.sc3,2,0.75,rep(0,13),log(1.5))

########################-------------------- Simulation Results for Scenario 3

######-------------------------- gamma=0.5 for strauss model

## Lasso method
Strauss.BIC.Lasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                        Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,Pr=.6)
Strauss.ERIC.Lasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                         Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,Pr=.6)

## Ridge method
Strauss.BIC.Ridge.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                        Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge",Pr=.6)
Strauss.ERIC.Ridge.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                         Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge",Pr=.6)

## Elastic Net method
Strauss.BIC.Enet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                       Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet",Pr=.6)
Strauss.ERIC.Enet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                        Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet",Pr=.6)

## Adaptive Lasso method
Strauss.BIC.ALasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                         Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="al",Pr=.6)
Strauss.ERIC.ALasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                          Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al",Pr=.6)

## Adaptive Elastic Net method 
Strauss.BIC.AEnet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                        Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet",Pr=.6)
Strauss.ERIC.AEnet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=.5,wdow=W,tr=trend.function3,
                                         Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet",Pr=.6)

save(Strauss.BIC.Lasso.Sc3,Strauss.ERIC.Lasso.Sc3,Strauss.BIC.Ridge.Sc3,Strauss.ERIC.Ridge.Sc3,Strauss.BIC.Enet.Sc3,
     Strauss.ERIC.Enet.Sc3,Strauss.BIC.ALasso.Sc3,Strauss.ERIC.ALasso.Sc3,Strauss.BIC.AEnet.Sc3,
     Strauss.ERIC.AEnet.Sc3, file="data/Strauss.Scenario3.RData")

######-------------------------- gamma=1 for Poisson model

## Lasso method
Pois.BIC.Lasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                     Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,Pr=.4)
Pois.ERIC.Lasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                      Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,Pr=.4)

## Ridge method
Pois.BIC.Ridge.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                     Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge",Pr=.4)
Pois.ERIC.Ridge.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                      Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge",Pr=.8)

## Elastic Net method
Pois.BIC.Enet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                    Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet",Pr=.4)
Pois.ERIC.Enet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                     Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet",Pr=.4)

## Adaptive Lasso method
Pois.BIC.ALasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                      Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="al",Pr=.4)
Pois.ERIC.ALasso.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                       Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al",Pr=.4)

## Adaptive Elastic Net 
Pois.BIC.AEnet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                     Qim=Qim3,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet",Pr=.4)
Pois.ERIC.AEnet.Sc3 <- Simulation.GPP(number.iterations=100,beta0=beta0.sc3,par.interact=1,wdow=W,tr=trend.function3,
                                      Qim=Qim3,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet",Pr=.4)

save(Pois.BIC.Lasso.Sc3, Pois.ERIC.Lasso.Sc3, Pois.BIC.Ridge.Sc3, Pois.ERIC.Ridge.Sc3, 
     Pois.BIC.Enet.Sc3, Pois.ERIC.Enet.Sc3, Pois.BIC.ALasso.Sc3, Pois.ERIC.ALasso.Sc3, 
     Pois.BIC.AEnet.Sc3, Pois.ERIC.AEnet.Sc3, file="data/Poisson.Scenario3.RData")

######-------------------------- gamma=1.5 for geyer model

## Lasso method
Geyer.BIC.Lasso.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                      tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,Pr=.25)
Geyer.ERIC.Lasso.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                       tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,Pr=.25)

## Ridge method
Geyer.BIC.Ridge.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                      tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="ridge",Pr=.25)
Geyer.ERIC.Ridge.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                       tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="ridge",Pr=.25)

## Elastic Net method
Geyer.BIC.Enet.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                     tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="enet",Pr=.25)
Geyer.ERIC.Enet.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                      tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="enet",Pr=.25)

## Adaptive Lasso method
Geyer.BIC.ALasso.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                       tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="al",Pr=.25)
Geyer.ERIC.ALasso.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                        tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                        f.dum=8,ns=4000,met="al",Pr=.25)

## Adaptive Elastic Net 
Geyer.BIC.AEnet.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                      tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="aenet",Pr=.25)
Geyer.ERIC.AEnet.Sc3 <- Simulation.GPP(Model="geyer",number.iterations=100,beta0=beta0.sc3,par.interact=1.5,wdow=W,
                                       tr=trend.function3,Qim=Qim3,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="aenet",Pr=.25)

save(Geyer.BIC.Lasso.Sc3, Geyer.ERIC.Lasso.Sc3, Geyer.BIC.Ridge.Sc3, Geyer.ERIC.Ridge.Sc3, 
     Geyer.BIC.Enet.Sc3, Geyer.ERIC.Enet.Sc3, Geyer.BIC.ALasso.Sc3, Geyer.ERIC.ALasso.Sc3, 
     Geyer.BIC.AEnet.Sc3, Geyer.ERIC.AEnet.Sc3, file="data/Geyer.Scenario3.RData")

