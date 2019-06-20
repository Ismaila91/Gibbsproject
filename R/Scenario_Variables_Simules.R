library(spatstat)
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

####-------------- Scenario I: Variables simulées

#####------------ (1). Fenêtre d'observation 1
W1 <- owin(c(0,250),c(0,125))
ScenI1 <- built_scenario1(W=W1, 500)
Theta.Init.StraussI1 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0),log(0.5))
Theta.Init.GeyerI1 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0),log(1.5))

######---------- gamma=0.5 pour strauss model

## Lasso method
Str.BIC.Lasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                        Qim=ScenI1$l2,int.mod=Strauss(12),ns=4000,f.dum=8)
Str.ERIC.Lasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                         Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000)

## Ridge method
Str.BIC.Ridge.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                        Qim=ScenI1$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge")
Str.ERIC.Ridge.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                         Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Str.BIC.Enet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                       Qim=ScenI1$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet")
Str.ERIC.Enet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                        Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Str.BIC.ALasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                         Qim=ScenI1$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="al")
Str.ERIC.ALasso.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                          Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net method 
Str.BIC.AEnet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                        Qim=ScenI1$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet")
Str.ERIC.AEnet.ScI1 <- Simulation.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3,
                                         Qim=ScenI1$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet")

save(Str.BIC.Lasso.ScI1,Str.ERIC.Lasso.ScI1,Str.BIC.Ridge.ScI1,Str.ERIC.Ridge.ScI1,Str.BIC.Enet.ScI1,
     Str.ERIC.Enet.ScI1,Str.BIC.ALasso.ScI1,Str.ERIC.ALasso.ScI1,Str.BIC.AEnet.ScI1,
     Str.ERIC.AEnet.ScI1, file="data/Strauss.ScenarioI1.RData")

######------------ gamma=1.5 for geyer model

## Lasso method
Gey.BIC.Lasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8)
Gey.ERIC.Lasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                       tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000)

## Ridge method
Gey.BIC.Ridge.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="ridge")
Gey.ERIC.Ridge.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                       tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Gey.BIC.Enet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                     tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="enet")
Gey.ERIC.Enet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Gey.BIC.ALasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                       tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="al")
Gey.ERIC.ALasso.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                        tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                        f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net 
Gey.BIC.AEnet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                      tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="aenet")
Gey.ERIC.AEnet.ScI1 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,
                                       tr=ScenI1$l3,Qim=ScenI1$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="aenet")

save(Gey.BIC.Lasso.ScI1, Gey.ERIC.Lasso.ScI1, Gey.BIC.Ridge.ScI1, Gey.ERIC.Ridge.ScI1, 
     Gey.BIC.Enet.ScI1, Gey.ERIC.Enet.ScI1, Gey.BIC.ALasso.ScI1, Gey.ERIC.ALasso.ScI1, 
     Gey.BIC.AEnet.ScI1, Gey.ERIC.AEnet.ScI1, file="data/Geyer.ScenarioI1.RData")



#####------------ (2). Fenêtre d'observation 2
W2 <- owin(c(0,500),c(0,250)) 
ScenI2 <- built_scenario1(W=W2, 2000) 
Theta.Init.StraussI2 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0),log(0.5))
Theta.Init.GeyerI2 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0),log(1.5))

######---------- gamma=0.5 pour strauss model

## Lasso method
Str.BIC.Lasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),ns=4000,f.dum=8)
Str.ERIC.Lasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000)

## Ridge method
Str.BIC.Ridge.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge")
Str.ERIC.Ridge.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Str.BIC.Enet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                    Qim=ScenI2$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet")
Str.ERIC.Enet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Str.BIC.ALasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="al")
Str.ERIC.ALasso.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                       Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net method 
Str.BIC.AEnet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                     Qim=ScenI2$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet")
Str.ERIC.AEnet.ScI2 <- Simulation.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3,
                                      Qim=ScenI2$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet")

save(Str.BIC.Lasso.ScI2,Str.ERIC.Lasso.ScI2,Str.BIC.Ridge.ScI2,Str.ERIC.Ridge.ScI2,Str.BIC.Enet.ScI2,
     Str.ERIC.Enet.ScI2,Str.BIC.ALasso.ScI2,Str.ERIC.ALasso.ScI2,Str.BIC.AEnet.ScI2,
     Str.ERIC.AEnet.ScI2, file="data/Strauss.ScenarioI2.RData")

######------------ gamma=1.5 for geyer model

## Lasso method
Gey.BIC.Lasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8)
Gey.ERIC.Lasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000)

## Ridge method
Gey.BIC.Ridge.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="ridge")
Gey.ERIC.Ridge.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Gey.BIC.Enet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                    tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="enet")
Gey.ERIC.Enet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                     f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Gey.BIC.ALasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="al")
Gey.ERIC.ALasso.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                       tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net 
Gey.BIC.AEnet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                     tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="aenet")
Gey.ERIC.AEnet.ScI2 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,
                                      tr=ScenI2$l3,Qim=ScenI2$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="aenet")

save(Gey.BIC.Lasso.ScI2, Gey.ERIC.Lasso.ScI2, Gey.BIC.Ridge.ScI2, Gey.ERIC.Ridge.ScI2, 
     Gey.BIC.Enet.ScI2, Gey.ERIC.Enet.ScI2, Gey.BIC.ALasso.ScI2, Gey.ERIC.ALasso.ScI2, 
     Gey.BIC.AEnet.ScI2, Gey.ERIC.AEnet.ScI2, file="data/Geyer.ScenarioI2.RData")


#####------------ (3). Fenêtre d'observation 3
W3 <- owin(c(0,1000),c(0,500))
ScenI3 <- built_scenario1(W=W3, 4000)
Theta.Init.StraussI2 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0),log(0.5))
Theta.Init.GeyerI2 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0),log(1.5))

######---------- gamma=0.5 pour strauss model

## Lasso method
Str.BIC.Lasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=8)
Str.ERIC.Lasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000)

## Ridge method
Str.BIC.Ridge.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="ridge")
Str.ERIC.Ridge.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Str.BIC.Enet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                    Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="enet")
Str.ERIC.Enet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Str.BIC.ALasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="al")
Str.ERIC.ALasso.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                       Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net method 
Str.BIC.AEnet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                     Qim=ScenI3$l2,int.mod=Strauss(12),ns=4000,f.dum=8,met="aenet")
Str.ERIC.AEnet.ScI3 <- Simulation.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                                      Qim=ScenI3$l2,int.mod=Strauss(12),inf.criteria="ERIC",f.dum=8,ns=4000,met="aenet")

save(Str.BIC.Lasso.ScI3,Str.ERIC.Lasso.ScI3,Str.BIC.Ridge.ScI3,Str.ERIC.Ridge.ScI3,Str.BIC.Enet.ScI3,
     Str.ERIC.Enet.ScI3,Str.BIC.ALasso.ScI3,Str.ERIC.ALasso.ScI3,Str.BIC.AEnet.ScI3,
     Str.ERIC.AEnet.ScI3, file="data/Strauss.ScenarioI3.RData")

######------------ gamma=1.5 for geyer model

## Lasso method
Gey.BIC.Lasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8)
Gey.ERIC.Lasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000)

## Ridge method
Gey.BIC.Ridge.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="ridge")
Gey.ERIC.Ridge.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="ridge")

## Elastic Net method
Gey.BIC.Enet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                    tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="enet")
Gey.ERIC.Enet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                     f.dum=8,ns=4000,met="enet")

## Adaptive Lasso method
Gey.BIC.ALasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="al")
Gey.ERIC.ALasso.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                       tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                       f.dum=8,ns=4000,met="al")

## Adaptive Elastic Net 
Gey.BIC.AEnet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                     tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,ns=4000,f.dum=8,met="aenet")
Gey.ERIC.AEnet.ScI3 <- Simulation.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,
                                      tr=ScenI3$l3,Qim=ScenI3$l2,int.mod=Geyer(12,1),satur=1,inf.criteria="ERIC",
                                      f.dum=8,ns=4000,met="aenet")

save(Gey.BIC.Lasso.ScI3, Gey.ERIC.Lasso.ScI3, Gey.BIC.Ridge.ScI3, Gey.ERIC.Ridge.ScI3, 
     Gey.BIC.Enet.ScI3, Gey.ERIC.Enet.ScI3, Gey.BIC.ALasso.ScI3, Gey.ERIC.ALasso.ScI3, 
     Gey.BIC.AEnet.ScI3, Gey.ERIC.AEnet.ScI3, file="data/Geyer.ScenarioI3.RData")

