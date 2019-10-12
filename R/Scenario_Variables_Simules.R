## ---- ScenarioI ----
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

#####------------ (1). Fenêtre d'observation 1
## ---- ScenarioI1 ----
library(spatstat)
W1 <- owin(c(0,250),c(0,125))
ScenI1 <- built_scenario1(W=W1, 500)

## ---- ThetaI1 ----
Theta.Init.StraussI1.02 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0-2),log(0.2))
Theta.Init.StraussI1.05 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0-2),log(0.5))
Theta.Init.GeyerI1 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0-2),log(1.5))
#Theta.Init.AreaI1 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0-2),log(1.2))

## ---- SimulationI1 ----
####--------- gamma=0.2 pour strauss model
Strauss.02_W1 <-  Sim.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.2,wdow=W1,tr=ScenI1$l3[W1],
                       Qim=ScenI1$l2,int.mod=Strauss(9.25),ns=500,f.dum=16)

####--------- gamma=0.5 pour strauss model
Strauss.05_W1 <-  Sim.GPP(number.iterations=500,beta0=ScenI1$l1,par.interact=.5,wdow=W1,tr=ScenI1$l3[W1],
                          Qim=ScenI1$l2,int.mod=Strauss(9.25),ns=500,f.dum=16)

####--------- gamma=1.5 pour geyer model
Geyer_W1 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenI1$l1,par.interact=1.5,wdow=W1,tr=ScenI1$l3[W1],
                       Qim=ScenI1$l2,int.mod=Geyer(9.25,1),satur=1,ns=500)

####--------- eta=1.2 pour area model
#Area_W1 <-  Sim.GPP(Model="areaint",number.iterations=500,beta0=ScenI1$l1,par.interact=1.2,wdow=W1,tr=ScenI1$l3[W1],
                #     Qim=ScenI1$l2,int.mod=AreaInter(9.25),ns=500)

save(Strauss.02_W1, Strauss.05_W1, Geyer_W1, file="data/ScenarioI_W1.RData")

#####------------ (2). Fenêtre d'observation 2
## ---- ScenarioI2 ----
library(spatstat)
W2 <- owin(c(0,500),c(0,250)) 
ScenI2 <- built_scenario1(W=W2, 2000)

## ---- ThetaI2 ----
Theta.Init.StraussI2.02 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0-2),log(0.2))
Theta.Init.StraussI2.05 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0-2),log(0.5))
Theta.Init.GeyerI2 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0-2),log(1.5))
#Theta.Init.AreaI2 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0-2),log(1.2))

## ---- SimulationI2 ----
Strauss.02_W2 <-  Sim.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.2,wdow=W2,tr=ScenI2$l3[W2],
                       Qim=ScenI2$l2,int.mod=Strauss(9.25),ns=2000,f.dum=16)

Strauss.05_W2 <-  Sim.GPP(number.iterations=500,beta0=ScenI2$l1,par.interact=.5,wdow=W2,tr=ScenI2$l3[W2],
                          Qim=ScenI2$l2,int.mod=Strauss(9.25),ns=2000,f.dum=16)

Geyer_W2 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenI2$l1,par.interact=1.5,wdow=W2,tr=ScenI2$l3[W2],
                     Qim=ScenI2$l2,int.mod=Geyer(9.25,1),satur=1,ns=2000)

#Area_W2 <-  Sim.GPP(Model="areaint",number.iterations=500,beta0=ScenI2$l1,par.interact=1.2,wdow=W2,tr=ScenI2$l3[W2],
                 #    Qim=ScenI2$l2,int.mod=AreaInter(9.25),ns=2000)

save(Strauss.02_W2, Strauss.05_W2, Geyer_W2, file="data/ScenarioI_W2.RData")

#####------------ (3). Fenêtre d'observation 3
## ---- ScenarioI3 ----
W3 <- owin(c(0,1000),c(0,500))
ScenI3 <- built_scenario1(W=W3, 4000)

## ---- ThetaI3 ----
Theta.Init.StraussI3.02 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0-2),log(0.2))
Theta.Init.StraussI3.05 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0-2),log(0.5))
Theta.Init.GeyerI3 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0-2),log(1.5))
#Theta.Init.AreaI3 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0-2),log(1.2))

## ---- SimulationI3 ----
Strauss.02_W3 <-  Sim.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.2,wdow=W3,tr=ScenI3$l3,
                       Qim=ScenI3$l2,int.mod=Strauss(9.25),ns=4000,f.dum=16)

Strauss.05_W3 <-  Sim.GPP(number.iterations=500,beta0=ScenI3$l1,par.interact=.5,wdow=W3,tr=ScenI3$l3,
                          Qim=ScenI3$l2,int.mod=Strauss(9.25),ns=4000,f.dum=16)

Geyer_W3 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.5,wdow=W3,tr=ScenI3$l3,
                     Qim=ScenI3$l2,int.mod=Geyer(9.25,1),satur=1,ns=4000)
#Area_W3 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.2,wdow=W3,tr=ScenI3$l3,
                 #    Qim=ScenI3$l2,int.mod=AreaInter(9.25),ns=4000)

save(Strauss.02_W3, Strauss.05_W3, Geyer_W3, file="data/ScenarioI_W3.RData")



save(Strauss.02_W1, Strauss.05_W1, Geyer_W1, Strauss.02_W2, Strauss.05_W2, Geyer_W2, Strauss.02_W3, Strauss.05_W3, 
     Geyer_W3, file="data/ScenarioI.RData")








