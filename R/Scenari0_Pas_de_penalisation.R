####--------------- Scénario 0: Pas de pénalisation

##------------- (1). Fenêtre d'observation numéro 1
## ---- Scenario01 ----
library(spatstat)
W1 <- owin(c(0,250),c(0,125))
Qim01 <- array(0, dim=c(101,201,2))
Qim01[,,1] <- bei.extra$elev$v
Qim01[,,2] <- bei.extra$grad$v
Qim01.cr <- Standardize.cov(Qim01,W1)
b01 <- round ( log((500/ integral( exp(2*Qim01.cr[[1]]+0.75*Qim01.cr[[2]]),W1) ) ),4)
tf01 <- exp(2*Qim01.cr[[1]]+0.75*Qim01.cr[[2]])[W1]

## ---- Theta01 ----
Theta.Init.Strauss01.02 <- c(b01,2,0.75,log(0.2))
Theta.Init.Strauss01.05 <- c(b01,2,0.75,log(0.5))
Theta.Init.Geyer01 <- c(b01,2,0.75,log(1.5))
#Theta.Init.Area01 <- c(b01,2,0.75,log(1.2))

## ---- Simulation01 ---- 
######---------- gamma=0.2 pour strauss model
Str01.02 <- Sim.GPP.Nopen(number.iterations=500,beta0=b01,wdow=W1,tr=tf01,Qim=Qim01,int.mod=Strauss(12),ns=500,f.dum=8)

######---------- gamma=0.5 pour strauss model
Str01.05 <- Sim.GPP.Nopen(number.iterations=500,beta0=b01,par.interact=.5,wdow=W1,tr=tf01,
                       Qim=Qim01,int.mod=Strauss(9.25),ns=500,f.dum=16)

######------------ gamma=1.5 for geyer model
Gey01 <- Sim.GPP.Nopen(Model="geyer",number.iterations=500,beta0=b01,par.interact=1.5,wdow=W1,
                                     tr=tf01,Qim=Qim01,int.mod=Geyer(9.25,1),satur=1,ns=500,f.dum=2)

######------------ eta=1.2 for Area model
#Area01 <- Sim.GPP.Nopen(Model="areaint",number.iterations=2000,beta0=b01,par.interact=1.2,wdow=W1,
                     #  tr=tf01,Qim=Qim01,int.mod=AreaInter(r=9.25),ns=500,f.dum=2)

##------------- (2). Fenêtre d'observation numéro 2
## ---- Scenario02 ----
library(spatstat)
W2 <- owin(c(0,500),c(0,250))
Qim02 <- array(0, dim=c(101,201,2))
Qim02[,,1] <- bei.extra$elev$v
Qim02[,,2] <- bei.extra$grad$v
Qim02.cr <- Standardize.cov(Qim02,W2)
b02 <- round ( log((2000/ integral( exp(2*Qim02.cr[[1]]+0.75*Qim02.cr[[2]]),W2) ) ),4)
tf02 <- exp(2*Qim02.cr[[1]]+0.75*Qim02.cr[[2]])[W2]

## ---- Theta02 ----
Theta.Init.Strauss02.02 <- c(b02,2,0.75,log(0.2))
Theta.Init.Strauss02.05 <- c(b02,2,0.75,log(0.5))
Theta.Init.Geyer02 <- c(b02,2,0.75,log(1.5))
#Theta.Init.Area02 <- c(b02,2,0.75,log(1.2))

## ---- Simulation02 ----
######---------- gamma=0.5 pour strauss model
Str02.02 <- Sim.GPP.Nopen(number.iterations=500,beta0=b02,par.interact=.2,wdow=W2,tr=tf02,
                        Qim=Qim02,int.mod=Strauss(9.25),ns=2000,f.dum=16)

######---------- gamma=0.5 pour strauss model
Str02.05 <- Sim.GPP.Nopen(number.iterations=500,beta0=b02,par.interact=.5,wdow=W2,tr=tf02,
                          Qim=Qim02,int.mod=Strauss(9.25),ns=2000,f.dum=16)

######------------ gamma=1.5 for geyer model
Gey02 <- Sim.GPP.Nopen(Model="geyer",number.iterations=500,beta0=b02,par.interact=1.5,wdow=W2,
                        tr=tf02,Qim=Qim02,int.mod=Geyer(9.25,1),satur=1,ns=2000,f.dum=2)

######------------ eta=1.2 for geyer model
#Area02 <- Sim.GPP.Nopen(Model="areaint",number.iterations=2000,beta0=b02,par.interact=1.2,wdow=W2,
                 #       tr=tf02,Qim=Qim02,int.mod=AreaInter(9.25),ns=2000,f.dum=2)

##------------- (3). Fenêtre d'observation numéro 3
## ---- Scenario03 ----
library(spatstat)
W3 <- owin(c(0,1000),c(0,500))
Qim03 <- array(0, dim=c(101,201,2))
Qim03[,,1] <- bei.extra$elev$v
Qim03[,,2] <- bei.extra$grad$v
Qim03.cr <- Standardize.cov(Qim03,W3)
b03 <- round ( log((4000/ integral( exp(2*Qim03.cr[[1]]+0.75*Qim03.cr[[2]]),W3) ) ),4)
tf03 <- exp(2*Qim03.cr[[1]]+0.75*Qim03.cr[[2]])

## ---- Theta03 ----
Theta.Init.Strauss03.02 <- c(b03,2,0.75,log(0.2))
Theta.Init.Strauss03.05 <- c(b03,2,0.75,log(0.5))
Theta.Init.Geyer03 <- c(b03,2,0.75,log(1.5))
#Theta.Init.Area03 <- c(b03,2,0.75,log(1.2))


## ---- Simulation03 ----
######---------- gamma=0.5 pour strauss model
Str03.02 <- Sim.GPP.Nopen(number.iterations=500,beta0=b03,par.interact=.2,wdow=W3,tr=tf03,
                        Qim=Qim03,int.mod=Strauss(9.25),ns=4000,f.dum=16)

######---------- gamma=0.5 pour strauss model
Str03.05 <- Sim.GPP.Nopen(number.iterations=500,beta0=b03,par.interact=.5,wdow=W3,tr=tf03,
                           Qim=Qim03,int.mod=Strauss(9.25),ns=4000,f.dum=16)

######------------ gamma=1.5 for geyer model
Gey03 <- Sim.GPP.Nopen(Model="geyer",number.iterations=500,beta0=b03,par.interact=1.5,wdow=W3,
                        tr=tf03,Qim=Qim03,int.mod=Geyer(9.25,1),satur=1,ns=4000,f.dum=2)

######------------ gamma=1.5 for geyer model
#Area03 <- Sim.GPP.Nopen(Model="areaint",number.iterations=2000,beta0=b03,par.interact=1.2,wdow=W3,
                     #   tr=tf03,Qim=Qim03,int.mod=AreaInter(9.25),ns=4000,f.dum=2)

save(Str01.02, Str01.05, Gey01, Str02.02, Str02.05, Gey02, Str03.02, Str03.05, Gey03, file="data/Scenario01.RData")
