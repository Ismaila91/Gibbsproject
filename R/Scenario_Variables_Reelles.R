## ---- ScenarioII ----
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

#####------------ (1). Fenêtre d'observation 1
## ---- ScenarioII1 ----
library(spatstat)
W1 <- owin(c(0,250),c(0,125))
ScenII1 <- built_scenario2(W=W1, 500)

## ---- ThetaII1
Theta.Init.StraussII1.02 <- c(ScenII1$l1,2,0.75,rep(0,ScenII1$l0-2),log(0.2))
Theta.Init.StraussII1.05 <- c(ScenII1$l1,2,0.75,rep(0,ScenII1$l0-2),log(0.5))
Theta.Init.GeyerII1 <- c(ScenII1$l1,2,0.75,rep(0,ScenII1$l0-2),log(1.5))
#Theta.Init.AreaI1 <- c(ScenI1$l1,2,0.75,rep(0,ScenI1$l0-2),log(1.2))

## ---- SimulationII1 ----
####--------- gamma=0.2 pour strauss model
StrII.02_W1 <-  Sim.GPP(number.iterations=500,beta0=ScenII1$l1,par.interact=.2,wdow=W1,tr=ScenII1$l3[W1],
                          Qim=ScenII1$l2,int.mod=Strauss(9.25),ns=500,f.dum=16)

####--------- gamma=0.5 pour strauss model
StrII.05_W1 <-  Sim.GPP(number.iterations=500,beta0=ScenII1$l1,par.interact=.5,wdow=W1,tr=ScenII1$l3[W1],
                          Qim=ScenII1$l2,int.mod=Strauss(9.25),ns=500,f.dum=16)

####--------- gamma=1.5 pour geyer model
GeyII_W1 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenII1$l1,par.interact=1.5,wdow=W1,tr=ScenII1$l3[W1],
                     Qim=ScenII1$l2,int.mod=Geyer(9.25,1),satur=1,ns=500)

####--------- eta=1.2 pour area model
#Area_W1 <-  Sim.GPP(Model="areaint",number.iterations=500,beta0=ScenI1$l1,par.interact=1.2,wdow=W1,tr=ScenI1$l3[W1],
#     Qim=ScenI1$l2,int.mod=AreaInter(9.25),ns=500)

save(StrII.02_W1, StrII.05_W1, GeyII_W1, file="data/ScenarioII_W1.RData")

#####------------ (2). Fenêtre d'observation 2
## ---- ScenarioII2 ----
library(spatstat)
W2 <- owin(c(0,500),c(0,250)) 
ScenII2 <- built_scenario2(W=W2, 2000)

## ---- ThetaII2 ----
Theta.Init.StraussII2.02 <- c(ScenII2$l1,2,0.75,rep(0,ScenII2$l0-2),log(0.2))
Theta.Init.StraussII2.05 <- c(ScenII2$l1,2,0.75,rep(0,ScenII2$l0-2),log(0.5))
Theta.Init.GeyerII2 <- c(ScenII2$l1,2,0.75,rep(0,ScenII2$l0-2),log(1.5))
#Theta.Init.AreaI2 <- c(ScenI2$l1,2,0.75,rep(0,ScenI2$l0-2),log(1.2))

## ---- SimulationII2 ----
StrII.02_W2 <-  Sim.GPP(number.iterations=500,beta0=ScenII2$l1,par.interact=.2,wdow=W2,tr=ScenII2$l3[W2],
                          Qim=ScenII2$l2,int.mod=Strauss(9.25),ns=2000,f.dum=16)

StrII.05_W2 <-  Sim.GPP(number.iterations=500,beta0=ScenII2$l1,par.interact=.5,wdow=W2,tr=ScenII2$l3[W2],
                          Qim=ScenII2$l2,int.mod=Strauss(9.25),ns=2000,f.dum=16)

GeyII_W2 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenII2$l1,par.interact=1.5,wdow=W2,tr=ScenII2$l3[W2],
                     Qim=ScenII2$l2,int.mod=Geyer(9.25,1),satur=1,ns=2000)

#Area_W2 <-  Sim.GPP(Model="areaint",number.iterations=500,beta0=ScenI2$l1,par.interact=1.2,wdow=W2,tr=ScenI2$l3[W2],
#    Qim=ScenI2$l2,int.mod=AreaInter(9.25),ns=2000)

save(StrII.02_W2, StrII.05_W2, GeyII_W2, file="data/ScenarioII_W2.RData")


#####------------ (3). Fenêtre d'observation 3
## ---- ScenarioII3 ----
W3 <- owin(c(0,1000),c(0,500))
ScenII3 <- built_scenario2(W=W3, 4000)

## ----- ThetaII3 ----
Theta.Init.StraussII3.02 <- c(ScenII3$l1,2,0.75,rep(0,ScenII3$l0-2),log(0.2))
Theta.Init.StraussII3.05 <- c(ScenII3$l1,2,0.75,rep(0,ScenII3$l0-2),log(0.5))
Theta.Init.GeyerII3 <- c(ScenII3$l1,2,0.75,rep(0,ScenII3$l0-2),log(1.5))
#Theta.Init.AreaI3 <- c(ScenI3$l1,2,0.75,rep(0,ScenI3$l0-2),log(1.2))

## ---- SimulationI3 ----
StrII.02_W3 <-  Sim.GPP(number.iterations=500,beta0=ScenII3$l1,par.interact=.2,wdow=W3,tr=ScenII3$l3,
                          Qim=ScenII3$l2,int.mod=Strauss(9.25),ns=4000,f.dum=16)

StrII.05_W3 <-  Sim.GPP(number.iterations=500,beta0=ScenII3$l1,par.interact=.5,wdow=W3,tr=ScenII3$l3,
                          Qim=ScenII3$l2,int.mod=Strauss(9.25),ns=4000,f.dum=16)

GeyII_W3 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenII3$l1,par.interact=1.5,wdow=W3,tr=ScenII3$l3,
                     Qim=ScenII3$l2,int.mod=Geyer(9.25,1),satur=1,ns=4000)
#Area_W3 <-  Sim.GPP(Model="geyer",number.iterations=500,beta0=ScenI3$l1,par.interact=1.2,wdow=W3,tr=ScenI3$l3,
#    Qim=ScenI3$l2,int.mod=AreaInter(9.25),ns=4000)

save(StrII.02_W3, StrII.05_W3, GeyII_W3, file="data/ScenarioII_W3.RData")


save(StrII.02_W1, StrII.05_W1, GeyII_W1, StrII.02_W2, StrII.05_W2, GeyII_W2, StrII.02_W3, StrII.05_W3, 
     GeyII_W3, file="data/ScenarioII.RData")


