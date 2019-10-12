W <- owin(c(0,1000),c(0,500))
Qim <- array(0, dim=c(101,201,2))
Qim[,,1] <- bei.extra$elev$v
Qim[,,2] <- bei.extra$grad$v
Qim.cr <- Standardize.cov(Qim,W)
beta0 <- round ( log((4000/ integral( exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]]),W) ) ),4)
tr.f <- exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]])
mod <- list(cif="strauss", par=list(beta=exp(beta0), gamma=.5, r=12), w=W, trend=tr.f)
X <- rmh(mod, start=list(n.start=4000), control=list(nrep=1e6), verbose=FALSE)
plot(X, main="")
 
m <- 500
nbr_pts <- rep(NA,m)
for(i in 1:m){
  X <- rmh(mod, start=list(n.start=4000), control=list(nrep=1e6), verbose=FALSE)
  nbr_pts[i] <- npoints(X)
} 

nbr_en_moy <- mean(nbr_pts)


mod1 <- list(cif="areaint", par=list(beta=exp(beta0), eta=2, r=12), w=W, trend=tr.f)
X1 <- rmh(mod, start=list(n.start=4000), control=list(nrep=1e6), verbose=FALSE)
mod2 <- list(cif="areaint", par=list(beta=exp(beta0), eta=0.02, r=12), w=W, trend=tr.f)
X2 <- rmh(mod, start=list(n.start=4000), control=list(nrep=1e6), verbose=FALSE)

