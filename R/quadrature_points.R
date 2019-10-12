## Standardize the covariates
Standardize.cov <- function(Qimage, Wdw){
  cov.list <- NULL
  Qmeans <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, mean)
  Qsd <- apply(matrix(Qimage, prod(dim(Qimage[,,1])), length(Qimage[1,1,])), 2, sd)
  for(k in 1:length(Qimage[1,1,])){
    cov.list[[k]] <- as.im((Qimage[,,k]-Qmeans[k])/Qsd[k], W=Wdw)
  }
  return(cov.list)
} 

library(spatstat)
W <- owin(c(0,1000),c(0,500))
Qim <- array(0, dim=c(101,201,50))
Qim[,,1] <- bei.extra$elev$v
Qim[,,2] <- bei.extra$grad$v
for(i in 3:50) { Qim[,,i] <- rnoise(rnorm,dimyx=c(101,201),w=W,mean=0,sd=1)$v }
Qim.cr <- Standardize.cov(Qim,W)
b0 <- round ( log((4000/ integral( exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]]),W) ) ),4) # we fix the intercept such that 
# we have 4000 points in average with the inhomogeneous Poisson point process.
trend.function <- exp(2*Qim.cr[[1]]+0.75*Qim.cr[[2]])


# Covariates list and names
cov.list <- NULL; lstnames <- NULL

# Standardize the covariates
Qmeans <- apply(matrix(Qim, prod(dim(Qim[,,1])), length(Qim[1,1,])), 2, mean)
Qsd <- apply(matrix(Qim, prod(dim(Qim[,,1])), length(Qim[1,1,])), 2, sd)
cov.list[[1]] <- as.im((Qim[,,1]-Qmeans[1])/Qsd[1], W=W)
lstnames <- 'x1'
for(k in 2:length(Qim[1,1,])){
  lstnames <- c(lstnames, paste('x', k, sep = ''))
  cov.list[[k]] <- as.im((Qim[,,k]-Qmeans[k])/Qsd[k], W=W)
}
names(cov.list) <- lstnames

# Fitting Gibbs 
mod <- list(cif="geyer", par=list(beta=exp(b0), gamma=1.5, r=12, sat=1), w=W, trend=trend.function)
X <- rmh(mod, start=list(n.start=4000), control=list(nrep=1e6), verbose=FALSE)
temp <- ppm(X ~ .  , data=cov.list, nd=70)
var.cov <- vcov(temp, fine=FALSE)
quad.ppm(temp)
