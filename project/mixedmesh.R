rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, TMB, ar.matrix)
set.seed(12345)

plot_mesh_sim <- function(x, proj){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8))
}


n <- 500 # number of observations on the grid
m <- 12 # number of time points
loc <- matrix(runif(n*2), n, 2) # simulate observed points
mesh <- inla.mesh.create(loc, refine=list(max.edge=0.05)) # create mesh
par(mfrow=c(1,1))
plot(mesh)
points(loc[,1], loc[,2], col="red", pch=20)

proj <- inla.mesh.projector(mesh)

sigma0 <-  .3   # Standard deviation
range0 <- 1. # Spatial range
kappa0 <- sqrt(8)/range0 # inla paramter transform
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) # inla parameter transform
rho <- .91 # tenporal autocorrelation

spde <- inla.spde2.matern(mesh) # create the spde from the mesh

# parameterize the spde and get the projected precision matrix
Q <- inla.spde2.precision(spde, theta=c(log(tau0), log(kappa0)))
x <- c(sim.AR(1, Q))



M <- length(proj$x)
DT <- data.table(x=signif(rep(proj$x, M), 4), y=signif(rep(proj$y, each=M), 4), 
                 obs=c(inla.mesh.project(proj, field=x)))
sampleN <- 15
ranges <- quantile(DT$obs, seq(0, 1, length.out=sampleN), na.rm=T)
replength <- nrow(DT[x > .33 & x < .66 & y > .33 & y < .66,])
DT[x > .33 & x < .66 & y > .33 & y < .66, obs:=ranges[rbinom(replength, sampleN-1, .5) + 1]]
ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
    lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8))
