rm(list=ls())
set.seed(123)
par(mfrow=c(1,1))
pacman::p_load(pracma, data.table, wpp2015, bayesTFR, knitr, bayesPop, dplyr,
               ggplot2, forecast, bayesLife)

# Question 1

## Pull data sources from 
data(e0F)
head(e0F)
hond_e0 <- c(as.matrix(subset(e0F, country == "Honduras")[,3:(ncol(e0F) - 1)]))
years <- seq(1950, 2010, 5)
N <- length(hond_e0)
plot(hond_e0)
plot(hond_e0[1:(N-1)], diff(hond_e0))

e0_proj <- function(params, e0=hond_e0){
    a <- exp(params[1:4])
    k <- exp(params[5])
    z <- exp(params[6])
    dl1 <- k/(1 + exp((-4.4)/a[2] * (e0 - a[1] - .5*a[2])))
    dl2 <- (z-k)/(1 + exp((-4.4)/a[4] * (e0 - sum(a[1:3] - .5 * a[4]))))
    return(dl1 + dl2)
}

objfunc1 <- function(params, e0=hond_e0){
    N <- length(e0)
    sigma <- exp(params[7])
    preds <- e0 + e0_proj(params, e0)
    return(sum(-1 * dnorm(e0[2:N], preds[1:(N-1)], sigma, TRUE)))
}

plot(hond_e0[1:(N-1)], diff(hond_e0))
vals <- log(c(1.08e-07, 62.77, 20.83, 9.83, 4.37, .87, .20))
preds <- e0_proj(vals, hond_e0)[1:(N-1)]
lines(hond_e0[1:(N-1)], preds, col="red")

Obj <- nlminb(vals, objfunc1)
Obj$message
preds <- e0_proj(Obj$par, hond_e0)

jpeg('~/Documents/Classes/statdemog/week7/fitteddle0.jpg')
plot(hond_e0[1:(N-1)], diff(hond_e0), xlab="Life Expectancy", 
     ylab="Change in Life Expectancy", main="Honduras Life Expectancy Females")
lines(hond_e0[1:(N-1)], preds[1:(N-1)], col="blue")
dev.off()

mu2015 <- (hond_e0 + preds)[length(preds)]
eps <- exp(Obj$par[7])
c(`2.5%`=mu2015 - 1.96*eps, mean=mu2015, `97.5%`=mu2015 + 1.96*eps)

x <- seq(mu2015 - 4*eps, 
         mu2015 + 4*eps, length=100)
hx <- dnorm(x, mean=mu2015, sd=eps)
jpeg('~/Documents/Classes/statdemog/week7/e0honddens.jpg')
plot(x, hx, "l", xlab="2015 Honduras Life Expectancy Forecast", ylab="Density")
dev.off()

# question 2
e0dir <- '~/Downloads/e0/sim03092016'
e0mcmc <- get.e0.mcmc(e0dir, burnin=1000)
preds <- get.e0.prediction(e0mcmc)

jpeg('~/Documents/Classes/statdemog/week7/bayesdlcompare.jpg')
par(mfrow=c(1,2))
e0.DLcurve.plot(e0mcmc, "Japan", pi=95)
e0.DLcurve.plot(e0mcmc, "Cambodia", pi=95)
par(mfrow=c(1,1))
dev.off()

jpeg('~/Documents/Classes/statdemog/week7/bayesparamcomparej.jpg')
par(mfrow=c(1,2))
e0.partraces.cs.plot("Japan", e0mcmc, par.names=c("k.c", "z.c"))
dev.off()
jpeg('~/Documents/Classes/statdemog/week7/bayesparamcomparec.jpg')
e0.partraces.cs.plot("Cambodia", e0mcmc, par.names=c("k.c", "z.c"))
par(mfrow=c(1,1))
dev.off()

pJ <- get.e0.parameter.traces.cs(e0mcmc$mcmc.list, burnin=1000,
                                 get.country.object("Japan", meta=e0mcmc$meta))
summary(pJ[,"k.c_c392"])
summary(pJ[,"z.c_c392"])

pC <- get.e0.parameter.traces.cs(e0mcmc$mcmc.list, burnin=1000,
                                 get.country.object("Cambodia", meta=e0mcmc$meta))
summary(pC[,"k.c_c116"])
summary(pC[,"z.c_c116"])


jpeg('~/Documents/Classes/statdemog/week7/bayese0compare.jpg')
par(mfrow=c(1,2))
e0.trajectories.plot(preds, "Japan")
e0.trajectories.plot(preds, "Cambodia")
par(mfrow=c(1,1))
dev.off()

countries <- c("Japan", "China", "Russian Federation", "Republic of Korea", 
               "Dem. People's Rep. of Korea")
Nsims <- ncol(get.e0.trajectories(preds, "Japan"))
best <- matrix(NA,nrow=0, ncol=nrow(get.e0.trajectories(preds, "Japan")))
worst <- matrix(NA,nrow=0, ncol=nrow(get.e0.trajectories(preds, "Japan")))
for(i in 1:Nsims){
    e0sims <- sapply(countries, function(x) c(get.e0.trajectories(preds, x)[,i]))
    best <- rbind(best, apply(e0sims[,1] > e0sims[,2:length(countries)], 1, all))
    worst <- rbind(worst, apply(e0sims[,1] < e0sims[,2:length(countries)], 1, all))
}

DT <- data.table(best=apply(best, 2, mean)[-1], worst=apply(worst, 2, mean)[-1], 
                 year=seq(2015, 2095, 5))
kable(DT, "markdown")

# question 3
e0dir <- '~/Downloads/e0/sim03092016'
tfrdir <- '~/Downloads/TFR/sim03092016/'
popdir <- '~/Downloads/pop/'

# pop.pred <- pop.predict(
#     end.year=2100, start.year=1950, present.year=2015,
#     wpp.year=2015, output.dir=popdir, nr.traj=50,
#     inputs=list(tfr.sim.dir=tfrdir,
#                 e0F.sim.dir=e0dir, 
#                 e0M.sim.dir="joint_"),
#     keep.vital.events=TRUE)

pop.pred <- get.pop.prediction(popdir)

par(mfrow=c(1,1))

jpeg('~/Documents/Classes/statdemog/week7/CNpopproj.jpg')
pop.trajectories.plot(pop.pred, country="Canada", sex="both", sum.over.ages=T)
dev.off()
jpeg('~/Documents/Classes/statdemog/week7/CNpopprojmale.jpg')
pop.trajectories.plot(pop.pred, country="Canada", sex="male", sum.over.ages=T)
dev.off()
jpeg('~/Documents/Classes/statdemog/week7/CNpopproj65plus.jpg')
pop.trajectories.plot(pop.pred, country="Canada", sex="both", age=14:27,  sum.over.ages=T)
dev.off()
jpeg('~/Documents/Classes/statdemog/week7/CNpopprojpsrm.jpg')
pop.trajectories.plot(pop.pred, country="Canada", sex="male",  age="psr")
dev.off()
jpeg('~/Documents/Classes/statdemog/week7/CNpopprojpsrf.jpg')
pop.trajectories.plot(pop.pred, country="Canada", sex="female",  age="psr")
dev.off()

nacode <- 905
pop.aggr <- pop.aggregate(pop.pred, regions=nacode, verbose=TRUE)
jpeg('~/Documents/Classes/statdemog/week7/NApopproj.jpg')
pop.trajectories.plot(pop.aggr, nacode, sum.over.ages=TRUE)
dev.off()
