# Homework Week 7
### Neal Marquez

1. Obtain the values of life expectancy at birth for Honduras for 1950-2015 from the 2015
World Population Prospects and write them out in your homework.

a)  Fit a version of the six-parameter double logistic gain model by nonlinear least
squares to the gains in life expectancy, assuming that the error variance remains
constant over time (this will just give one set of double logistic parameter values
and an estimated error variance).

The double logistic model uses two logistic equations in order to capture the change between time periods of measured life expectancy and anticipates a period of high increase and then stabilizes at a lower value. The model is fit using the `nlminb` optimizer in `R` in order to fit the 6 parameters of the model as well as the error variance term. The the fit is appropriate for the data despite the model having correlated parameters and being unable to converge. Note that the functional form of the model is as follows

$$
g(\mathcal{l}, \theta) \sim \frac{z-k}{1 + exp(\frac{-2ln(9)}{\Delta_4}(\mathcal{l} - \sum_{i=1}^{3} \Delta_i - .5\Delta_4))} + \frac{k}{1 + exp(\frac{-2ln(9)}{\Delta_2}(\mathcal{l} - \Delta_1 - .5\Delta_2))}
$$

Where $\mathcal{l}$ is the life expectancy of a given country and $g(\mathcal{l}, \theta)$ produces the expected increase by the next time period which in our case is five year intervals. By adding an error variance we get the following formula

$$
\mathcal{l}_{t} - \mathcal{l}_{t+1} \sim \mathcal{N}(g(\mathcal{l}_t, \theta), \sigma)
$$

b) Plot the observed gains against their fitted values, and comment on the fit.

![](/home/nmarquez//Documents/Classes/statdemog/week7/fitteddle0.jpg '')

The terms in the model have a degree of colinearity and can reach multiple solutions to provide similar goodness of fits which lead to the optimizer finalizing at relative convergence. Despite this the model accurately describes the data and seems to catch the latest trend well.
The final estimates for the parameters $\Delta_{1,2,3,4}$, $k$, $z$ and $\sigma$ are

```
1.07999999202428e-07, 62.7697153840251, 20.8251164413887, 9.83464768360668, 4.37345511908253, 0.870946870772529, 0.240632705628447
```

c) Find the predictive distribution of Honduras e0 for 2015-2020 conditional on this model, analytically or by simulation. Plot the distribution and give its median and a 95% prediction interval.

Since the model follows a normal distribution we can get the first years out of sample mean and median by plugging in the estimated parameters for the mean function $g(\mathcal{l}, \theta)$ and by adding and subtracting $1.96 \sigma$ from the mean we can get the 95% confidence intervals which for $\mathcal{l}_{2015}$ is

```
   2.5%     median    97.5%
75.78505 76.25669 76.72833
```

Below is the density plot for the distribution.

![](/home/nmarquez//Documents/Classes/statdemog/week7/e0honddens.jpg '')

2. A fully converged simulation for Bayesian modeling and projection of female life expectancy at birth is available at
http://www.stat.washington.edu/raftery/Stat593/Homework/e0simPAA16.tgz
This consists of three MCMC chains, each of length 160,000, thinned by 50, and 1,000 project trajectories for both female and male life expectancy. After unpacking you will find a README file that contains the code used to generate the simulation. 

a) Use the get.e0.mcmc and get.e0.prediction functions to obtain the MCMC and
prediction objects, respectively.

The data from above was downloaded and manipulated using the `bayesLife` package in `R`.

b) Pick any two countries, and assess the double logistic fit for each of the two countries

![](/home/nmarquez//Documents/Classes/statdemog/week7/bayesdlcompare.jpg '')

The model fits shown above are for the life expectancy double logistic curves of the countries Japan and Cambodia. They show the data in black, the median and quantiles in red and draws from the posterior in grey. Japan’s trajectory is more traditional of what we expect for the double logistic model and the curve fits all data points well. Although Cambodia has outliers in the data and has not yet hit the second stage stabilization point it is still able to provide a double logistic shape because of the strength it pools from similar countries. It should be noted that the plots above are without sampling the error variance.

c) Compare the fitted model for the two countries. In which one has life expectancy been rising faster.

Cambodia is closer to its point of highest five year gains than Japan and has seen greater more recent improvements in life expectancy. While this may be the case the plots below show that the model predicts the gains in Cambodia life expectancy will stabilize at a value slightly lower than Japan’s stabilization point, seen in the $z$ plot. The summary statistics for Japans $z$ value are as follows

```
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
0.001528 0.437848 0.530408 0.504240 0.596603 0.652953

```
While for Cambodia they are

```
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
0.008826 0.426964 0.520116 0.493281 0.588284 0.652947

```

![](/home/nmarquez//Documents/Classes/statdemog/week7/bayesparamcomparej.jpg '')
![](/home/nmarquez//Documents/Classes/statdemog/week7/bayesparamcomparec.jpg '')

d) Compare their projections, taking account of uncertainty.

![](/home/nmarquez//Documents/Classes/statdemog/week7/bayese0compare.jpg '')

Looking at the median value for both countries by 2100 Japan is likely to still have a greater life expectancy than Cambodia. It should be noted however that the uncertainty is much wider in the Cambodia’s forecast due to it being much lower in the current life expectancy trajectory of the double logistic curve while japan is at its stabilization point.

e) Select one of the two countries you picked. Assuming conditional independence of
life expectancy gain between countries given the model parameters (a reasonable assumption), compute the probability that the female e0 of your country will be larger/smaller than the female of all your neighboring countries in all future time periods (until 2100).

For this exercise I chose Japan as my reference country and its neighbors that I choose were the following 

"China", "Russian Federation", "Republic of Korea", "Dem. People's Rep. of Korea".

Because of the inclusion of North Korea which has an extremely low e0 trajectory Japan was never shown to be the worst country of its neighbors. However the high trajectory of South Korea made it so that by 2100 the probability that Japan was the best country in terms of high life expectancy value was down to below 30%. The probabilities by year are as follows

|  best| worst| year|
|-----:|-----:|----:|
| 0.957|     0| 2015|
| 0.768|     0| 2020|
| 0.612|     0| 2025|
| 0.530|     0| 2030|
| 0.463|     0| 2035|
| 0.425|     0| 2040|
| 0.386|     0| 2045|
| 0.356|     0| 2050|
| 0.334|     0| 2055|
| 0.333|     0| 2060|
| 0.315|     0| 2065|
| 0.294|     0| 2070|
| 0.289|     0| 2075|
| 0.293|     0| 2080|
| 0.297|     0| 2085|
| 0.280|     0| 2090|
| 0.289|     0| 2095|

3. Using the converged TFR and e0 simulations provided, generate probabilistic population projections for all countries.

a) For Canada, generate and plot probabilistic projections of the following quantities to 2100

i. total population
ii. total male population
iii. total population over 65
iv. the potential support ratio, defined as the number of people aged 20–64 divided by the number of people 65 and over

![](/home/nmarquez//Documents/Classes/statdemog/week7/CNpopproj.jpg '')
![](/home/nmarquez//Documents/Classes/statdemog/week7/CNpopprojmale.jpg '')
![](/home/nmarquez//Documents/Classes/statdemog/week7/CNpopproj65plus.jpg '')
![](/home/nmarquez//Documents/Classes/statdemog/week7/CNpopprojpsrm.jpg '')
![](/home/nmarquez//Documents/Classes/statdemog/week7/CNpopprojpsrf.jpg '')

The potential support ratio for Japan is broken down by males and females above. The potential support ratio is a reflection of the age structure of the country in question and higher numbers usually reflect a population that is high in TFR with higher fertility rates. As countries move through phase 3 TFR the potential support ratio should have strong declines. In the plots above we see this with Japan who is well into the phase 3 TFR process and which shows a steep declines in the potential support ratio bot in the past and in the future predictions.

b) Generate and plot a probabilistic projection of the total population of North America to 2100.

The aggregation of countries population projections can be done in the `bayesPop` library for `R`
By using the `pop.aggregate` function and giving the function a particular region code. In our case North America is `905`. The population forecast for both sexes all ages is shown below.

![](/home/nmarquez//Documents/Classes/statdemog/week7/NApopproj.jpg '')


## Code Appendix

```R
rm(list=ls())
set.seed(123)
par(mfrow=c(1,1))
pacman::p_load(pracma, data.table, wpp2015, bayesTFR, knitr, bayesPop, dplyr,
               ggplot2, forecast, bayesLife)

# Question 1

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
```