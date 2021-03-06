# Homework Week 6  
### Neal Marquez  

1. Redo Homework 3, Question 3(b) once again. This time, solve it both
analytically (expressing the predictive distribution as a normal distribution,
conditional on the estimated parameters), and by simulation. Show your working
clearly.  

The Lee Carter model was fit two times, once using a Least Squares approach and another time using an SVD approach. For each of the models values for $a_x$,$b_x$, $\mu_k$, and $k_t$ were obtained and used to determine the mean forecast as well as $\sigma_k$, the error for the autoregressive $k$ terms and $\sigma_{\epsilon}$ the overall error in the model. Predictions were made for Mexican females aged 75-80 for the 2015-2020 time period. Uncertainty was generated by simulating data using 10000 simulations for the two error terms. The model was then forecasted to 2015 as follows

$log(m_{75,2015}) = a_{75} + b_{75} * k_{2015} + \epsilon_{75,2015}$  
$k_{2015} = k_{2010} + \mu_k + \epsilon_{2015}$  

The uncertainty for 2015 estimates of mortality can also be calculated analytically. The first source of error comes from the $k$ term which has the following distribution.

$k_{2015} \sim \mathcal{N}(k_{2010} + \mu_k, \sigma_k)$

the $k_{2015}$ term is then multiplied by the scalar value for $b_{75}$ which alters the variance for the distribution such that the resulting variance is $b_{75}^2 \sigma_k^2$. This resulting distribution is normal and can be combined with the variance from the overall error term which is also normal and the variances can be added such that the final variance of the estimates are now $b_{75}^2 \sigma_k^2 + \sigma^2_{\epsilon}$. This total variance can be converted to a standard deviation and 1.96 times the standard deviation can then be added and subtracted from the mean values in order to get the 2.5% and the 97.5% quantiles.

The values generated by the two methods are shown below, first for the least squares approach and then for the svd approach, where ymin and ymax are the intervals that were simulated and ymin2 and ymax2 are the intervals that were created analytically. Mxt is the median mortality rate for women in mexico age 75 in 2015.

### least squares method
|     mxt|      ymin|      ymax|     ymin2|     ymax2|
|---------:|---------:|---------:|---------:|---------:|
| 0.0366498| 0.0321386| 0.0416743| 0.0322563| 0.0416417|

### svd method
|     mxt|      ymin|      ymax|     ymin2|     ymax2|
|---------:|---------:|---------:|---------:|---------:|
| 0.0371044| 0.0328375| 0.0419152| 0.0327884| 0.0419885|

Below are the plots for the least squares model and the svd model respectively.

![](/home/nmarquez/Documents/Classes/statdemog/week6/lsleecarter.jpg "")
![](/home/nmarquez/Documents/Classes/statdemog/week6/svdleecarter.jpg "")

2. Obtain the values of TFR for Honduras for 1950-2015 from the 2015 World Population Prospects.

a. Fit a version of the double logistic decline model by nonlinear least squares to
these data, assuming that the error variance remains constant over time (this
will just give one set of double logistic parameter values and an estimated error
variance). Note: Honduras is still in Phase II up to 2015.

The double logistic model uses two logistic equations in order to capture the change between time periods of measured total fertility rate as a country moves from phase two to three. The model is fit using the `nlminb` optimizer in `R` in order to fit the 5 parameters of the model as well as the error variance term. Because Honduras is still in Phase 2 the model was unable to converge on a final answer as the parameters $\Delta_{1,2,3,4}$ dictate the time at which phase 3 starts and without data from this phase the model fails to choose a single set of parameters to converge on. Nevertheless the fit is appropriate for the data. Note that the functional form of the model is as follows

$$
g(f, \theta) \sim \frac{-d}{1 + exp(\frac{-2ln(9)}{\Delta_1}(f - \sum_{i=2}^{4} \Delta_i + .5\Delta_1))} + \frac{d}{1 + exp(\frac{-2ln(9)}{\Delta_3}(f - \Delta_4 + .5\Delta_3))}
$$

Where f is the total fertility rate of a given country and $g(f, \theta)$ produces the expected decrease by the next time period which in our case is five year intervals. By adding an error variance we get the following formula

$$
f_{t} - f_{t+1} \sim \mathcal{N}(g(f_t, \theta), \sigma)
$$

b. Plot the observed declines against their fitted values, and comment on the fit.

![](/home/nmarquez/Documents/Classes/statdemog/week6/tfrdiff.jpg "")  

The fit for the first logistic curve fits the increase in phase 2 declines well, however, the model uses the second logistic regression to capture one data point at the end of the TFR time series. The final fitted values were very different depending on the selected start values. The final estimates for the parameters $\Delta_{12,3,4}$, $d$ and $\sigma$ are

```
0.6656653 3.2326534 0.1940434 2.8104614 0.6024593 0.1005825
```
resepctively.

c. Find the predictive distribution of Honduras TFR for 2015-2020 conditional on
this model, analytically or by simulation. Plot the distribution and give its median
and a 95% prediction interval.

Since the model follows a normal distribution we can get the first years out of sample mean and median by plugging in the estimated parameters for the mean function $g(f, \theta)$ and by adding and subtracting $1.96 \sigma$ from the mean we can get the 95% confidence intervals which for $f_{2015}$ is

```
   2.5%     median    97.5%
2.272828 2.469970 2.667112
```

Below is the density plot for the distribution.

![](/home/nmarquez/Documents/Classes/statdemog/week6/tfrhonddens.jpg "")

3. Obtain the values of TFR for the Netherlands for 1950-2015 from the 2015 World
Population Prospects.

a. Identify the period in which the Netherlands entered Phase III of the fertility
Model.

Phase III of the fertility model begins when a country has two five year increases of fertility which are below a total tfr value of 2.  For the Netherlands this year is 1985 which can be seen in the plot below

![](/home/nmarquez/Documents/Classes/statdemog/week6/ntlts.jpg "")

b.  Fit a (non-Bayesian) AR(1) model to the Phase III data, estimating the long-term
mean, autoregressive parameter, and error variance.

AR1 models follow the following predictive distribution

$$
X_t \sim \mathcal{N}(c + \rho X_{t-1}, \sigma)
$$

Where the long term mean of the process can be calculated as follows

$$
\mu = \frac{c}{1-\rho}
$$

$\rho$ and $\sigma^2$ represent the autoregressive parameter and the error variance
respectively. Autoregressive integrated moving average (ARIMA) models, which include the AR1 model, can be estimated in `R` using the forecast package. An AR1 model corresponds to an ARIMA model of class (1,0,0). The values for the long term mean, autoregressive parameter, and error variance are

```
1.6563  0.7682  0.0048
```

respectively. A plot of the forecasts to 2060 are shown below with the 80% and 95% prediction intervals.

![](/home/nmarquez/Documents/Classes/statdemog/week6/ntltsforecast.jpg "")

c) Find the predictive distribution of Honduras TFR for 2015-2020 conditional on
this model, analytically or by simulation. Plot the distribution and give its median
and a 95% prediction interval.

Since this model follows a normal distribution we can calculate the uncertainty for the 2015 time period analytically by taking the mean/median which is $c + \rho * tfr_{2010}$ and adding and subtracting $1.96 \sigma$ in order to get the upper and lower bounds which are as follows

```
2.5%     median 97.5%
1.594047 1.729517 1.864988
```

Below is the density plot for the distribution.

![](/home/nmarquez/Documents/Classes/statdemog/week6/ntldens2015.jpg "")



4. Download a fully converged simulation containing three MCMC chains from phase II, each of
length 62,000, and three MCMC chains of Phase III, each of length 70,000, both thinned
by 30, and 1,000 projection trajectories for all countries.

a. Use the get.tfr.mcmc, get.tfr3.mcmc, and get.tfr.prediction functions to obtain
the MCMC (II and III) and the prediction objects, respectively.

MCMC traces can be obtained and manipulated in `R` by using the `get.tfr.mcmc` and `get.tfr3.mcmc` functions which are a part of `bayesTFR`. The model's predictions can be obtained by using the `get.tfr.predictions` function from the same package. Predictions use the tfr2 and tfr3 models in order to generate probabilistic forecasts for countries where phase III has not completed.

b. Using the converged simulation, pick two countries that have current TFR larger
than 2.1 and assess the double logistic fit for each of the two countries.

Below are the fitted double logistic curve plots for TFR decrease as a function of TFR level for the countries Guatemala and Cambodia which in 2010 both had  TFR values of greater than 2.1 by UN population estimates. The model uses a hierarchical structure in order to pool strength from adjacent countries and their patterns of TFR phase completion for the $\Delta$ and $d$ parameters. You can see that even though Cambodia observes a strong outlier in the data the fit is still what we would typically view as a normal fertility transition pattern. Note that the median is shown in the solid red line and the 95% confidence intervals with a dotted red line. The yellow dots represent data and black lines are single draws of the posterior.

![](/home/nmarquez/Documents/Classes/statdemog/week6/bayestfrdlc.jpg "")

c. Compare the fitted model for the two countries: in which one has fertility been
declining faster?

Both Guatemala and Cambodia have experienced rapid declines in the near past. It appears that both are just past the peaks of their declines and the model fits show that cambodia has a higher upper bound for expected decline, median value of `0.0967`, than Guatemala, median value of `0.0876`, when accounting for 2000 simulations of burn in for the chains. The density plots for each of the parameters are shown below. 

![](/home/nmarquez/Documents/Classes/statdemog/week6/bayestfrdens1.jpg "") 
![](/home/nmarquez/Documents/Classes/statdemog/week6/bayestfrdens2.jpg "")  

d. Compare their projections, taking account of uncertainty.

![](/home/nmarquez/Documents/Classes/statdemog/week6/projfor.jpg "")

Because Cambodia has a lower value for TFR it appears to be further along the TFR transition and more often begins phase III before Guatemala. The uncertainty for the point estimates for the two countries look very similar and by 2100 it appears as if both countries have stabilized at their phase III constant value.  

### Code Appendix 

```
rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, bayesTFR, knitr, bayesPop, dplyr,
               ggplot2, forecast)

# question 1
data(mxF)
DF <- subset(as.data.table(mxF), country == "Mexico", 
             select=c("age", "2005-2010"))
setnames(DF, c("age", "asmr"))
DF$age <- as.numeric(as.character(DF$age))
mxt <- log(as.matrix(subset(mxF, country == "Mexico")[,4:16]))
ax <- apply(mxt, 1, mean)
kt <- apply(mxt - ax, 2, mean)
meanadjxt <- mxt - ax 

bx <- apply(meanadjxt, 1, function(x) lm(x ~ 0 + kt)$coefficients[[1]])
nu <- (kt[length(kt)] - kt[1]) / (length(kt) - 1)
sigma_k <- sd(kt[2:length(kt)] - kt[1:(length(kt) -1)])
eps <- sd(mxt - (bx %*% t(kt) + ax))

kt_forecast <- c(kt, kt[length(kt)] + nu)


DFforecast <- data.table(lnmxt=c(bx %*% t(kt_forecast) + ax), 
                         age=rep(DF$age, length(kt_forecast)),
                         year=rep(seq(1950, 2015, 5), each=nrow(DF)),
                         model="ls")


# simulate uncertainty 
sims <- 10000
simvals <- ax[17] + bx[17] * rnorm(sims, kt[length(kt)] + nu, sigma_k) + 
    rnorm(sims, sd=eps)
quantile(simvals, .975)
quantile(simvals, .025)

DFforecast[,ymin:= lnmxt]
DFforecast[,ymax:= lnmxt]
DFforecast[year == 2015 & age == 75, ymin:= quantile(simvals, .025)]
DFforecast[year == 2015 & age == 75, ymax:= quantile(simvals, .975)]

# analytically get uncertainty
toterr <- sqrt(sigma_k**2 * bx[17]**2 + eps**2)
DFforecast[,ymin2:= lnmxt]
DFforecast[,ymax2:= lnmxt]
DFforecast[year == 2015 & age == 75, ymin2:= lnmxt - 1.96 * toterr]
DFforecast[year == 2015 & age == 75, ymax2:= lnmxt + 1.96 * toterr]


jpeg('~/Documents/Classes/statdemog/week6/lsleecarter.jpg')
ggplot(DFforecast[age==75,], aes(x=year, y=lnmxt, color=age, group=age)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
    labs(title="LS Lee-Carter")
dev.off()

DFforecast[year == 2015 & age == 75,list(lnmxt, ymin, ymax, ymin2, ymax2)] %>%
    exp %>% kable(format="markdown")



lm.t.x <- t(mxt)
Y <- sweep(lm.t.x, 2, ax)             #mean centered data
Y.svd <- svd(Y)                       #returns U %*% diag(d) %*% t(V) = Y
bxsvd <- Y.svd$v[,1]
b1sign <- sign(bxsvd[1])
ktsvd <- Y.svd$d[1]*Y.svd$u[,1]
bxsvd <- bxsvd*b1sign 
ktsvd <- ktsvd*b1sign
scalefactor <- length(bxsvd) / sum(bxsvd)
ktsvd <- ktsvd / scalefactor
bxsvd <- bxsvd * scalefactor
nusvd <- (ktsvd[length(ktsvd)] - ktsvd[1]) / (length(ktsvd) - 1)
sigma_ksvd <- sd(ktsvd[2:length(ktsvd)] - ktsvd[1:(length(ktsvd) -1)])
ktsvd_forecast <- c(ktsvd, ktsvd[length(ktsvd)] + nusvd)
epssvd <- sd(mxt - (bxsvd %*% t(ktsvd) + ax))


DFforecastsvd <- data.table(lnmxt=c(bxsvd %*% t(ktsvd_forecast) + ax), 
                            age=rep(DF$age, length(ktsvd_forecast)),
                            year=rep(seq(1950, 2015, 5), each=nrow(DF)),
                            model="svd")

sims <- 10000
simvalssvd <- ax[17] + bxsvd[17] * rnorm(sims, ktsvd[length(ktsvd)] + nusvd, sigma_ksvd) + 
    rnorm(sims, sd=epssvd)

DFforecastsvd[,ymin:= lnmxt]
DFforecastsvd[,ymax:= lnmxt]
DFforecastsvd[year == 2015 & age == 75, ymin:= quantile(simvalssvd, .025)]
DFforecastsvd[year == 2015 & age == 75, ymax:= quantile(simvalssvd, .975)]
DFforecastsvd[year == 2015 & age == 75,]

# analytically get uncertainty
toterrsvd <- sqrt(sigma_ksvd**2 * bxsvd[17]**2 + epssvd**2)
DFforecastsvd[,ymin2:= lnmxt]
DFforecastsvd[,ymax2:= lnmxt]
DFforecastsvd[year == 2015 & age == 75, ymin2:= lnmxt - 1.96 * toterrsvd]
DFforecastsvd[year == 2015 & age == 75, ymax2:= lnmxt + 1.96 * toterrsvd]
DFforecastsvd[year == 2015 & age == 75,]

jpeg('~/Documents/Classes/statdemog/week6/svdleecarter.jpg')
ggplot(DFforecastsvd[age==75,], aes(x=year, y=lnmxt, color=age, group=age)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
    labs(title="SVD Lee Carter")
dev.off()

DFforecastsvd[year == 2015 & age == 75,list(lnmxt, ymin, ymax, ymin2, ymax2)] %>%
    exp %>% kable(format="markdown")

# TFR Hondurus 
rm(list=ls())

data(tfr)
hond_tfr <- c(as.matrix(subset(tfr, country == "Honduras")[,3:(ncol(tfr) - 1)]))
plot(hond_tfr)
hond_diff <- -1 * diff(hond_tfr)
plot(hond_tfr[1:length(hond_diff)], hond_diff)


tfr_proj1 <- function(params, tfr=hond_tfr){
    a <- exp(params[1:4])
    d <- exp(params[5])
    dl1 <- d / (1 + exp(-(tfr - a[2])/a[1]))
    dl2 <- d / (1 + exp(-(tfr - a[4])/a[3]))
    return(dl1 - dl2)
}

tfr_proj2 <- function(params, tfr=hond_tfr){
    a <- exp(params[1:4])
    d <- exp(params[5])
    dl1 <- -d/(1 + exp((-2 * log(9))/a[1] * (tfr - sum(a[2:4] + a[1]/2))))
    dl2 <- d/(1 + exp((-2 * log(9))/a[3] * (tfr - a[4] - .5* a[3])))
    return(dl1 + dl2)
}

objfunc1 <- function(params, tfr=hond_tfr){
    N <- length(tfr)
    sigma <- exp(params[6])
    preds <- tfr_proj2(params, tfr)
    return(sum(-1 * dnorm(tfr[2:N], (tfr - preds)[1:(N-1)], sigma, TRUE)))
}

Obj <- nlminb(log(c(1.1, 3.2, 1.0, 1.0, 0.6, 1)), objfunc1)
Obj$message
preds <- tfr_proj2(Obj$par)
par(mfrow=c(1,1))
plot(hond_tfr, preds)
years <- seq(1950, 2010, 5)

jpeg('~/Documents/Classes/statdemog/week6/tfrest.jpg')
plot(years, hond_tfr, xlab="Year", ylab="TFR", main="TFR estimates")
lines(years + 5, (hond_tfr - preds), col="red")
dev.off()

line_fit <- seq(2, 7.5, .1)
m <- length(line_fit)

jpeg('~/Documents/Classes/statdemog/week6/tfrdiff.jpg')
plot(hond_tfr[1:(length(hond_tfr)-1)], 
     -1 * diff(hond_tfr), xlim=rev(range(hond_tfr)), 
     xlab="TFR", ylab=expression(TFR[t] - TFR[t+1]),
      main=paste0("LS Fit of change in TFR"))
lines(line_fit, tfr_proj2(Obj$par, line_fit), col="red")
dev.off()

tail(hond_tfr)
mu2015 <- (hond_tfr - preds)[length(preds)]
eps <- exp(Obj$par[6])
c(`2.5%`=mu2015 - 1.96*eps, mean=mu2015, `97.5%`=mu2015 + 1.96*eps)

x <- seq(mu2015 - 4*eps, 
         mu2015 + 4*eps, length=100)
hx <- dnorm(x, mean=mu2015, sd=eps)
jpeg('~/Documents/Classes/statdemog/week6/tfrhonddens.jpg')
plot(x, hx, "l", xlab="2015 Honduras TFR Forecast", ylab="Density")
dev.off()

# TFR Netherlands
rm(list=ls())
data(tfr)
ntfr <- c(as.matrix(subset(tfr, country == "Netherlands")[,3:(ncol(tfr) - 1)]))
years <- seq(1950, 2010, 5)

for(i in 1:(length(ntfr)-2)){
    if(ntfr[i] < 2){
        if(ntfr[i+2] > ntfr[i+1] & ntfr[i+1] > ntfr[i]){
            phase3 <- i+1
            break
        }
    }
}

cols <- c("red", "green")[(1:length(ntfr) >= phase3) + 1]


jpeg('~/Documents/Classes/statdemog/week6/ntlts.jpg')
plot(years, ntfr, pch=20, col=cols, xlab="Year", ylab="TFR",
     main="Netherlands Time Series")
legend("topright", legend=c("Not Phase 3", "Phase 3"), col=c("red", "green"),
       pch=20)
dev.off()

ntfr3 <- ts(ntfr[phase3:length(ntfr)], start=years[phase3], 
            end=years[length(years)], frequency = 1/5)
plot(ntfr3)

m1 <- Arima(ntfr3, c(1,0,0))
summary(m1)
preds <- forecast(m1,10)

jpeg('~/Documents/Classes/statdemog/week6/ntltsforecast.jpg')
plot(preds, main="Netherlands Forecast", ylab="TFR", xlab="Year")
dev.off()

mu2015 <- preds$mean[1]
eps <- m1$sigma2**.5
c(`2.5%`=mu2015 - 1.96*eps, mean=mu2015, `97.5%`=mu2015 + 1.96*eps)
x <- seq(mu2015 - 4*eps, 
         mu2015 + 4*eps, length=100)
hx <- dnorm(x, mean=mu2015, sd=eps)

jpeg('~/Documents/Classes/statdemog/week6/ntldens2015.jpg')
plot(x, hx, "l", xlab="2015 Netherlands TFR Forecast", ylab="Density")
dev.off()

# bayes tfr compare

rm(list = ls())

tfdir <- '~/Downloads/TFR/sim03092016'
m2 <- get.tfr.mcmc(tfdir)
m3 <- get.tfr3.mcmc(tfdir)
pred <- get.tfr.prediction(tfdir)



jpeg('~/Documents/Classes/statdemog/week6/bayestfrdlc.jpg')
par(mfrow=c(1,2))
DLcurve.plot(m2, "Guatemala", col=c("yellow", "red", "black"), burnin=2000,
             pi=95, show.legend=FALSE)
DLcurve.plot(m2, "Cambodia", col=c("yellow", "red", "black"), burnin=2000,
             pi=95, show.legend=FALSE)
dev.off()

jpeg('~/Documents/Classes/statdemog/week6/bayestfrdens1.jpg')
tfr.pardensity.cs.plot("Guatemala",m2, par.names="d")
dev.off()
jpeg('~/Documents/Classes/statdemog/week6/bayestfrdens2.jpg')
tfr.pardensity.cs.plot("Cambodia",m2, par.names="d")
dev.off()

get.tfr.parameter.traces.cs(m2$mcmc.list, 
                            get.country.object("Guatemala", meta=m2$meta), 
                            par.names = "d", burnin = 100, thin=30) %>% median

get.tfr.parameter.traces.cs(m2$mcmc.list, 
                            get.country.object("Cambodia", meta=m2$meta), 
                            par.names = "d", burnin = 100, thin=30) %>% median


jpeg('~/Documents/Classes/statdemog/week6/projfor.jpg')
par(mfrow=c(1,2))
tfr.trajectories.plot(pred, country="Guatemala")
tfr.trajectories.plot(pred, country="Cambodia")
dev.off()
```
