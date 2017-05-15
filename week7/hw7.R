rm(list=ls())
set.seed(123)
pacman::p_load(pracma, data.table, wpp2015, bayesTFR, knitr, bayesPop, dplyr,
               ggplot2, forecast)

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
    preds <- e0_proj(params, e0)[1:(N-1)]
    return(sum(-1 * dnorm(diff(e0), preds, sigma, TRUE)))
}

plot(hond_e0[1:(N-1)], diff(hond_e0))
vals <- log(c(.77, 44.97, 20.21, .42, 3.93, 1.0, .4))
preds <- e0_proj(vals, hond_e0)[1:(N-1)]
lines(hond_e0[1:(N-1)], preds, col="red")

Obj <- nlminb(vals, objfunc1)
Obj$message
preds <- e0_proj(Obj$par, hond_e0)

jpeg('~/Documents/Classes/statdemog/week7/fitteddl.jpg')
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