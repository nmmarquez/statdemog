rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, ggplot2, knitr, bayesPop, dplyr)

# question 2
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
