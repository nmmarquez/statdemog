rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, ggplot2, knitr, bayesPop, mlmRev)

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


jpeg('~/Documents/Classes/statdemog/week5/lsleecarter.jpg')
ggplot(DFforecast[age==75,], aes(x=year, y=lnmxt, color=age, group=age)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
    labs(title="LS Lee-Carter")
dev.off()



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

# analytically get uncertainty
toterrsvd <- sqrt(sigma_ksvd**2 * bxsvd[17]**2 + epssvd**2)
DFforecastsvd[,ymin2:= lnmxt]
DFforecastsvd[,ymax2:= lnmxt]
DFforecastsvd[year == 2015 & age == 75, ymin2:= lnmxt - 1.96 * toterrsvd]
DFforecastsvd[year == 2015 & age == 75, ymax2:= lnmxt + 1.96 * toterrsvd]
DFforecastsvd[year == 2015 & age == 75,]

jpeg('~/Documents/Classes/statdemog/week5/svdleecarter.jpg')
ggplot(DFforecastsvd[age==75,], aes(x=year, y=lnmxt, color=age, group=age)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
    labs(title="SVD Lee Carter")
dev.off()

DFobs <- data.table(lnmxt=c(mxt), age=rep(DF$age, length(ktsvd)),
                    year=rep(seq(1950, 2010, 5), each=nrow(DF)), model="data",
                    ymin=NA, ymax=NA)
ggplot(DFobs, aes(x=year, y=lnmxt, color=age, group=age)) + geom_line()

DFfff <- rbindlist(list(DFobs, DFforecast))

jpeg('~/Documents/Classes/statdemog/week5/allleecarterls.jpg')
ggplot(DFfff, aes(x=year, y=lnmxt, color=age, 
                  group=interaction(age, model), linetype=model)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) +
    labs(title="Lee Carter LS against Data", y="Log Female Mortality in Mexico")
dev.off()