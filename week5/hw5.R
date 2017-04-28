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

DFforecast[,ymin:=lnmxt - 1.96*sigma_k]
DFforecast[year != 2015, ymin:=lnmxt]
DFforecast[,ymax:=lnmxt + 1.96*sigma_k]
DFforecast[year != 2015, ymax:=lnmxt]

jpeg('~/Documents/Classes/statdemog/week5/lsleecarter.jpg')
ggplot(DFforecast, aes(x=year, y=lnmxt, color=age, group=age)) + geom_line() +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
    labs(title="LS Lee-Carter")
dev.off()


lm.t.x <- t(mxt)
Y <- sweep(lm.t.x, 2, ax);             #mean centered data
Y.svd <- svd(Y);                       #returns U %*% diag(d) %*% t(V) = Y
bxsvd <- Y.svd$v[,1];
b1sign <- sign(bxsvd[1]);
ktsvd <- Y.svd$d[1]*Y.svd$u[,1];
bxsvd <- bxsvd*b1sign; 
ktsvd <- ktsvd*b1sign;
nusvd <- (ktsvd[length(ktsvd)] - ktsvd[1]) / (length(ktsvd) - 1)
sigma_ksvd <- sd(ktsvd[2:length(ktsvd)] - ktsvd[1:(length(ktsvd) -1)])
ktsvd_forecast <- c(ktsvd, ktsvd[length(ktsvd)] + nusvd)
epssvd <- sd(mxt - (bxsvd %*% t(ktsvd) + ax))


DFforecastsvd <- data.table(lnmxt=c(bxsvd %*% t(ktsvd_forecast) + ax), 
                            age=rep(DF$age, length(ktsvd_forecast)),
                            year=rep(seq(1950, 2015, 5), each=nrow(DF)),
                            model="svd")

DFforecastsvd[,ymin:=lnmxt - 1.96*sigma_ksvd]
DFforecastsvd[year != 2015, ymin:=lnmxt]
DFforecastsvd[,ymax:=lnmxt + 1.96*sigma_ksvd]
DFforecastsvd[year != 2015, ymax:=lnmxt]

jpeg('~/Documents/Classes/statdemog/week5/svdleecarter.jpg')
ggplot(DFforecastsvd, aes(x=year, y=lnmxt, color=age, group=age)) + geom_line() +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
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

DFfff2 <- rbindlist(list(DFobs, DFforecastsvd))

jpeg('~/Documents/Classes/statdemog/week5/allleecartersvd.jpg')
ggplot(DFfff2, aes(x=year, y=lnmxt, color=age, 
                   group=interaction(age, model), linetype=model)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) +
    labs(title="Lee Carter SVD against Data", y="Log Female Mortality in Mexico")
dev.off()


jpeg('~/Documents/Classes/statdemog/week5/allleecarter75.jpg')
ggplot(subset(DFfff, age == 75), aes(x=year, y=lnmxt, color=age, 
                                     group=interaction(age, model), linetype=model)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) +
    labs(title="Lee Carter Compare")
dev.off()

kable(subset(DFfff, year == 2015), format="markdown")
kable(subset(DFfff2, year == 2015), format="markdown")


# question 3
mlm <- lmer(math ~ 1 + (1|schoolid), data=subset(egsingle, year==.5))
mlm@beta

params <- c(mlm@beta, as.data.frame(summary(mlm)$varcor)$sdcor)
names(params) <- c("mu_alpha", "sigma_alpha", "sigma_epsilon")
params

sd(summary(subset(egsingle, year==.5)$math))