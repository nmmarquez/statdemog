rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, ggplot2, knitr, HPbayes, demogR)

data(mxF)
DF <- subset(as.data.table(mxF), country == "Mexico", 
             select=c("age", "2005-2010"))
setnames(DF, c("age", "asmr"))
DF$age <- as.numeric(as.character(DF$age))

#2a

gmp <- lm(log(asmr) ~ age, data=subset(DF, age >= 50))
summary(gmp)

gmp_predict <- function(p, age){
    alp <- exp(p[1])
    bet <- p[2]
    if (length(p) == 2){
        yhat <-alp * exp(bet * age)
    }
    else{
        gam <- p[3]
        yhat <- alp * exp(bet * age) + gam
    }
    return(yhat)
}

gmpfunc <- function(p, data=subset(DF, age >= 50)){
    yhat <- gmp_predict(p, data$age)
    return(mean((log(data$asmr) - log(yhat))**2)**.5)
}

gmpp <- optim(c(0,0), gmpfunc)
c(gmpp$par[1], gmpp$par[2])
gmp$coefficients
gmp_predict(gmpp$par, subset(DF, age >= 50)$age)
exp(predict(gmp))
subset(DF, age >= 50)$asmr

gmpm <- optim(c(gmpp$par, 0), gmpfunc)
gmpm
gmp_predict(gmpm$par, subset(DF, age >= 50)$age)

# 2b
DF[,qx:=1 - (1-asmr)**5]
DF[age == 0, qx:=1 - (1-asmr)**1]
DF[age == 1, qx:=1 - (1-asmr)**4]
DF$type <- "data"
DF2 <- copy(DF)
DF2$type <- "Gompertz"
DF2[,asmr:= gmp_predict(gmpp$par, age)]
DF2[,qx:=1 - (1-asmr)**5]
DF2[age == 0, qx:=1 - (1-asmr)**1]
DF2[age == 1, qx:=1 - (1-asmr)**4]
DF3 <- copy(DF)
DF3$type <- "Gompertz-Makeham"
DF3[,asmr:= gmp_predict(gmpm$par, age)]
DF3[,qx:=1 - (1-asmr)**5]
DF3[age == 0, qx:=1 - (1-asmr)**1]
DF3[age == 1, qx:=1 - (1-asmr)**4]

DFall <- rbindlist(list(DF, DF2, DF3))

jpeg('~/Documents/Classes/statdemog/week3/comparegomp.jpg')
ggplot(data=DFall, aes(age, log(asmr), color=type)) + geom_point() + 
    labs(title="Gompertz Compare")
dev.off()

#2c

hppred <- function(p, age){
    p[1]**((age + p[2])**-p[3]) +
        p[4] * exp(-p[5]*(log(age) - log(p[6]))**2) +
        p[7] * p[8]**age
}


hpfunc <- function(p, data=DF){
    yhat <- hppred(p, data$age)
    return(mean((logit(data$qx) - yhat)**2)**.5)
}

plot(logit(DF$qx))

hpls <- optim(c(0.033, 0.932, 0.204, 0.105, 3.821, 41.458, 0.001, 1.08), hpfunc)
hpls

DF4 <- copy(DF)
DF4$type <- "HP"
DF4[,qx:= inv.logit(hppred(hpls$par, age))]
DFHP <- rbindlist(list(DF, DF4))

ggplot(DFHP, aes(age, qx, color=type)) + geom_point()

jpeg('~/Documents/Classes/statdemog/week3/hpplot.jpg')
ggplot(DFHP, aes(age, logit(qx), color=type)) + geom_point() + 
    labs(title="HP Fit")
dev.off()

# 2d
N <- nrow(DF) - 2
rmse <- apply(cdmltw(sex = "F")$nmx, 1, function(x) 
    mean((x[1:N] - DF$asmr[1:N])**2)**.5)

mlt <- which(rmse == min(rmse))

mltqx <- logit(cdmltw(sex = "F")$nqx[mlt,1:N])
ltqx <- logit(DF$qx[1:N])

mltlm <- lm(ltqx ~ mltqx)
summary(mltlm)

DF5 <- copy(DF)
DF5$type <- "Model LT"
DF5[1:N, qx:= inv.logit(predict(mltlm))]
DF5[(N+1):nrow(DF5), qx:=NA]
DFHP <- rbindlist(list(DF, DF4, DF5))

ggplot(DFHP, aes(age, qx, color=type)) + geom_point()

jpeg('~/Documents/Classes/statdemog/week3/mltplot.jpg')
ggplot(subset(DFHP, type != "HP"), aes(age, logit(qx), color=type)) + 
    geom_point() + labs(title="Model Table Fit")
dev.off()

DFtotal <- rbindlist(list(DF, DF3, DF4, DF5))
ggplot(DFtotal, aes(age, qx, color=type)) + geom_point()

jpeg('~/Documents/Classes/statdemog/week3/allmodelcompare.jpg')
ggplot(DFtotal, aes(age, logit(qx), color=type)) + geom_point() + 
    labs(title="Model Compare")
dev.off()