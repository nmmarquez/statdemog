rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, ggplot2, knitr, bayesPop)

# question 1  

asmig <- subset(age.specific.migration(years=2010)$female, name == "Mexico")

leslie_project <- function(les, v, proj=1, migrate=0){
    for(i in 1:proj){
        migpop <- v * migrate
        v <- c(les %*% v) + migpop
    }
    newpop <- v
    return(newpop)
}

data(mxF, pop, percentASFR, tfr)
DF <- subset(as.data.table(mxF), country == "Mexico", 
             select=c("age", "2005-2010"))
setnames(DF, c("age", "asmr"))

# we need to combine the first two age groups
DF$age <- as.numeric(as.character(DF$age))
DF[age == 0, asmr:= (DF$asmr[1] + 4 * DF$asmr[2]) / 5]
DF <- subset(DF, age != 1)
DF$pop <- subset(popF, country == "Mexico")$`2005`

ASFR <- subset(percentASFR, country=="Mexico")$`2005` * 
    subset(tfr, country == "Mexico")$`2005-2010` / 500

DF[,asfr:=0]
DF[age >=15 & age < 50,asfr:=ASFR]
DF[,asmigr:=asmig[["2005-2010"]] / pop ]

rm(list=c("mxF", "pop", "popF", "popM", "tfr", "percentASFR", "ASFR"))

jpeg('~/Documents/Classes/statdemog/week2/migstructure.jpg')
ggplot(data=DF, aes(x=age, y=asmigr)) + geom_area() + 
    labs(y="Migration Rate per Person", title="Migration Structure")
dev.off()

#b
fert <- DF$asfr
age.int <- 5
srb <- 1.05
surv <- c(1, (1 - DF$asmr)**5)
k <- 1/(1 + srb) * surv[1] * 0.5
dbl.fert <- age.int * fert + c(age.int * fert[-1], 0) * surv[-1]
k * dbl.fert

#c
n.age.grps <- length(fert)
n.surv <- length(surv)
lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)
lesM[1, ] <- k * dbl.fert
lesM[2:n.age.grps, 1:(n.age.grps - 1)] <- diag(surv[-c(1, 
                                                       n.surv)])
lesM[n.age.grps, n.age.grps] <- surv[n.surv]

for(i in 1:nrow(lesM)){
    for(j in 1:ncol(lesM)){
        if (lesM[i,j] != 0){
            cat(paste0("(", i, ",", j, "): ", sprintf("%05f", lesM[i,j])))
            cat("\n")
        }
    }
}

#d
leslie_project(lesM, DF$pop, 1)

#e
leslie_project(lesM, DF$pop, 3)

DFcompare <- data.table(age=DF$age, pop=leslie_project(lesM, DF$pop, 1),
                        year=2010, migration=F)
DFcompare2 <- copy(DFcompare)
DFcompare2$pop <- leslie_project(lesM, DF$pop, 3)
DFcompare2[,year:=2020]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/compareyear.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=as.factor(year))) + geom_point(alpha=.5) + 
    labs(title="Population Projections No Migration")
dev.off()

#f
migrate <- DF$asmigr
migrate

#g
leslie_project(lesM, DF$pop, 1, migrate)
leslie_project(lesM, DF$pop, 3, migrate)

DFcompare <- data.table(age=DF$age, pop=leslie_project(lesM, DF$pop, 1, migrate),
                        year=2010, migration=F)
DFcompare2 <- copy(DFcompare)
DFcompare2$pop <- leslie_project(lesM, DF$pop, 3, migrate)
DFcompare2[,year:=2020]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/compareyear2.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=as.factor(year))) + geom_point(alpha=.5) + 
    labs(title="Population Projections With Migration")
dev.off()

#h
DFcompare <- data.table(age=DF$age, pop=leslie_project(lesM, DF$pop, 1),
                        year=2010, migration=F)
DFcompare2 <- copy(DFcompare)
DFcompare2$pop <- leslie_project(lesM, DF$pop, 1, migrate)
DFcompare2[,migration:=T]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/comparemig2010.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=migration)) + geom_point(alpha=.5) + 
    labs(title="Population Projections 2010")
dev.off()

DFcompare <- data.table(age=DF$age, pop=leslie_project(lesM, DF$pop, 3),
                        year=2020, migration=F)
DFcompare2 <- copy(DFcompare)
DFcompare2$pop <- leslie_project(lesM, DF$pop, 3, migrate)
DFcompare2[,migration:=T]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/comparemig2020.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=migration)) + geom_point(alpha=.5) + 
    labs(title="Population Projections 2020")
dev.off()



# question 2
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

# Question 3

mxt <- log(as.matrix(subset(mxF, country == "Mexico")[,4:16]))
ax <- apply(mxt, 1, mean)
kt <- apply(mxt - ax, 2, mean)
meanadjxt <- mxt - ax 

bx <- apply(meanadjxt, 1, function(x) lm(kt ~ 0 + x)$coefficients[[1]])
nu <- (kt[length(kt)] - kt[1]) / (length(kt) - 1)
sigma_k <- sd(kt[2:length(kt)] - kt[1:(length(kt) -1)])
eps <- sd(mxt - (bx %*% t(kt) + ax))

kt_forecast <- c(kt, kt[length(kt)] + nu)
ktsvd_forecast <- c(ktsvd, ktsvd[length(ktsvd)] + nusvd)
DFforecast <- data.table(lnmxt=c(bx %*% t(kt_forecast) + ax), 
                         age=rep(DF$age, length(kt_forecast)),
                         year=rep(seq(1950, 2015, 5), each=nrow(DF)),
                         model="ls")

DFforecast[,ymin:=lnmxt - 1.96*sigma_k]
DFforecast[year != 2015, ymin:=lnmxt]
DFforecast[,ymax:=lnmxt + 1.96*sigma_k]
DFforecast[year != 2015, ymax:=lnmxt]

jpeg('~/Documents/Classes/statdemog/week3/lsleecarter.jpg')
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

jpeg('~/Documents/Classes/statdemog/week3/svdleecarter.jpg')
ggplot(DFforecastsvd, aes(x=year, y=lnmxt, color=age, group=age)) + geom_line() +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) + 
    labs(title="SVD Lee Carter")
dev.off()

DFobs <- data.table(lnmxt=c(mxt), age=rep(DF$age, length(ktsvd)),
                    year=rep(seq(1950, 2010, 5), each=nrow(DF)), model="data",
                    ymin=NA, ymax=NA)
ggplot(DFobs, aes(x=year, y=lnmxt, color=age, group=age)) + geom_line()

DFfff <- rbindlist(list(DFobs, DFforecast))

jpeg('~/Documents/Classes/statdemog/week3/allleecarterls.jpg')
ggplot(DFfff, aes(x=year, y=lnmxt, color=age, 
                  group=interaction(age, model), linetype=model)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) +
    labs(title="Lee Carter LS against Data")
dev.off()

DFfff2 <- rbindlist(list(DFobs, DFforecastsvd))

jpeg('~/Documents/Classes/statdemog/week3/allleecartersvd.jpg')
ggplot(DFfff2, aes(x=year, y=lnmxt, color=age, 
                  group=interaction(age, model), linetype=model)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) +
    labs(title="Lee Carter SVD against Data")
dev.off()


jpeg('~/Documents/Classes/statdemog/week3/allleecarter75.jpg')
ggplot(subset(DFfff, age == 75), aes(x=year, y=lnmxt, color=age, 
                  group=interaction(age, model), linetype=model)) + 
    geom_line() + geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=age), alpha=.5) +
    labs(title="Lee Carter Compare")
dev.off()
