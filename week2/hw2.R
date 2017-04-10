rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, ggplot2, knitr, wppExplorer)

DF <- data.table(age=1:3, pop=c(18, 17, 14) * 1000, 
                 fert=c(0, .9, .2), surv=c(.65, .75, .15))

calc_fert <- function(DF, crude=T){
    if(!crude){
        return(sum(DF$fert * DF$pop / DF$pop))
    }
    return(sum(DF$fert * DF$pop) / sum(DF$pop))
}

leslie_project <- function(les, v, proj=1, migrate=0){
    for(i in 1:proj){
        migpop <- v * migrate
        v <- c(les %*% v) + migpop
    }
    newpop <- v
    return(newpop)
}

# a
calc_fert(DF, T)

# b
calc_fert(DF, F)

#c
L <- matrix(0, 3, 3)
L[1,] <- DF$fert * DF$surv
L[2,1] <- DF$surv[1]
L[3,c(2,3)] <- DF$surv[2:3]

colnames(L) <- 1:3
kable(L, "markdown")

#d
leslie_project(L, DF$pop)

#e
leslie_project(L, DF$pop, 10)

#f 
DF2 <- copy(DF)
DF2$pop <- leslie_project(L, DF$pop, 10)

calc_fert(DF2)
calc_fert(DF2, FALSE)

#g
lamb <- eigen(L)$values[1]
log(lamb)

#h
u <- eigen(L)$vectors[,1]
u / sum(u)

#f
u_ <- eigen(t(L))$vectors[,1]
u_ / u_[1] # relative reproductive value


# 2 

# a
rm(list=setdiff(ls(), "leslie_project"))
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

rm(list=c("mxF", "pop", "popF", "popM", "tfr", "percentASFR", "ASFR"))

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
DFcompare2[,year:=2015]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/compareyear.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=as.factor(year))) + geom_point(alpha=.5) + 
    labs(title="Population Projections No Migration")
dev.off()

#f
migrate <- subset(wpp.indicator("migrate", sex="F"), charcode=="MX" & Year == 2005)$value / 1000
migrate

#g
leslie_project(lesM, DF$pop, 1, migrate)
leslie_project(lesM, DF$pop, 3, migrate)

DFcompare <- data.table(age=DF$age, pop=leslie_project(lesM, DF$pop, 1, migrate),
                        year=2010, migration=F)
DFcompare2 <- copy(DFcompare)
DFcompare2$pop <- leslie_project(lesM, DF$pop, 3, migrate)
DFcompare2[,year:=2015]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/compareyear2.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=as.factor(year))) + geom_point(alpha=.5) + 
    labs(title="Population Projections With Migration")
dev.off()

#h
DFcompare <- data.table(age=DF$age, pop=leslie_project(lesM, DF$pop, 3),
                        year=2020, migration=F)
DFcompare2 <- copy(DFcompare)
DFcompare2$pop <- leslie_project(lesM, DF$pop, 3, migrate)
DFcompare2[,migration:=T]
DFcompare<-rbind(DFcompare, DFcompare2)

jpeg('~/Documents/Classes/statdemog/week2/comparemig.jpg')
ggplot(DFcompare, aes(x=age, y=pop, color=migration)) + geom_point(alpha=.5) + 
    labs(title="Population Projections 2020")
dev.off()
