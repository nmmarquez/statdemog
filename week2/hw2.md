# Homework Week 2 Stat 593  
### Neal Marquez 


1) Consider a one-sex closed population with 3 age groups. Their initial
population, age- specific fertility rates (in expected births per individual in
the next period) and survivalprobabilities are as follows:  

Sections  

    a) Find the crude birth rate for this population in the next period.  

The crude birth rate for the entire population is .369 births per person-year.  

    b) Find the total fertility rate for this population.  

The total fertility rate for the population that is fertile is 1.1 births per 
person year.

    c) Write down the Leslie matrix for this population.  

|    1|     2|   3+|
|----:|-----:|----:|
| 0.00| 0.675| 0.03|
| 0.65| 0.000| 0.00|
| 0.00| 0.750| 0.15|

    d) Project the population by age forward one period.  

`[11895, 11700, 14850]`  

    e) Project the population by age forward 10 periods.  

`[392.2886, 366.7122, 524.3282]`  

    f) Find the crude birth rate and TFR for this population 10 periods into the future.  

The crude birth rate will decrease to 0.3388895 while the TFR will remqain the 
same because it is independent of age distribution. 

    g) Find the asymptotic rate of increase of the population.  

The population is decreasing at a rate of 0.3815704 every annual cycle.  

    h) Find the stable age distribution of the population.  

The stable age distribution is `[0.3037604, 0.2891735, 0.4070661]`.  

    i) Find the reproductive value of individuals in each age group.  

The relative reproductive values are `[1.00000000, 1.05044361, 0.05630754]`  

2. This question continues Question 2 of Homework 1.  

Sections  

    a) From the UN’s 2015 World Population Prospects, extract the UN estimates 
    of the age-specific fertility rates in Mexico in 2005-2010, and of the 
    female population by age in 2005-2010.  
    
See code attached below.  

    b) For each age group, calculate 5F~x, the expected number of live female 
    births perwoman per five-year period.  
    
The expected number of live female births for each age group are  
```
(0.00, 0.00, 0.08, 0.24, 0.32, 0.26, 0.16, 0.07, 0.01, 0.00, 0.00, 0.00, 
0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
```  

    c) Using these numbers and the mortality rates from Homework 1, form and 
    write out the Leslie matrix for this population.  
    
The non zero values are 

```
(row, col): value
(1,3): 0.086657
(1,4): 0.250306
(1,5): 0.326896
(1,6): 0.269822
(1,7): 0.163942
(1,8): 0.069453
(1,9): 0.014382
(1,10): 0.002371
(2,1): 0.978062
(3,2): 0.998353
(4,3): 0.998709
(5,4): 0.998288
(6,5): 0.997652
(7,6): 0.996973
(8,7): 0.995980
(9,8): 0.994295
(10,9): 0.991438
(11,10): 0.986806
(12,11): 0.979396
(13,12): 0.967662
(14,13): 0.949345
(15,14): 0.920998
(16,15): 0.877819
(17,16): 0.813603
(18,17): 0.721625
(19,18): 0.597372
(20,19): 0.444033
(21,20): 0.275736
(21,21): 0.077760
```

    d)  Assuming that fertility and mortality rates stay constant over time, 
    and that net migration is zero at all ages, project the population forward 
    one period from 2005, to 2010.  
    
See next plot.  

    e)  Under the same assumptions, project the population forward 15 years, to 2020
    
![](/home/nmarquez/Documents/Classes/statdemog/week2/compareyear.jpg "")  

    g) Assuming that age-specific net migration rates will stay constant to 
    2020, project the Mexican population forward from 2005 to 2010 and to 2020.

![](/home/nmarquez/Documents/Classes/statdemog/week2/compareyear2.jpg "")  

    h)  Compare your population projections for 2010 and 2020 with and without 
    migration.
    
![](/home/nmarquez/Documents/Classes/statdemog/week2/comparemig.jpg "")  

## Appendix  

```{R}
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

```