rm(list=ls())
pacman::p_load(pracma, data.table, wpp2015, ggplot2, knitr)

# data cleaning
data(mxF)
DF <- subset(as.data.table(mxF), country == "Mexico", 
             select=c("age", "2005-2010"))
setnames(DF, c("age","nMx") )
str(DF)
DF[,age:=as.numeric(as.character(age))]

# plot the pdf
x <- seq(0, 10, .1)

fx <- function(x){
    exp(x**3 / -x) * x**2
}

jpeg('~/Documents/Classes/statdemog/week1/pdfX.jpg')
qplot(x, fx(x),geom = "line")
dev.off()

# hazard function
Hx <- function(x){
    (x^3) / (3)
}

# survival function
Sx <- function(x){
    exp(-Hx(x))
}

# life expectancy
e0 <- gamma(1/3) / 3^(2/3)
e0

# life expectancy at age 2
gammainc(8/3, 1/3)[["uppinc"]] / 3^(2/3) / Sx(2)

# 2a plot dat ain reg space and log space
jpeg('~/Documents/Classes/statdemog/week1/gomp.jpg')
ggplot(data=DF, aes(x=age, y=nMx)) + geom_point()
dev.off()

jpeg('~/Documents/Classes/statdemog/week1/loggomp.jpg')
ggplot(data=DF, aes(x=age, y=log(nMx))) + geom_point()
dev.off()

# Life Table code

DF[,n:=c(1,4, rep(5, nrow(DF) - 3), Inf)]
DF[,nqx:= 1-exp(-n * nMx)]
DF[age==max(age), nqx:=1]
DF[,npx:=1-nqx]
DF[,lx:=100000]
DF[,ndx:=lx*nqx]
for(i in 2:nrow(DF)){
    prevd <- DF[i-1,ndx]
    prevl <- DF[i-1,lx]
    DF[i, lx:= prevl - prevd]
    DF[,ndx:=lx*nqx]
}

DF[,nLx:=0]
currl <- DF[1,lx]
nextl <- DF[2,lx]
DF[1,nLx:=.3 * currl + .7 * nextl]
for(i in 2:(nrow(DF) - 1)){
    currl <- DF[i,lx]
    nextl <- DF[i+1,lx]
    DF[i, nLx:= n*.5*(nextl + currl)]
}

# pulled from ref table 
DF[nrow(DF),nLx:=lx * 5]

DF[,Tx:=0]
fullnLx <- DF$nLx
for(i in 1:nrow(DF)){
    DF[i,Tx:=sum(fullnLx[i:nrow(DF)])]
}

DF[,ex:=Tx/lx]

kable(subset(DF, select=-n), format="markdown", digits = 3)
