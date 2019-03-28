rm(list=ls())
library(forecast)
library(wpp2015)
library(tidyverse)

## Pull data sources from wpp2015
data(mxF)
## Only look at Mexico
DF <- subset(as_tibble(mxF), country == "Mexico", 
             select=c("age", "2005-2010"))
# get the age group names
DF$age <- as.numeric(as.character(DF$age))
# pull the mxt data for all years relevant to this analysis
# 1950-2015 in 5 year groups
mxt <- log(as.matrix(subset(mxF, country == "Mexico")[,4:16]))
row.names(mxt) <- DF$age
colnames(mxt) <- seq(1950, 2010, by=5)
ax <- apply(mxt, 1, mean) # mean within ages across time
meanadjxt <- mxt - ax # demeaned within ages
kt <- apply(meanadjxt, 2, mean)  # mean within time across ages after demeaning 

# use the linear model approach to see how to modify kt to get better age
# specific models
bx <- apply(meanadjxt, 1, function(x) lm(x ~ 0 + kt)$coefficients[[1]])
# trend of random walk
nu <- (kt[length(kt)] - kt[1]) / (length(kt) - 1)
# variance of random walk
sigma_k <- sd(kt[2:length(kt)] - kt[1:(length(kt) -1)])
# random noise component
eps <- sd(mxt - (bx %*% t(kt) + ax))

# forecast time periods 10
nF <- 10
kt_forecast <- c(kt, kt[length(kt)] + cumsum(rep(nu, nF)))
ran_var <- rep(rep(eps^2, length(kt_forecast)), nrow(DF))
for_var <- rep(
    cumsum(c(rep(0, length(kt)), rep(sigma_k^2, nF))),
    each=nrow(DF)) * rep(bx^2, length(kt_forecast))

total_err <- sqrt(ran_var + for_var)

DFforecast <- tibble(
    lnmxt=c(bx %*% t(kt_forecast) + ax), 
    age=rep(DF$age, length(kt_forecast)),
    year=rep(seq(1950, 2010 + 5*nF, 5), each=nrow(DF)),
    model="ls") %>%
    mutate(err=total_err) %>%
    mutate(lo=lnmxt - 1.96*err, hi=lnmxt + 1.96*err)

DFforecast %>%
    mutate_at(c("lnmxt", "lo", "hi"), function(x) exp(x) * 1000) %>%
    ggplot(aes(x=year, y=lnmxt, ymin=lo, ymax=hi, group=age)) +
    geom_line(aes(color=age)) +
    geom_ribbon(aes(fill=age), alpha=.4) +
    theme_classic() +
    scale_fill_distiller(palette = "Spectral") +
    scale_color_distiller(palette = "Spectral") +
    coord_trans(y="log") +
    labs(x="Year", y="Mortality Rate per 1000", color="Age", fill="Age",
         title="Traditional Lee-Carter Estimation: Mexico") +
    scale_y_continuous(breaks=c(.1, 1,10,100, 400))

# lets create a generalized function for forecasting data with an
# arbitrary number of eigen values using the SVD approach.
# ei is the number of eigenvalues and nF is the number of forecast years

lcProj <- function(mxt, ei=1, nF=10){
    a_ <- nrow(mxt)
    y_ <- ncol(mxt)
    ages <- as.numeric(row.names(mxt))
    years <- seq(1950, 2010 + 5*nF, 5)
    lm.t.x <- t(mxt)
    Y <- sweep(lm.t.x, 2, ax)             #mean centered data
    Y.svd <- svd(Y)                       #returns U %*% diag(d) %*% t(V) = Y
    bx_ <- matrix(0, nrow=a_, ncol=ei)
    kt_ <- matrix(0, nrow=y_, ncol=ei)
    kt_forecast <- matrix(0, nrow=y_+nF, ncol=ei)
    trend_err <- matrix(0, nrow=a_, ncol=y_+nF)
    ax <- apply(mxt, 1, mean)
    pred_mxt <- sapply(1:y_, function(z) ax)
    pred_mxt_forecast <- sapply(1:(y_+nF), function(z) ax)
    for(i in 1:ei){
        bx_[,i] <- Y.svd$v[,i]
        bxsign <- sign(bx_[1,i])
        kt_[,i] <- Y.svd$d[i]*Y.svd$u[,i]
        scalefactor <- length(bx_[,i]) / sum(bx_[,i])
        kt_[,i] <- kt_[,i] / scalefactor
        bx_[,i] <- bx_[,i] * scalefactor
        pred_mxt <- pred_mxt + (bx_[,i] %*% t(kt_[,i]))
        m_ <- Arima(kt_[,i], order=c(0,1,0), include.drift=T)
        kt_forecast[,i] <- c(kt_[,i], as.vector(forecast(m_, nF)$mean))
        kt_err <- c(rep(0, y_), 1:nF * m_$sigma2)
        trend_err <- trend_err + c(bx_[,i]^2 %*% t(kt_err))
        pred_mxt_forecast <- pred_mxt_forecast + 
            (bx_[,i] %*% t(kt_forecast[,i]))
    }
    
    bxDF <- tibble(
        age=rep(ages, ei),
        bx=rep(1:ei, each=a_),
        value=c(bx_)) %>%
        spread("bx", "value", sep="_")
    ktDF <- tibble(
        year=rep(years, ei),
        kt=rep(1:ei, each=y_+nF),
        value=c(kt_forecast)) %>%
        spread("kt", "value", sep="_")
    
    total_err <- (trend_err + sd(mxt - pred_mxt)^2)^.5
    
    svdDFforecast <- tibble(
        lnmxt=c(pred_mxt_forecast), 
        age=rep(ages, y_+nF),
        year=rep(years, each=a_),
        model=paste0("svd: ", ei)) %>%
        mutate(err=c(total_err)) %>%
        mutate(lo=lnmxt - 1.96*err, hi=lnmxt + 1.96*err) %>%
        left_join(bxDF, by="age") %>%
        left_join(ktDF, by="year")
    
    svdDFforecast
}

# lets try forecasting with more than one eigen values and see the differnce
# in the forecasts
svdDFforecast <- bind_rows(lapply(1:4, function(i){
    lcProj(mxt, i)}))

svdDFforecast %>%
    mutate_at(c("lnmxt", "lo", "hi"), function(x) exp(x) * 1000) %>%
    ggplot(aes(x=year, y=lnmxt, ymin=lo, ymax=hi, group=age)) +
    geom_line(aes(color=age)) +
    geom_ribbon(aes(fill=age), alpha=.4) +
    theme_classic() +
    scale_fill_distiller(palette = "Spectral") +
    scale_color_distiller(palette = "Spectral") +
    coord_trans(y="log") +
    labs(x="Year", y="Mortality Rate per 1000", color="Age", fill="Age",
         title="SVD Lee-Carter Estimation: Mexico") +
    scale_y_continuous(breaks=c(.1, 1,10,100, 400)) +
    facet_wrap(~model)

svdDFforecast %>%
    select_at(c("age", grep("bx_", names(.), value = T))) %>%
    na.omit() %>%
    unique() %>%
    gather("bx", "value", -age) %>%
    ggplot(aes(x=age, y=value)) +
    geom_point() +
    theme_classic() +
    facet_wrap(~bx, scales="free_y") +
    labs(x="Age", y="Bx", title="Age Effects for First Four Eigen Values")

# Note that how we forecast the kt's is a modelers choice. Right now we are 
# forecasting each separately using a RW w/drift is maybe not the best choice
# but its sufficient here

svdDFforecast %>%
    select_at(c("year", grep("kt_", names(.), value = T))) %>%
    na.omit() %>%
    unique() %>%
    gather("kt", "value", -year) %>%
    mutate(Forecast=year>2010) %>%
    ggplot(aes(x=year, y=value, color=Forecast)) +
    geom_point() +
    theme_classic() +
    facet_wrap(~kt, scales="free_y") +
    labs(x="Year", y="Kt", title="Time Effects for First Four Eigen Values")
