rm(list=ls())
pacman::p_load(INLA, ar.matrix, data.table, ggplot2, clusterPower)
set.seed(123)
source("~/Documents/Classes/statdemog/project/utils.R")

# create sets of model parameters to test
model_list <- list(
    list(sigma0=NULL, betas=c(1,-2), cov_cor=NULL, p=0., mixed_corr=FALSE),
    list(sigma0=.3, betas=c(1,-2), cov_cor=NULL, p=0., mixed_corr=FALSE),
    list(sigma0=.3, betas=c(1,-2), cov_cor=NULL, p=.2, mixed_corr=FALSE),
    list(sigma0=.3, betas=c(1,-2,1), cov_cor=.8, p=.2, mixed_corr=FALSE),
    list(sigma0=.3, betas=c(1,-2), cov_cor=NULL, p=.2, mixed_corr=TRUE))

names(model_list) <- c("basic", "geo-temp", "geo-temp-zib", "spatcov", "nonstat")

# build the data sets
data_obj_list <- lapply(model_list, function(x) lapply(1:10, function(i)
    build_data(100, 12, sigma0=x$sigma0, betas=x$betas, cov_cor=x$cov_cor,
               p=x$p, mixed_corr=x$mixed_corr)))

res <- lapply(data_obj_list, function(dl) lapply(dl, function(ds)
    lapply(c("binomial", "zeroinflatedbinomial1"), function(fam)
        run_model(ds, fam))))

save(list=ls(), file="~/Documents/Classes/statdemog/project/results.Rdata")
