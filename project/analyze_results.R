rm(list=ls())
pacman::p_load(INLA, ar.matrix, data.table, ggplot2, clusterPower, knitr, pROC)
setwd("~/Documents/Classes/statdemog/project/")
load("./results.Rdata")
set.seed(123)

betaDT <- as.data.table((
    sapply(res, function(sim) sapply(sim, function(i) sapply(i, function(m)
        m$summary.fixed[2,"mean"]))) + 2)**2)
betaDT[,model:=rep(c("Binom", "ZIB"), nrow(betaDT)/2)]
kable(betaDT[,lapply(.SD,mean),by=model], "markdown")

# observe the missing values in the data
for(sim_num in 1:length(data_obj_list)){
    for (i in 1:length(data_obj_list[[1]])){
        N <- nrow(data_obj_list[[sim_num]][[i]]$data)
        if(sim_num >= 3){
            data_obj_list[[sim_num]][[i]]$data$full_obs <- 
                rbinom(N, 1, data_obj_list[[sim_num]][[i]]$data$praw*.8)
        }
        else{
            data_obj_list[[sim_num]][[i]]$data$full_obs <- 
                rbinom(N, 1, data_obj_list[[sim_num]][[i]]$data$praw)
        }
        data_obj_list[[sim_num]][[i]]$data[!is.na(obs), full_obs:=obs]
    }
}

eval_auc <- function(data_obj, rez){
    return(auc(as.integer(data_obj$data[is.na(obs), full_obs]),
               rez$phat[is.na(data_obj$data$obs)])[1])
}

eval_latent <- function(data_obj, rez){
    return(mean((c(data_obj$reff) - rez$summary.random$i$mean)**2)**.5)
}

aucDT <- as.data.table((
    sapply(1:length(data_obj_list), function(i) 
        sapply(1:length(data_obj_list[[1]]), function(j) sapply(1:2, function(m)
        eval_auc(data_obj_list[[i]][[j]], res[[i]][[j]][[m]]))))))
setnames(aucDT, names(aucDT), names(data_obj_list))
aucDT[,model:=rep(c("Binom", "ZIB"), nrow(aucDT)/2)]
kable(aucDT[,lapply(.SD,mean),by=model], "markdown")

latDT <- as.data.table((
    sapply(1:length(data_obj_list), function(i) 
        sapply(1:length(data_obj_list[[1]]), function(j) sapply(1:2, function(m)
            eval_latent(data_obj_list[[i]][[j]], res[[i]][[j]][[m]]))))))
setnames(latDT, names(latDT), names(data_obj_list))
latDT[,model:=rep(c("Binom", "ZIB"), nrow(latDT)/2)]
kable(latDT[,lapply(.SD,mean),by=model], "markdown")

jpeg("./sim1true.jpg")
plot_field(data_obj_list[[1]][[5]])
dev.off()
jpeg("./sim1binom.jpg")
plot_result_mesh(data_obj_list[[1]][[5]], res[[1]][[5]][[1]])
dev.off()
jpeg("./sim1zib.jpg")
plot_result_mesh(data_obj_list[[1]][[5]], res[[1]][[5]][[2]])
dev.off()

jpeg("./sim2true.jpg")
plot_field(data_obj_list[[2]][[5]])
dev.off()
jpeg("./sim2binom.jpg")
plot_result_mesh(data_obj_list[[2]][[5]], res[[2]][[5]][[1]])
dev.off()
jpeg("./sim2zib.jpg")
plot_result_mesh(data_obj_list[[2]][[5]], res[[2]][[5]][[2]])
dev.off()

jpeg("./sim5true.jpg")
plot_field(data_obj_list[[5]][[5]])
dev.off()
jpeg("./sim5binom.jpg")
plot_result_mesh(data_obj_list[[5]][[5]], res[[5]][[5]][[1]])
dev.off()
jpeg("./sim5zib.jpg")
plot_result_mesh(data_obj_list[[5]][[5]], res[[5]][[5]][[2]])
dev.off()

jpeg("./sim4true.jpg")
plot_field(data_obj_list[[4]][[5]])
dev.off()
jpeg("./sim4binom.jpg")
plot_result_mesh(data_obj_list[[4]][[5]], res[[4]][[5]][[1]])
dev.off()
jpeg("./sim4zib.jpg")
plot_result_mesh(data_obj_list[[4]][[5]], res[[4]][[5]][[2]])
dev.off()
