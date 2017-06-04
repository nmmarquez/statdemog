pacman::p_load(INLA, ar.matrix, data.table, ggplot2, clusterPower)

mesh_to_dt <- function(x, proj, time){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), time=time, 
                     obs=c(inla.mesh.project(proj, field=x)))
    DT
}

create_obs_mesh <- function(n, sigma0, range0){
    loc <- matrix(runif(n*2), n, 2) # simulate observed points
    mesh <- inla.mesh.create(loc, refine=list(max.edge=1.)) # create mesh
    mesh$loc.obs <- loc # save the observed locations
    mesh$proj <- inla.mesh.projector(mesh) # project mesh
    mesh$spde <- inla.spde2.matern(mesh) # make spde
    mesh$kappa0 <- sqrt(8)/range0 # inla paramter transform
    mesh$tau0 <- 1/(sqrt(4*pi)*mesh$kappa0*sigma0) # inla parameter transform
    mesh$Q <- inla.spde2.precision(mesh$spde, # make the precision matrix
                                   theta=c(log(mesh$tau0), log(mesh$kappa0)))
    return(mesh)
}

plot_mesh <- function(mesh){
    # plot a mesh and the points
    plot(mesh)
    points(mesh$loc.obs[,1], mesh$loc.obs[,2], col="red", pch=20)
}


build_data <- function(n, m, sigma0=.3, range0=1., rho=.95, cov_cor=NULL, 
                       betas=c(), X=NULL, mixed_corr=FALSE, p=0., miss=.2){
    sigma0_ <- ifelse(is.null(sigma0), .3, sigma0)
    mesh <- create_obs_mesh(n, sigma0_, range0) # create the mesh 
    Q <- kronecker(Q.AR1(m, 1, rho), mesh$Q) # combine time and space Q's
    x <- matrix(data=c(sim.AR(1, Q)), nrow=mesh$n, ncol=m)
    x <- x - mean(x) # ensure the simulation is zero centered
    if(mixed_corr){
        mixset <- apply(mesh$loc[,1:2] > .25 & mesh$loc[,1:2] < .75, 1, all)
        x[mixset,] <- sapply(1:m, function(i) 
            runif(sum(mixset), range(x)[1], range(x)[2]))
    }
    ran_ <- c(x[mesh$idx$loc,])
    if(is.null(sigma0)){
        x <- x * 0
        ran_ <- ran_ * 0
    }
    if(length(betas) == 0){
        betas_ <- 0
    }
    else{
        betas_ <- betas
    }
    if(is.null(X)){
        X <- sapply(1:length(betas_), function(x) runif(n*m))
        X[,1] <- 1
        if(!is.null(cov_cor)){
            X[,length(betas_)] <- rnorm(n*m, mean=cov_cor * ran_)
        }
    }
    fix_ <- rowSums(sapply(1:length(betas_), function(i) betas_[i] * X[,i]))
    yraw <- fix_ + ran_
    # obs data 
    p1 <- expit(yraw)
    yobs <- rbinom(n*m, size= 1, prob=p1)
    yobs[ rbinom(n*m, size=1, prob=p) == 1 ] <- 0
    yobs[ rbinom(n*m, size=1, prob=miss) == 1 ] <- NA
    DT <- data.table(yraw=yraw, praw=p1, obs=yobs)
    DT[,xcoord:=rep(mesh$loc.obs[,1], m)]
    DT[,ycoord:=rep(mesh$loc.obs[,2], m)]
    DT[,time:=rep(1:m, each=nrow(mesh$loc.obs))]
    if(length(betas) != 0){
        for(i in 1:ncol(X)){
            DT[[paste0("X", i)]] <- X[,i]
        }
    }
    return(list(data=DT, mesh=mesh, betas=betas, reff=x, m=m))
}

plot_field <- function(data_obj){
    # plot the true proj field from the data obj
    M <- length(data_obj$mesh$proj$x)
    dat_list <- lapply(1:data_obj$m, function(i) 
        data.table(x=rep(data_obj$mesh$proj$x, M), 
                   y=rep(data_obj$mesh$proj$y, each=M), time=i,
                   obs=c(inla.mesh.project(data_obj$mesh$proj, 
                                           field=data_obj$reff[,i]))))
    DT <- rbindlist(dat_list)
    ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8)) + 
        facet_wrap(~time)
}

run_model <- function (data_obj, family=c("binomial", "zeroinflatedbinomial1")){
    iset <- inla.spde.make.index('i', n.spde=data_obj$mesh$spde$n.spde, 
                                 n.group=data_obj$m)
    A <- inla.spde.make.A(mesh=data_obj$mesh, 
                          loc=cbind(data_obj$data$xcoord, data_obj$data$ycoord), 
                          group=data_obj$data$time) 
    
    sdat <- inla.stack(tag='stdata', data=list(y=data_obj$data$obs), 
                       A=list(A,1,1),  effects=list(iset, X2=data_obj$data$X2,
                                                    X1=data_obj$data$X1))
    h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))
    formulae <- y ~ 0 + X1 + X2 +
        f(i, model=data_obj$mesh$spde, group=i.group, 
          control.group=list(model='ar1', hyper=h.spec)) 
    prec.prior <- list(prior='pc.prec', param=c(1, 0.01))
    if(family=="binomial"){
        contfam <- list()
    }
    else{
        contfam <- list(hyper=list(theta=prec.prior))
    }
    print(system.time(res <- inla(formulae,  data=inla.stack.data(sdat),
                                  control.predictor=list(compute=TRUE, 
                                                         A=inla.stack.A(sdat),
                                                         link=1),
                                  control.family=contfam,
                                  family = family[1])))
    phat <- res$summary.fitted.values[1:nrow(data_obj$data),"mean"]
    if(!(family == "binomial")){
        vname <- "zero-probability parameter for zero-inflated binomial_1"
        phat <- phat * res$summary.hyperpar[vname, 1]
    }
    res$phat <- phat
    return(res)
}

plot_result_mesh <- function(data_obj, res){
    sets <- lapply(1:data_obj$m, function(i) 
        ((data_obj$mesh$n * (i-1)) + 1):(data_obj$mesh$n * i))
    inlalist <- lapply(1:data_obj$m, function(i) 
        mesh_to_dt(res$summary.random$i$mean[sets[[i]]], data_obj$mesh$proj, i))
    inlaDT <- rbindlist(inlalist)
    ggplot(inlaDT, aes(x, y, z=obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8)) + 
        facet_wrap(~time)
}

