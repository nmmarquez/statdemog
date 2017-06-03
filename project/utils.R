pacman::p_load(INLA, ar.matrix, data.table, ggplot2)

create_obs_mesh <- function(n, sigma0, range0){
    loc <- matrix(runif(n*2), n, 2) # simulate observed points
    mesh <- inla.mesh.create(loc, refine=list(max.edge=1.)) # create mesh
    mesh$loc.obs <- loc
    mesh$proj <- inla.mesh.projector(mesh)
    mesh$spde <- inla.spde2.matern(mesh) 
    mesh$kappa0 <- sqrt(8)/range0 # inla paramter transform
    mesh$tau0 <- 1/(sqrt(4*pi)*mesh$kappa0*sigma0) # inla parameter transform
    mesh$Q <- inla.spde2.precision(mesh$spde, 
                                   theta=c(log(mesh$tau0), log(mesh$kappa0)))
    return(mesh)
}

plot_mesh <- function(mesh){
    plot(mesh)
    points(mesh$loc.obs[,1], mesh$loc.obs[,2], col="red", pch=20)
}


build_data <- function(n, m, sigma0=.3, range0=1., rho=.95, cov_cor=NULL, 
                       betas=c(), X=NULL, mixed_corr=FALSE){
    if(!is.null(sigma0)){
        mesh <- create_obs_mesh(n, sigma0, range0)
        Q <- kronecker(Q.AR1(m, 1, rho), mesh$Q)
        x <- matrix(data=c(sim.AR(1, Q)), nrow=mesh$n, ncol=m)
        if(mixed_corr){
            mixset <- apply(mesh$loc[,1:2] > .25 & mesh$loc[,1:2] < .75, 1, all)
            x[mixset,] <- sapply(1:m, function(i) 
                runif(sum(mixset), range(x)[1], range(x)[2]))
        }
        ran_ <- c(x[mesh$idx$loc,])
    }
    else{
        mesh <- NULL
        x <- rep(0, n *m)
        ran_ <- x
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
    DT <- data.table(yraw=yraw)
    if(length(betas) != 0){
        for(i in 1:ncol(X)){
            DT[[paste0("X", i)]] <- X[,i]
        }
    }
    return(list(data=DT, mesh=mesh, betas=betas, reff=x, m=m))
}

plot_field <- function(data_obj){
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

data_obj <- build_data(100, 12, betas=c(1,2,3), mixed_corr=TRUE)
plot_mesh(data_obj$mesh)
plot_field(data_obj)
