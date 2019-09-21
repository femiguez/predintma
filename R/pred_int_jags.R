#' Prediction intervals using JAGS
#' 
#' @title pred_int_jags
#' @param formula standard R formula response ~ trial
#' @param data data frame
#' @param df degrees of freedom for the t-distribution. If 'NA' it will be treated as a nuisance parameter
#' @param level probability level
#' @param quiet whether to display progress of sampling
#' @param return type of information to return either 'pdi' just the prediction interval or 'all' which include the output from the model and Gelman diagnostics
#' @return prediction interval (or more detail)
#' @export
#' 

pred_int_jags <- function(formula, data, df=1, level = 0.95, 
                          quiet = FALSE, return = c("pdi","all")){
  
  
  return <- match.arg(return)
  
  var.names <- all.vars(as.formula(formula))
  dat <- data[,var.names]
  names(dat)  <- c("y","trial")
  
  alpha <- 1 - level
  alpha.d2 <- alpha/2
  lowq <- alpha.d2
  upperq <- 1 - alpha.d2
  
  progress.bar <- "text"
  if(quiet) progress.bar <- "none"
  
  datl <- list(y = dat$y, trial = dat$trial,
              t.n = length(unique(dat$trial)), 
              n = nrow(dat), df = df)
  
  init <- list(mu = 0, tau = 0.01, tau.trial = 0.01)
  
  if(!is.na(df)){
    modelstring="
    model {
    # Single intercept model likelihood
    for (i in 1:n) {
    y[i]~dnorm(mu + b[trial[i]],tau)
    }
    
    # trial effect 
    for(j in 1:t.n){
    b[j] ~ dt(0, tau.trial, df)
    }
    
    # priors
    mu ~ dnorm(0,0.00001) # intercept prior
    tau ~ dgamma(0.0001,0.0001) ## tau is the residual precision
    sigma <- 1.0/sqrt(tau)
    tau.trial ~ dgamma(0.0001, 0.0001)
    sigma.trial <- 1.0/sqrt(tau.trial)
    # generate predictions 
    beff ~ dt(0, tau.trial,df)
    pred <- mu + beff
  }"
    mdl=jags.model(textConnection(modelstring), data=datl, inits=init, 
                   n.chains = 2, quiet = quiet)
    update(mdl,n.iter = 15000, n.burnin=10000, progress.bar = progress.bar)
    output=coda.samples(model=mdl, 
                        variable.names=c("mu","pred","sigma","sigma.trial"), 
                        n.iter=30000, thin=10, progress.bar = progress.bar)
    geld <- gelman.diag(output)
    mc1 <- output[[2]]
    pdi <- quantile(as.data.frame(mc1)$pred, probs = c(0.5, lowq, upperq))
  }
  
  if(is.na(df)){
    modelstring="
    model {
    # Single intercept model likelihood
    for (i in 1:n) {
    y[i]~dnorm(mu + b[trial[i]],tau)
    }
    
    # trial effect 
    for(j in 1:t.n){
    b[j] ~ dt(0, tau.trial, df)
    }
    
    # priors
    mu ~ dnorm(0,0.00001) # intercept prior
    tau ~ dgamma(0.0001,0.0001) ## tau is the residual precision
    df ~ dunif(1,10)
    sigma <- 1.0/sqrt(tau)
    tau.trial ~ dgamma(0.0001, 0.0001)
    sigma.trial <- 1.0/sqrt(tau.trial)
    # generate predictions 
    beff ~ dt(0, tau.trial,df)
    pred <- mu + beff
  }"
    mdl=jags.model(textConnection(modelstring), data=datl, inits=init, 
                   n.chains = 2, quiet = quiet)
    update(mdl,n.iter = 15000, n.burnin=10000, quiet = quiet)
    output=coda.samples(model=mdl, 
                        variable.names=c("mu","pred","sigma","sigma.trial","df"), 
                        n.iter=30000, thin=10, quiet = quiet)
    geld <- gelman.diag(output)
    mc1 <- output[[2]]
    pdi <- quantile(as.data.frame(mc1)$pred, probs = c(0.5, lowq, upperq))
  }
  if(return == "pdi"){
    ans <- pdi
  }
  if(return == "all"){
    ans <- list(pdi=pdi, output = output, geld = geld)
  }
  return(ans)
}