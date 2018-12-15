#' Fit Bayesian Dose Response Model Using rJAGS
#'
#' @param dat trial data_frame
#' @param modtext BUGS formatted text of model compatibale with analysis data
#' @param modvars vector of model variables to sample from model
#' @param burn number of burn in samples
#' @param mcmcSamps number of MCMC samples for inference
#' @param chains number of sample chains
#' @param type response data type c("bin","cont")
#' @param trace_plots indicator if traceplots should be output
#' @return MCMC samples (as MCMC list) from modtext on modvars
#' @examples
#' Fit bayesian sigmoid EMAX model:
#'
#' fitmod(sample_data,modtext="sigEmax",modvars="dr_prob",burn=10000,mcmcSamps=10000,chains=3)
#'
#' @export
#' @import rjags
#'
fitmod <- function(dat,
                   modtext,
                   modvars=c("prob"),
                   burn,
                   mcmcSamps=1000,
                   chains=1,
                   type = "bin",
                   trace_plots = FALSE) {
  init_seed <- list(.RNG.seed = round(runif(1, 1, 10000)), .RNG.name = "base::Wichmann-Hill")
  model <- suppressWarnings(jags.model(textConnection(modtext), data = dat_jags, n.chains = chains,
                                       n.adapt = 1000, quiet = TRUE))
  update(model, burn = burn, progress.bar = "none")
  samps <- coda.samples(model, n.iter = mcmcSamps, variable.names = modvars, progress.bar = "none")
  if (trace_plots) plot(samps)
  return(samps)
}

#' Creates BUGS model for linear dose combination model
#'
#' @param data_type response data type, args = c("binomial")
#' @examples
#' Dichotomous linear model:
#' linear()
#'
#' @export
#'
linear <- function(data_type){
  if(data_type == "binomial"){
    mod_text <- "
    model{
      for(i in 1:N){
        y[i]~dbin(p[i],n[i])
        logit(p[i]) <- theta0+theta1*d1[i]+theta2*d2[i]+theta12*d1[i]*d2[i]
      }
      theta0 ~ dnorm(mu0,1/sd0^2)
      theta1 ~ dnorm(mu1,1/sd1^2)
      theta2 ~ dnorm(mu2,1/sd2^2)
      theta12 ~ dnorm(mu3,1/sd3^2)
      # probs
      for(i in 1:N){
        prob[i] <- ilogit(theta0+theta1*d1[i]+theta2*d2[i]+theta12*d1[i]*d2[i])
      }
    }
    "
  }
}
