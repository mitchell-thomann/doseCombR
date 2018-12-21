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
  ## MT: add call to load model text in jags
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
#' @param datatype character variable for response data type, args = c(\strong{"binomial"})
#' @param prior list of parameters in model.  For this model, must contain lists: \strong{theta0}, 
#'   \strong{theta1}, \strong{theta2}, \strong{theta12}.  Each list has values \emph{pmean} 
#'   and \emph{psd} for prior means and standard deviations respectively.
#'     
#' @examples
#' Dichotomous linear combination model:
#' linear("binomial", prior = list(theta0 = list(0,10), theta1 = list(0,10), 
#'        theta2 = list(0,10), theta12 = list(0,10)))
#'        
#' @export
#'
#'
linear <- function(datatype,
                   prior = list(theta0=list(pmean,psd),
                                theta1=list(pmean,psd),
                                theta2=list(pmean,psd),
                                theta12=list(pmean,psd))){
  if(data_type == "binomial"){
    mod_text <- paste0("
    model{
      for(i in 1:N){
        y[i]~dbin(p[i],n[i])
        logit(p[i]) <- theta0+theta1*d1[i]+theta2*d2[i]+theta12*d1[i]*d2[i]
      }
      theta0 ~ dnorm(",theta0$pmean,",1/",theta0$psd,"^2)
      theta1 ~ dnorm(",theta1$pmean,",1/",theta1$psd,"^2)
      theta2 ~ dnorm(",theta2$pmean,",1/",theta2$psd,"^2)
      theta12 ~ dnorm(",theta12$pmean,",1/",theta12$psd,"^2)
      # probs
      for(i in 1:N){
        prob[i] <- ilogit(theta0+theta1*d1[i]+theta2*d2[i]+theta12*d1[i]*d2[i])
      }
    }
    ")
  }
  return(mod_text)
}
