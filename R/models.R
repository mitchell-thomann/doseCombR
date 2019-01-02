#' Generate random dose-response data from monotherapy sigmoid EMAX model
#'
#' @param type data type, only supports "binomial" currently
#' @param dose vector of monotherapy investigational doses
#' @param n vector of sample sizes by dose
#' @param e0 placebo effect in EMAX model, logit scale for binomial data
#' @param emax emax effect in EMAX model, logit scale for binomial data
#' @param ed50 ed50 parameter in EMAX model
#' @param hill slope parameter in EMAX model
#'
#' @return tibble of dose response data & plot of true & observed effect
#'
#' @importFrom dplyr data_frame
#' @import ggplot2
#'
gendr_emax <- function(type = "binomial",
                       dose = c(0, 2, 5, 10, 15, 30),
                       n = rep(50, 6),
                       e0 = -2,
                       emax = 2,
                       ed50 = 7.5,
                       hill = 2) {
  emax_dr <- e0 + emax * dose^hill / (ed50^hill + dose^hill)
  if (type == "binomial") {
    emax_prob <- exp(emax_dr) / (1 + exp(emax_dr))
    dat_true <- data_frame(dose = dose, resp = emax_prob)
    resp <- rbinom(n = length(n), size = n, prob = emax_prob)
    dat_dr <- data_frame(dose = dose, n = n, resp = resp, respmean = resp / n)
  }
  p <- ggplot(data = dat_true, aes(x = dose, y = resp)) + geom_point() + geom_line() +
    geom_point(data = dat_dr, aes(x = dose, y = respmean), colour = "red") +
    ylab("Response") + xlab("Dose") + ylim(0, 1) + scale_x_continuous(breaks = dose)
  return(list(data = dat_dr, plot = p))
}

#' Linear combination function
#' 
#' @param type data type, "binomial" or "continuous"
#' @param dose_1 combination dose 1
#' @param dose_2 combination dose 2
#' @param e0 placebo effect, logit scale for binomial data
#' @param theta1 monotherapy slope for dose1, logit scale for binomial data
#' @param theta2 monotherapy slope for dose 2, logit scale for binomial data
#' @param theta12 interaction term, < 0 for antagonism, 0 for additivity, >0 for synergism
#' 
#' @return single combination probability
#' 
linearcomb <- function(type, dose_1, dose_2, e0, theta1, theta2, theta12){
  lc <- e0 + dose_1*theta1 + dose_2*theta2 + dose_1*dose_2*theta12
  if(type == "binomial") lc <- exp(lc)/(1+exp(lc))
  return(lc)
}

#' Generate random combination dose response data using linear models
#'
#' @param type data type, only supports "binomial" currently
#' @param dose_1 vector of monotherapy investigational doses for drug 1
#' @param dose_2 vector of monotherapy investigational doses for drug 2
#' @param n_total total number of subjects, (currently) assumes equivalent number of subjects per 
#'   dose combination
#' @param e0 placebo effect, logit scale for binomial data
#' @param theta1 monotherapy slope for dose1, logit scale for binomial data
#' @param theta2 monotherapy slope for dose 2, logit scale for binomial data
#' @param theta12 interaction term, < 0 for antagonism, 0 for additivity, >0 for synergism
#'
#' @return tibble of dose response data
#'
#' @importFrom dplyr data_frame
#' @importFrom tidyr gather
#' @import plotly
#'
gendr_linearcomb <- function(type = "binomial",
                             dose_1 = c(0, 2, 5, 10, 15, 30),
                             dose_2 = c(0, 1, 2, 4),
                             n_total = 240,
                             e0 = -2.5,
                             theta1 = 0.25,
                             theta2 = 0.5,
                             theta12 = 0) {
  dat_doses <- data_frame(dose1 = rep(dose_1, length(dose_2)),
                          dose2 = rep(dose_2, each = length(dose_1)))
  dat_doses$n <- rep(n_total/(length(dose_1)*length(dose_2)), nrow(dat_doses))
  dat_dr <- dat_doses %>% mutate(true_probability = linearcomb(dose_1 = dose1, dose_2 = dose2,
                                                   e0 = e0, theta1 = theta1, theta2 = theta2,
                                                   theta12 = theta12, type = "binomial"),
                                 resp = rbinom(size = n, prob = probability, n = n()),
                                 observed_probability = resp/n)
  dat_plot <- dat_dr %>% gather("true_probability", "observed_probability", key = "group", value = "mean")
  p <- plot_ly(dat_plot, x = ~dose1, y = ~dose2, z = ~mean, color = ~group) %>% add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Dose 1'),
                        yaxis = list(title = 'Dose 2'),
                        zaxis = list(title = 'Value')))
  return(list(data = dat_dr, plot = p))
}

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
#'
#' @return MCMC samples (as MCMC list) from modtext on modvars
#'
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
                   modvars = c("prob"),
                   burn,
                   mcmcSamps = 1000,
                   chains = 1,
                   type = "bin",
                   trace_plots = FALSE) {
  ## MT: add call to load model text in bugs
  init_seed <- list(.RNG.seed = round(runif(1, 1, 10000)), .RNG.name = "base::Wichmann-Hill")
  model <- suppressWarnings(jags.model(textConnection(modtext),
    data = dat_jags, n.chains = chains,
    n.adapt = 1000, quiet = TRUE
  ))
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
                   prior = list(
                     theta0 = list(pmean, psd),
                     theta1 = list(pmean, psd),
                     theta2 = list(pmean, psd),
                     theta12 = list(pmean, psd)
                   )) {
  if (data_type == "binomial") {
    mod_text <- paste0("
    model{
      for(i in 1:N){
        y[i]~dbin(p[i],n[i])
        logit(p[i]) <- theta0+theta1*d1[i]+theta2*d2[i]+theta12*d1[i]*d2[i]
      }
      theta0 ~ dnorm(", theta0$pmean, ",1/", theta0$psd, "^2)
      theta1 ~ dnorm(", theta1$pmean, ",1/", theta1$psd, "^2)
      theta2 ~ dnorm(", theta2$pmean, ",1/", theta2$psd, "^2)
      theta12 ~ dnorm(", theta12$pmean, ",1/", theta12$psd, "^2)
      # probs
      for(i in 1:N){
        prob[i] <- ilogit(theta0+theta1*d1[i]+theta2*d2[i]+theta12*d1[i]*d2[i])
      }
    }
    ")
  }
  return(mod_text)
}
