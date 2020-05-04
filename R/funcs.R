##' Get parameters for Beta distribution
##'
##' This function finds the two shape parameters for the Beta distribution of a random variable between 0 and 1.
##' Note that the function uses a different method depending on the arguments supplied. The three methods are:
##' \enumerate{
##' \item When the mean (\code{mu}) and variance (\code{sigma}) are supplied, the solution is found analytically.
##' \item When observed probabilities (\code{probs}) at each quantile (\code{quantiles}) are given, the
##' solution is found by minimizing the Sum of the Squared Errors (SSE) using the Nelder-Mead optimization
##' algorithm. Note that the fitting algorithm performs best when the five standard quantiles are supplied (Min, 25th, Median, 75th, Max).
##' \item When only observed probabilities (\code{probs}) are supplied, the function uses Maximum Likelihood Estimation (MLE).
##' }
##'
##' @param mu scalar giving the mean \eqn{\mu}
##' @param sigma scalar giving the variance \eqn{\sigma^2}
##' @param quantiles vector of quantiles for which proportions are observed. Expects: c('min', '25th', 'median', '75th', 'max').
##' @param probs vector of observed probabilities or proportions
##'
##' @return A list containing the two shape parameters of the Beta distribution
##'
##' @author John Giles
##'
##' @example R/examples/get_beta_params.R
##'
##' @family beta_params
##'
##' @export
##'


get_beta_params <- function(
  mu=NULL,
  sigma=NULL,
  quantiles=NULL,
  probs=NULL
) {

  if (all(!is.null(mu), !is.null(sigma), is.null(quantiles), is.null(probs))) {

    message('Calculating Beta distribution parameters analytically from mean (mu) and variance (sigma)')

    shape1 <- ((1-mu) / sigma - 1/mu) * mu^2

    return(list(shape1=shape1, shape2=shape1 * (1 / mu-1)))

  } else if (all(is.null(mu), is.null(sigma), !is.null(quantiles), !is.null(probs))) {

    if(!(identical(length(quantiles), length(probs)))) {
      stop("Dimensions of 'quant' and 'probs' arguments must match.")
    }

    message('Calculating Beta distribution parameters from quantiles using sums of squares')

    fit_beta <- function(x,           # vector of shape and rate parameters for beta distribution
                         quantiles,   # the quantiles for which proportions are observed
                         probs        # the observed proportions
    ) {
      sum((qbeta(quantiles, x[1], x[2]) - probs)^2)
    }

    suppressWarnings(
      params <- optim(par=c(1,1), # initialize with flat beta distribution
                      fn=fit_beta,
                      quantiles=quantiles,
                      probs=probs,
                      method='Nelder-Mead',
                      control=list(maxit=1e5))$par
    )

    return(list(shape1=params[1], shape2=params[2]))

  } else if (all(is.null(mu), is.null(sigma), is.null(quantiles), !is.null(probs))) {

    message('Calculating Beta distribution parameters from probabilities using maximum likelihood')
    message(paste('n =', length(probs), sep=' '))

    suppressWarnings(
      params <- MASS::fitdistr(probs, 'beta', list(shape1=2, shape2=2), control=list(maxit=1e5))
    )

    return(list(shape1=as.numeric(params$estimate[1]),
                shape2=as.numeric(params$estimate[2])))

  } else {

    stop('Arguments must be only: mu & sigma | quantiles & probs | probs only')
  }
}


##' Get Beta parameters for age groups
##'
##' This function finds the two shape parameters for the Beta distribution for aggregated age groups given individual-level
##' vaccination data. The individual-level data is comprised of the age of the individual (\code{age}) and a binary indicator of vaccination status
##' (\code{vacc}). The function first uses a Binomial GAM to estimate the mean (mu) and variance (sigma) of the proportion vaccinated for each
##' of the defined age groups and then finds the shape parameters of the Beta distribution analytically using \code{\link{get_beta_params}}.
##'
##' Note that a Binomial GLM is used if the number of age groups is 2 or less.
##'
##' @param age a vector giving the age at vaccination of each individual
##' @param vacc a binary vector indicating the vaccination status of each individual
##' @param breaks a scalar or vector of age group breaks (passed to \code{\link{.bincode}}). If NULL (default), calculates Beta parameters for all data together.
##'
##' @return A dataframe containing the sample size, mu, sigma, and Beta distribution paramters for each age group
##'
##' @author John Giles
##'
##' @example R/examples/get_beta_params_age.R
##'
##' @family beta_params
##'
##' @export
##'


get_beta_params_age <- function(
  age=NULL,     # age of individual
  vacc=NULL,    # binary indicator of vaccination status
  breaks=NULL   # if beta params are to be calculated for age groups, give vector of breaks (see .bincode: a numeric vector of two or more cut points, sorted in increasing order)
) {

  tmp <- cbind(age, vacc)
  tmp <- tmp[complete.cases(tmp),]
  message(paste('Complete observations: n =', nrow(tmp), 'of', length(age), sep=' '))
  age <- tmp[,1]; vacc <- tmp[,2]

  if (is.null(breaks)) breaks <- c(0, ceiling(max(age)))
  if (breaks[1] != 0) breaks <- c(0, breaks)
  if (max(breaks) < max(age)) breaks <- c(breaks, ceiling(max(age)))
  breaks <- sort(breaks)

  age_group <- .bincode(age, breaks, right=FALSE)
  age_label <- factor(age_group, labels=stringr::str_c(breaks[-length(breaks)], "-", breaks[-1]))

  message(paste('Estimating mean (mu) and variance (sigma) for', nlevels(age_label), 'age groups', sep=' '))
  fit <- tryCatch({

    suppressWarnings(
      mgcv::gam(vacc ~ s(age_group, bs='tp', k=nlevels(age_label)-1), family='binomial')
    )

  }, error = function (e) {

    glm(vacc ~ age_group, family='binomial')
  })

  suppressWarnings(
    pred <- predict(fit, newdata=data.frame(age_group=sort(unique(age_group))), type='response', se.fit=TRUE)
  )

  params <- get_beta_params(mu=pred$fit, sigma=pred$se.fit)

  return(
    data.frame(age=levels(age_label),
               n=as.vector(table(age_group)),
               mu=pred$fit,
               sigma=pred$se.fit,
               shape1=params$shape1,
               shape2=params$shape2,
               row.names=NULL)
  )
}


##' Calculate vaccine doses received with routine vaccination
##'
##' This function calculates the proportion immune using conditional probabilities. The method
##' assumes that vaccination events are dependent, where individuals that have recieved the first
##' dose are the most likely to recieve the second dose and those that have received both the first
##' and second doses are the most likely to receive the third.
##'
##' When \code{v3 = NULL}, the function uses the simpler two dose method.
##'
##' @param v1 a scalar giving the proportion vaccinated with first routine immunization
##' @param v2 a scalar giving the proportion vaccinated with second routine immunization
##' @param v3 a scalar giving the proportion vaccinated with third routine immunization (default = NULL)
##'
##' @return A dataframe containing the relative proportions of the population that have recieved 0, 1, 2, or 3 doses
##'
##' @author John Giles
##'
##' @example R/examples/calc_doses.R
##'
##' @family prop_vacc
##'
##' @export
##'

calc_doses <- function(
  v1,  # Proportion vaccinated with first campaign
  v2,  # proportion vaccinated with second campaign
  v3=NULL   # proportion vaccinated with third campaign
){

  if (!all(v1 >= 0 & v1 <= 1, v2 >= 0 & v2 <= 1)) stop('Arguments must be between 0 and 1')

  if (is.null(v3)) {

    if (v2 > v1) {
      d_12 <- 0           # Dropout rate from 1 to 2
      s_12 <- v2 - v1     # Surplus of dose 2 not given to ppl who received dose 1
    } else {
      d_12 <- (v1 - v2)/v1
      s_12 <- 0
    }

    # Conditional probability terms for two doses
    p_1_2 <- v1*(1-d_12)          # Pr( v2 | v1 )
    p_1_n2 <- v1*d_12             # Pr( notv2 | v1 )

    if (v2 <= v1) { # Pr( v2 | not v1 )
      p_n1_2 <- 0
    } else {
      p_n1_2 <- (1-v1)*(v2-v1)
    }

    if (v2 <= v1) { # Pr( notv2 | notv1 )
      p_n1_n2 <- (1-v1)
    } else {
      p_n1_n2 <- (1-v1)*(1-(v2-v1))
    }

    den <- sum(p_1_2, p_1_n2, p_n1_2, p_n1_n2)

    return(data.frame(doses=0:2,
                      prop=c(p_n1_n2/den,
                             (p_1_n2 + p_n1_2)/den,
                             p_1_2/den)))

  } else if (!is.null(v3)) {

    if(!(v3 >= 0 & v3 <= 1)) stop('Arguments must be between 0 and 1')

    if (v2 > v1) {
      d_12 <- 0           # Dropout rate from 1 to 2
      s_12 <- v2 - v1     # Surplus of dose 2 not given to ppl who received dose 1
    } else {
      d_12 <- (v1 - v2)/v1
      s_12 <- 0
    }

    p_1_2 <- v1*(1-d_12)     # Pr( v2 | v1 )

    if (v3 >= p_1_2) {
      d_23 <- 0            # Dropout rate from 2 doses to 3
      s_23 <- v3 - p_1_2   # Surplus of dose 3 not given to ppl that already received 2 doses
    } else {
      d_23 <- (p_1_2 - v3)/p_1_2
      s_23 <- 0
    }

    # Conditional probability terms for three doses

    p_1_2_3 <- p_1_2*(1-d_23)
    p_1_2_n3 <- p_1_2*d_23

    if (v3 <= p_1_2) {
      p_1_n2_3 <- 0
    } else {
      p_1_n2_3 <- v1*d_12*s_23
    }

    if (v2 <= v1) {
      p_n1_2_3 <- 0
    } else if (v3 <= p_1_2) {
      p_n1_2_3 <- 0
    } else {
      p_n1_2_3 <- (1-v1)*s_12*s_23
    }

    if (v3 <= p_1_2) {
      p_1_n2_n3 <- v1*d_12
    } else {
      p_1_n2_n3 <- v1*d_12*(1-s_23)
    }

    if (v2 <= v1) {
      p_n1_2_n3 <- 0
    } else if (v3 <= p_1_2) {
      p_n1_2_n3 <- 0
    } else {
      p_n1_2_n3 <- (1-v1)*s_12*(1-s_23)
    }

    if (v3 <= p_1_2) {
      p_n1_n2_3 <- 0
    } else if (v3 > p_1_2 & v2 <= v1) {
      p_n1_n2_3 <- (1-v1)*(s_23 - v1*d_12 - (1-v1)*(v2-v1))
    } else if (v3 > p_1_2 & v2 > v1) {
      p_n1_n2_3 <- (1-v1)*s_12*(s_23 - v1*d_12 - (1-v1)*(v2-v1))
    }

    if (v3 <= p_1_2 & v2 <= v1) {
      p_n1_n2_n3 <- (1-v1)
    } else if (v3 <= p_1_2 & v2 > v1) {
      p_n1_n2_n3 <- (1-v1)*(1-s_12)
    } else if (v3 > p_1_2 & v2 <= v1) {
      p_n1_n2_n3 <- (1-v1)*(1-s_23)
    } else if (v3 > p_1_2 & v2 > v1) {
      p_n1_n2_n3 <- (1-v1)*(1-s_12)*(1-s_23)
    }

    den <- sum(p_1_2_3, p_1_2_n3, p_1_n2_3, p_n1_2_3,
               p_1_n2_n3, p_n1_2_n3, p_n1_n2_3, p_n1_n2_n3)

    return(data.frame(doses=0:3,
                      prop=c(p_n1_n2_n3/den,
                             sum(p_1_n2_n3, p_n1_2_n3, p_n1_n2_3)/den,
                             sum(p_1_2_n3, p_1_n2_3, p_n1_2_3)/den,
                             p_1_2_3/den)))
  }
}


##' Calculate vaccine doses received with routine vaccination and SIA campaign
##'
##' This function calculates the proportion immune by calculating the conditional probability of routine
##' vaccination with two- or three-dose routine coverage and one SIA campaign. The method
##' assumes that vaccination events are dependent, where individuals that have recieved the first
##' dose are the most likely to recieve the second dose and those that have received both the first
##' and second doses are the most likely to receive the third. Receipt of a dose through the SIA campaign
##' is dependent upon vaccination with at least one dose prior to the SIA campaign.
##'
##' When \code{v3 = NULL}, the function uses the simpler two dose method.
##'
##' @param v1 a scalar giving the proportion vaccinated with first routine immunization
##' @param v2 a scalar giving the proportion vaccinated with second routine immunization
##' @param v3 a scalar giving the proportion vaccinated with third routine immunization (default = NULL)
##' @param S a scalar giving the proportion vaccinated with supplemental campaign
##'
##' @return A dataframe containing the relative proportions of the population that have recieved 0, 1, 2, 3, or 4 doses
##'
##' @author John Giles
##'
##' @example R/examples/calc_doses_SIA.R
##'
##' @family prop_vacc
##'
##' @export
##'

calc_doses_SIA <- function(
  v1,        # Proportion vaccinated with first routine immunization
  v2,        # proportion vaccinated with second routine immunization
  v3=NULL,   # proportion vaccinated with third routine immunization
  S          # proportion vaccinated with supplemental campaign
){

  if (!all(v1 >= 0 & v1 <= 1, v2 >= 0 & v2 <= 1, S >= 0 & S <= 1)) stop('Arguments must be between 0 and 1')

  p_prior <- sum(calc_doses(v1=v1, v2=v2)[1:2,2])

  if (is.null(v3)) {

    if (v2 >= v1) {
      d_12 <- 0           # Dropout rate from 1 to 2
    } else {
      d_12 <- (v1 - v2)/v1
    }

    if (S >= p_prior) {
      d_S <- 0
    } else {
      d_S <- (p_prior - S)/p_prior
    }

    # Conditional probability terms for 3 doses
    p_1_2_S <- v1*(1-d_12)*(1-d_S)

    # Conditional probability terms for 2 doses
    p_1_2_nS <- v1*(1-d_12)*d_S


    p_1_n2_S <- v1*d_12*(1-d_S)
    p_1_n2_nS <- v1*d_12*d_S


    if (v2 <= v1) {

      p_n1_2_S <- 0
      p_n1_2_nS <- 0
      p_n1_n2_S <- (1-v1)*(1-d_S)
      p_n1_n2_nS <- (1-v1)*d_S

    } else {

      p_n1_2_S <- (1-v1)*(v2-v1)*(1-d_S)
      p_n1_2_nS <- (1-v1)*(v2-v1)*d_S
      p_n1_n2_S <- (1-v1)*(1-(v2-v1))*(1-d_S)
      p_n1_n2_nS <- (1-v1)*(1-(v2-v1))*d_S
    }

    den <- sum(p_1_2_S, p_1_2_nS, p_1_n2_S, p_1_n2_nS,
               p_n1_2_S, p_n1_2_nS, p_n1_n2_S, p_n1_n2_nS)

    return(data.frame(doses=0:3,
                      prop=c(p_n1_n2_nS/den,
                             sum(p_1_n2_nS, p_n1_2_nS, p_n1_n2_S)/den,
                             sum(p_1_2_nS, p_1_n2_S, p_n1_2_S)/den,
                             p_1_2_S/den)))

  } else if (!is.null(v3)) {

    if(!(v3 >= 0 & v3 <= 1)) stop('Arguments must be between 0 and 1')

    p_prior <- sum(calc_doses(v1=v1, v2=v2, v3=v3)[1:3,2])

    if (v2 >= v1) {
      d_12 <- 0           # Dropout rate from 1 to 2
      s_12 <- v2 - v1     # Surplus of dose 2 not given to ppl who received dose 1
    } else {
      d_12 <- (v1 - v2)/v1
      s_12 <- 0
    }

    p_1_2 <- v1*(1-d_12)     # Pr( v2 | v1 )

    if (v3 >= p_1_2) {
      d_23 <- 0            # Dropout rate from 2 doses to 3
      s_23 <- v3 - p_1_2   # Surplus of dose 3 not given to ppl that already received 2 doses
    } else {
      d_23 <- (p_1_2 - v3)/p_1_2
      s_23 <- 0
    }

    if (S >= p_prior) {
      d_S <- 0              # Dropout rate from routine immunization to SIA
    } else {
      d_S <- (p_prior - S)/p_prior
    }

    # Possible ways to receive 4 doses

    p_1_2_3_S <- p_1_2*(1-d_23)*(1-d_S)

    # Possible ways to receive 3 doses

    p_1_2_3_nS <- p_1_2*(1-d_23)*d_S
    p_1_2_n3_S <- p_1_2*d_23*(1-d_S)

    if (v3 <= p_1_2) {
      p_1_n2_3_S <- 0

    } else {
      p_1_n2_3_S <- v1*d_12*s_23*(1-d_S)
    }

    if (v2 <= v1) {
      p_n1_2_3_S <- 0
    } else if (v3 <= p_1_2) {
      p_n1_2_3_S <- 0
    } else {
      p_n1_2_3_S <- (1-v1)*s_12*s_23*(1-d_S)
    }

    # Possible ways to receive 2 doses

    p_1_2_n3_nS <- p_1_2*d_23*d_S

    if (v3 <= p_1_2) {
      p_1_n2_3_nS <- 0

    } else {
      p_1_n2_3_nS <- v1*d_12*s_23*d_S
    }

    if (v2 <= v1) {
      p_n1_2_3_nS <- 0
    } else if (v3 <= p_1_2) {
      p_n1_2_3_nS <- 0
    } else {
      p_n1_2_3_nS <- (1-v1)*s_12*s_23*d_S
    }

    if (v3 <= p_1_2) {
      p_1_n2_n3_S <- v1*d_12*(1-d_S)
    } else {
      p_1_n2_n3_S <- v1*d_12*(1-s_23)*(1-d_S)
    }

    if (v2 <= v1) {
      p_n1_2_n3_S <- 0
    } else if (v3 <= p_1_2) {
      p_n1_2_n3_S <- 0
    } else {
      p_n1_2_n3_S <- (1-v1)*s_12*(1-s_23)*(1-d_S)
    }

    if (v3 <= p_1_2) {
      p_n1_n2_3_S <- 0
    } else if (v3 > p_1_2 & v2 <= v1) {
      p_n1_n2_3_S <- (1-v1)*(s_23 - v1*d_12 - (1-v1)*(v2-v1)) * (1-d_S)
    } else if (v3 > p_1_2 & v2 > v1) {
      p_n1_n2_3_S <- (1-v1)*s_12*(s_23 - v1*d_12 - (1-v1)*(v2-v1)) * (1-d_S)
    }

    # Possible ways to receive 1 dose

    if (v3 <= p_1_2) {
      p_1_n2_n3_nS <- v1*d_12*d_S
    } else {
      p_1_n2_n3_nS <- v1*d_12*(1-s_23)*d_S
    }

    if (v2 <= v1) {
      p_n1_2_n3_nS <- 0
    } else if (v3 <= p_1_2) {
      p_n1_2_n3_nS <- 0
    } else {
      p_n1_2_n3_nS <- (1-v1)*s_12*(1-s_23)*d_S
    }

    if (v3 <= p_1_2) {
      p_n1_n2_3_nS <- 0
    } else if (v3 > p_1_2 & v2 <= v1) {
      p_n1_n2_3_nS <- (1-v1)*(s_23 - v1*d_12 - (1-v1)*(v2-v1)) * d_S
    } else if (v3 > p_1_2 & v2 > v1) {
      p_n1_n2_3_nS <- (1-v1)*s_12*(s_23 - v1*d_12 - (1-v1)*(v2-v1)) * d_S
    }

    if (v3 <= p_1_2 & v2 <= v1) {
      p_n1_n2_n3_S <- (1-v1)*(1-d_S)
    } else if (v3 <= p_1_2 & v2 > v1) {
      p_n1_n2_n3_S <- (1-v1)*(1-s_12)*(1-d_S)
    } else if (v3 > p_1_2 & v2 <= v1) {
      p_n1_n2_n3_S <- (1-v1)*(1-s_23)*(1-d_S)
    } else if (v3 > p_1_2 & v2 > v1) {
      p_n1_n2_n3_S <- (1-v1)*(1-s_12)*(1-s_23)*(1-d_S)
    }

    # Possible ways to receive zero doses

    if (v3 <= p_1_2 & v2 <= v1) {
      p_n1_n2_n3_nS <- (1-v1)*d_S
    } else if (v3 <= p_1_2 & v2 > v1) {
      p_n1_n2_n3_nS <- (1-v1)*(1-s_12)*d_S
    } else if (v3 > p_1_2 & v2 <= v1) {
      p_n1_n2_n3_nS <- (1-v1)*(1-s_23)*d_S
    } else if (v3 > p_1_2 & v2 > v1) {
      p_n1_n2_n3_nS <- (1-v1)*(1-s_12)*(1-s_23)*d_S
    }

    den <- sum(p_1_2_3_S, p_1_2_n3_S, p_1_n2_3_S, p_n1_2_3_S,
               p_1_n2_n3_S, p_n1_2_n3_S, p_n1_n2_3_S, p_n1_n2_n3_S,
               p_1_2_3_nS, p_1_2_n3_nS, p_1_n2_3_nS, p_n1_2_3_nS,
               p_1_n2_n3_nS, p_n1_2_n3_nS, p_n1_n2_3_nS, p_n1_n2_n3_nS)

    return(data.frame(doses=0:4,
                      prop=c(p_n1_n2_n3_nS/den,
                             sum(p_1_n2_n3_nS, p_n1_2_n3_nS, p_n1_n2_3_nS, p_n1_n2_n3_S)/den,
                             sum(p_1_2_n3_nS, p_1_n2_3_nS, p_n1_2_3_nS,
                                 p_1_n2_n3_S, p_n1_2_n3_S, p_n1_n2_3_S)/den,
                             sum(p_n1_2_3_S, p_1_n2_3_S, p_1_2_n3_S, p_1_2_3_nS)/den,
                             p_1_2_3_S/den)))
  }
}



##' Calculate the total proportion vaccinated by routine vaccination
##'
##' This function calculates the proportion immune due to vaccination given the proportion vaccinated at
##' each routine activity and the effectiveness of each dose. The function assumes that vaccination events
##' are dependent by default, where individuals that have recieved the first
##' dose are the most likely to recieve the second dose and those that have received both the first
##' and second doses are the most likely to receive the third. The function uses either the two-dose
##' or three-dose method based on the length of \code{V}.
##'
## Under the assumption of dependence, the total
## proportion vaccinated is:
## \deqn{
## p_{\text{vacc}} = \sum_j \text{efficacy}_j \times \Pr(\text{dose}_j)
## }{
##  p_vacc = \sum_j effectiveness_j * Pr(dose_j)
## }
## Under the assumption of independence, the total proportion vaccinated is:
##  \deqn{
## p_{\text{vacc}} = 1 - ( \prod_j effectiveness_j * \Pr(not vaccinated by V_j) )
## }{
##  p_vacc = 1 - ( \prod_j effectiveness_j * Pr(not vaccinated by V_j) )
## }
##'
##' @param V a vector giving the proportion vaccinated for up to three routine immunization activities
##' @param effectiveness scalar or vector giving the vaccine effectiveness for each number of doses
##' @param independent logical indicating if receipt of routine vaccine dose is depends on the number
##' of prior doses received (default = FALSE)
##'
##' @return A scalar giving the total proportion of the population immune due to vaccination
##'
##' @author John Giles
##'
##' @example R/examples/calc_prop_vacc.R
##'
##' @family prop_vacc
##'
##' @export
##'

calc_prop_vacc <- function(
  V,
  effectiveness,
  independent=FALSE
){

  if (!(length(V) == length(effectiveness))) stop('length of V and effectiveness must be equal')

  if (independent == FALSE) {

    if (length(V) == 2) { # Two dose vaccine

      p <- calc_doses(v1=V[1], v2=V[2])

      return(sum(c(0, effectiveness) * p[,2]))

    } else if (length(V) == 3) { # 3 dose vaccine

      p <- calc_doses(v1=V[1], v2=V[2], v3=V[3])

      return(sum(c(0, effectiveness) * p[,2]))

    }
  } else if (independent == TRUE) {

    if (length(V) == 2) {

      return(1 - ((1-effectiveness[1]*V[1]) * (1-effectiveness[2]*V[2])))

    } else if (length(V) == 3) { # 3 dose vaccine

      return(1 - ((1-effectiveness[1]*V[1]) * (1-effectiveness[2]*V[2]) * (1-effectiveness[3]*V[3])) )

    }
  }
}



##' Calculate the total proportion vaccinated by routine vaccination and SIA campaign
##'
##' This function calculates the proportion immune due to vaccination given the proportion vaccinated at
##' each routine activity plus one SIA campaign and the effectiveness of each dose. The function assumes
##' that vaccination events are dependent by default, where individuals that have recieved the first
##' dose are the most likely to recieve the second dose and those that have received both the first
##' and second doses are the most likely to receive the third. receipt of dose in SIA campaign is dependent on
##' having any number of prior doses. The function uses either the two-dose plus SIA
##' or three-dose plus SIA method based on the length of \code{V}.
##'
## Under the assumption of dependence, the total
## proportion vaccinated is:
## \deqn{
## p_{\text{vacc}} = \sum_{j} \text{efficacy}_j \times \Pr(\text{dose}_j)
## }{
##  p_vacc = \sum_{j} effectiveness_j * Pr(dose_j)
## }
## Under the assumption of independence, the total proportion vaccinated is:
##  \deqn{
## p_{\text{vacc}} = 1 - ( \prod_{j} effectiveness_j * \Pr(not vaccinated by V_j) )
## }{
##  p_vacc = 1 - ( \prod_{j} effectiveness_j * Pr(not vaccinated by V_j) )
## }
##'
##' Length of effectiveness must be equal to the number of vaccination activities
##'
##' @param V a vector giving the proportion vaccinated for up to three routine immunization activities
##' @param S a scalar giving the proportion vaccinated with SIA campaign
##' @param effectiveness scalar or vector giving the vaccine effectiveness for each number of doses
##' @param independent logical indicating if receipt of routine vaccine dose is depends on the number
##' of prior doses received and receipt of SIA dose depends on at least one prior dose (default = FALSE)
##'
##' @return A scalar giving the total proportion of the population immune due to vaccination
##'
##' @author John Giles
##'
##' @example R/examples/calc_prop_vacc.R
##'
##' @family prop_vacc
##'
##' @export
##'

calc_prop_vacc_SIA <- function(
  V,
  S,
  effectiveness,
  independent=FALSE
){

  if (!(length(c(V,S)) == length(effectiveness))) stop('length of c(V,S) and effectiveness must be equal')

  if (independent == FALSE) {

    if (length(V) == 2) { # Two dose vaccine ans SIA

      p <- calc_doses_SIA(v1=V[1], v2=V[2], S=S)

      return(sum(c(0, effectiveness) * p[,2]))

    } else if (length(V) == 3) { # 3 dose vaccine and SIA

      p <- calc_doses_SIA(v1=V[1], v2=V[2], v3=V[3], S=S)

      return(sum(c(0, effectiveness) * p[,2]))

    }
  } else if (independent == TRUE) {

    if (length(V) == 2) { # 2 dose vaccine and SIA

      return(1 - ((1-effectiveness[1]*V[1]) * (1-effectiveness[2]*V[2]) * (1-effectiveness[3]*S)))

    } else if (length(V) == 3) { # 3 dose vaccine and SIA

      return(1 - ((1-effectiveness[1]*V[1]) * (1-effectiveness[2]*V[2]) * (1-effectiveness[3]*V[3]) * (1-effectiveness[4]*S)) )

    }
  }
}

