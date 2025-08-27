#' @title Regra de Encolhimento Bayesiana com Priori Logística
#' 
#' @usage logis_shrink(d, a, s, t)
#' 
#' @param d coeficientes empiricos.
#' @param a parâmetro alpha.
#' @param s desvio padrão.
#' @param tau parâmetro tau da logística.
#' 
#' @details
#' Seja \eqn{\theta} uma variável aleatória com distribuição logística, então,
#' sua densidade é dada por:
#' \deqn{g(\theta; \tau) = \frac{\exp\{-\frac{\theta}{\tau}\}}
#' {\tau(1 + \exp\{-\frac{\theta}{\tau}\})^2} \mathcal{I}_{\mathbb{R}}^{(\theta)}}
#' onde \eqn{\tau > 0} é um parâmetro da distribuição
#' 
#' @references
#' Sousa, A.R.S. (2020). Bayesian wavelet shrinkage with logistic prior.
#' Communications in Statistics - Simulation and Computation.

logis_shrink <- function(d, a, s, tau) {
  u <- rnorm(10000)
  delta <- vector(mode='double')
  for(i in 1:length(d)) {
    int1 <- mean((s*u + d[i]) * dlogis(s*u + d[i], scale=tau))
    int2 <- mean(dlogis(s*u + d[i], scale=tau))
    delta[i] <- (1-a) * int1/(a * dnorm(d[i], sd=s)/s + (1-a) * int2)
  }
  return(delta)
}
