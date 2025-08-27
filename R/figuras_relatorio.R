# Figuras do relatório da FAPESP

# Packages
require(wavethresh)
require(fda.simu)
require(fda.usc)

data('tecator')

source('R/logis_shrink.R')

# Figuras -----------------------------------------------------------------

# Figura 1: Funções ondaleta e escala
{
  pdf(file='figuras/OndaletasEscala.pdf', height=5.2, width=7.8)
  par(mfrow=c(3,2), las=1, cex=0.9, mgp=c(2.6,1,0), mar=c(4,3.7,2,2.5))
  # N = 1
  draw.default(filter.number=1, family='DaubExPhase', main=expression(N==1),
               sub='', ylab=expression(psi(x)), enhance=F)
  draw.default(filter.number=1, family='DaubExPhase', main=expression(N==1),
               sub='', ylab=expression(phi(x)), enhance=F, scaling.function=T)
  # N = 2
  draw.default(filter.number=2, family='DaubExPhase', main=expression(N==2),
               sub='', ylab=expression(psi(x)), enhance=F)
  draw.default(filter.number=2, family='DaubExPhase', main=expression(N==2),
               sub='', ylab=expression(phi(x)), enhance=F, scaling.function=T)
  # N = 8
  draw.default(filter.number=8, family='DaubExPhase', main=expression(N==8),
               sub='', ylab=expression(psi(x)), enhance=F)
  draw.default(filter.number=8, family='DaubExPhase', main=expression(N==8),
               sub='', ylab=expression(phi(x)), enhance=F, scaling.function=T)
  dev.off()
}


# Figura 4: Regra de shrinkage bayesiana
{
  x <- seq(-5.95, 5.95, 0.05)
  a <- c(0.2, 0.5, 0.7, 0.9)
  t <- c(5, 20, 60, 180)
  
  ## Figura 4 (a): Variando t para alpha=0.8
  pdf(file='figuras/ShrinakgeLogistica1.pdf', height=4.3, width=5.75)
  par(mar=c(4.3, 4.3, 1.5, 1), cex=1.3, las=1)
  plot(x=1, type='n', xlab='d', ylab=expression(delta(d)),
       xlim=c(-5.5,5.5), ylim=c(-5.5,5.5))
  # abline(v=0, h=0, lty=1)
  for (i in 1:length(t))
    lines(x, logis_shrink(x, a=0.8, s=1, t[i]), col=i, lwd=2.5)
  legend('topleft',  bty='n', lwd=2, col=1:4, cex=1,
         legend=c(expression(tau == 5), expression(tau == 20),
                  expression(tau == 60), expression(tau == 180)))
  dev.off()
  
  ## Figura 4 (b): Variando alpha para t=50
  pdf(file='figuras/ShrinakgeLogistica2.pdf', height=4.3, width=5.75)
  par(mar=c(4.3, 4.3, 1.5, 1), cex=1.3, las=1)
  plot(x=1, type='n', xlab='d', ylab=expression(delta(d)),
       xlim=c(-5.5,5.5), ylim=c(-5.5,5.5))
  # abline(v=0, h=0, lty=1)
  for (i in 1:length(a))
    lines(x, logis_shrink(x, a[i], s=1, tau=50), col=i, lwd=2.5)
  legend('topleft', bty='n', lwd=2, col=1:4, cex=1,
         legend=c(expression(alpha == 0.2), expression(alpha == 0.5),
                  expression(alpha == 0.7), expression(alpha == 0.9)))
  dev.off()
}


# Figura 6: amostra, coeficientes empíricos e estimados e função recuperada
{
  ## Figura 6 (a): amostra
  pdf(file='figuras/ExemploAmostra.pdf', height=3.8, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  set.seed(282829)
  f <- DJ.EX()$bumps
  e <- rnorm(1024, 0, 7/5)
  y <- f + e
  
  plot(y, lwd=2, type='l', ylab='y', xlab='x', ylim=c(-2.5,57))
  lines(f, col='blue')
  dev.off()

  ## Figura 6 (b): coeficientes empíricos
  pdf(file='figuras/ExemploCoefEmp.pdf', height=3.8, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  d <- wd(y, filter.number=1, family='DaubExPhase')
  plot(d, scaling='by.level', main='', sub='',
       ylab='Nível de Resolução', xlab='Translação')
  dev.off()
  
  ## Figura 6 (c): coeficientes estimados
  pdf(file='figuras/ExemploCoefEst.pdf', height=3.8, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  d_hat <- threshold(d, type='soft', policy='sure')
  plot(d_hat, scaling='by.level', main='', sub='',
       ylab='Nível de Resolução', xlab='Translação')
  dev.off()
  sum(d_hat$D == 0)  # num. de coeficientes nulos

  ## Figura 6 (d): função recuperada
  pdf(file='figuras/ExemploRecup.pdf', height=3.8, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  f_hat <- wr(d_hat)
  plot(f_hat, lwd=2, type='l', ylab='y', xlab='x', ylim=c(-2.5,57))
  lines(f, col='blue')
  dev.off()
}


# Figura 7, 8 
{
  # definindo parâmetros
  set.seed(282829)
  M <- 1024
  x <- seq(1, M)/M
  alpha <- matrix(c(DJ.EX(M)$bumps, DJ.EX(M)$doppler), nrow=M)
  L <- ncol(alpha)
  I <- 10
  y1 <- runif(I)
  y <- matrix(c(y1, 1-y1), nrow=2)
  e <- matrix(rnorm(M*I, sd=7/4), nrow=M)
  A <- alpha %*% y + e
    
# Figura 7: bumps e doppler
  ## Figura 7 (a): Bumps
  # pdf('figuras/bumps.pdf', height=3.8, width=6.75)
  # par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  # plot(x, alpha[,1], type='l', ylab='y')
  # dev.off()
  
  ## Figura 7 (b): Doppler
  # pdf('figuras/doppler.pdf', height=3.8, width=6.75)
  # par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  # plot(x, alpha[,2], type='l', ylab='y')
  # dev.off()


# Figura 7: Amostra agregada
  pdf('figuras/ExemploAmostraFuncional.pdf', height=3, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1, las=1)
  plot(1, xlab='x', ylab='y', ylim=range(A), xlim=c(0,1), type='n')
  for (i in 1:I) lines(x, A[,i], col=i)
  dev.off()
  
# Figura 8: Funções recuperadas
  A_hat <- desagrega(A, y, policy='sure')
  ## Figura 8 (a): bumps recuperada
  pdf('figuras/ExemploFunRecup1.pdf', height=4, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  plot(x, A_hat[,1], type='l', ylab='y', xlab='x', lwd=2,
       ylim=range(c(A_hat[,1], alpha[,1])))
  lines(x, alpha[,1], col='blue')
  dev.off()
  
  ## Figura 8 (b): doppler recuperada
  pdf('figuras/ExemploFunRecup2.pdf', height=4, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.15, las=1)
  plot(x, A_hat[,2], type='l', ylab='y', xlab='x', lwd=2,
       ylim=range(c(A_hat[,2], alpha[,2])))
  lines(x, alpha[,2], col='blue')
  dev.off()
}



# Figura 9: Erro positivo
{
  # definindo parâmetros
  set.seed(282829)
  M <- 64; x <- seq(1, M)/M
  beta <- 144/49; lambda <- 48/49  # erro
  tau <- 5; alpha <- 0.9           # priori
  
  # gerando amostra
  f <- DJ.EX(M)$heavi
  e <- matrix(rgamma(M, shape=beta, rate=lambda), nrow=M)
  y <- f + e

## Figura 9 (a): Amostra
  pdf('figuras/ExemploErroPositivoAmostra.pdf', height=3.8, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.1, las=1)
  plot(x, y, type='b', ylim=c(-14, 15), pch=20)
  lines(x, f, type='b', col='blue', pch=20)
  legend('topright', bty='n', legend=c('Amostra', 'Heavisine'),
         col=c('black', 'blue'), lwd=2)
  dev.off()
  
  # dwt
  d <- wd(y, filter.number=5, family='DaubExPhase')
  p <- gera_ponto(d=d, filter.number=5, family='DaubExPhase')
  # post_gamma(p, d, beta, lambda, tau=5, alpha=0.8)
  delta_d <- ram_gamma_s(p, S_1=NULL, d, n_ite=100000, alpha=alpha, tau=tau,
                         beta=beta, lambda=lambda, gamma=2/3)
  delta_d_bt <- colMeans(delta_d$theta[seq(1000, 100000, 5),])
  
  # IDWT
  f_hat <- GenW(M, filter.number=5, family='DaubExPhase') %*% delta_d_bt
  

## Figura 9 (b): Função estimada
  pdf('figuras/ExemploErroPositivoEstimativa.pdf', height=3.8, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.1, las=1)
  plot(x, f_hat, ylab='y', type='b', pch=20)
  lines(x, f, type='b', col='blue', pch=20)
  legend('topright', bty='n', legend=c('Função estimada', 'Heavisine'),
         col=c('black', 'blue'), lwd=2)
  dev.off()
}


# Figura 10 e 11
{
  # Tecator
  fun_tec <- tecator$absorp.fdata[,21:84]
  x <- tecator$absorp.fdata$argvals[21:84]
  pesos <- t(tecator$y)

## Figura 10: Amostra de curvas agregadas do Tecator
  pdf('figuras/AplicacaoTecatorAmostra.pdf', height=3.1, width=6.75)
  par(mar=c(4.5, 4.7, 1, 1))
  plot(fun_tec, ylab='Absorbância', xlab='Comprimento de onda (nm)', main='')
  dev.off()
  
  # Bayesiano
  D <- apply(fun_tec$data, MARGIN=1, \(x) wd(x, filter.number=5, family='DaubExPhase'))
  D_shrink_bayes <- sapply(D, \(d)
            c(accessC(d, lev=0),
              logis_shrink(d$D, 0.9, s=mad(accessD(d, lev=nlevelsWT(d)-1)), tau=3)))
  Gamma_hat_bayes <-  D_shrink_bayes %*% t(pesos) %*% solve(pesos %*% t(pesos))
  alpha_hat_bayes <- GenW(n=64, filter.number=5, family='DaubExPhase') %*% Gamma_hat_bayes
  
  # Frequentista
  # D_shrink <- sapply(D, \(d) {
  #   d_hat <- threshold(d)
  #   c(accessC(d_hat, level=0), d_hat$D)
  # })
  # Gamma_hat <- D_shrink %*% t(pesos) %*% solve(pesos %*% t(pesos))
  # alpha_hat <- GenW(n=64, filter.number=5, family='DaubExPhase') %*% Gamma_hat
  
## Figura 11: Funções constituintes
  # Figura 11 (a): componente de gordura
  pdf('figuras/AplicacaoGorduraEst.pdf', height=4.7, width=6.5)
  par(mar=c(5.3, 6.2, 1, 1), mgp=c(4,.75,0), cex.lab=1.9, cex.axis=1.3, las=1)
  plot(x, alpha_hat_bayes[,1], type='l', lwd=2,
       ylab=expression(hat(alpha)[1](t)), xlab='Comprimento de onda (nm)')
  dev.off()
  
  # Figura 11 (b): componente de água
  pdf('figuras/AplicacaoAguaEst.pdf', height=4.7, width=6.5)
  par(mar=c(5.3, 6.2, 1, 1), mgp=c(4,.75,0), cex.lab=1.9, cex.axis=1.3, las=1)
  plot(x, alpha_hat_bayes[,2], type='l', lwd=2,
       ylab=expression(hat(alpha)[2](t)), xlab='Comprimento de onda (nm)')
  dev.off()
  
  # Figura 11 (c): componente de proteína
  pdf('figuras/AplicacaoProteinaEst.pdf', height=4.7, width=6.5)
  par(mar=c(5.3, 6.2, 1, 1), mgp=c(4,.75,0), cex.lab=1.9, cex.axis=1.3, las=1)
  plot(x, alpha_hat_bayes[,3], type='l', lwd=2,
       ylab=expression(hat(alpha)[3](t)), xlab='Comprimento de onda (nm)')
  dev.off()
}
