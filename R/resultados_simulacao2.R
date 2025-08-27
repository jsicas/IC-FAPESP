# Gráficos e tabelas da simulação 2, cujos cenários são:
#
#   - M = 32 e 64 (pontos por curva);
#   - L = 2 e 4 (quantidade de funções componentes);
#   - SNR = 3 (beta = 36/49 e lambda = 18/49) e 7 (beta = 4 e lambda = 2);
#   - I = 50 (número de curvas);
#   - alpha = 0.75 e tau = 5 (parâmetros da mistura da priori);
#   - n_ite = 50000 (tamanho da cadeia);
#   - rep = 400 (quantidade de replicações).

# packges
require(ggplot2)
require(dplyr)
require(tidyr)
require(xtable)

# carregando simulação
load('R/simulacao2_FAPESP.RData')

# Os objetos seguem a sintaxe: simu_SNR_M_L

L2 <- rbind.data.frame(cbind(snr=3, M=32, L=2, simu_3_32_2),
                       cbind(snr=7, M=32, L=2, simu_7_32_2),
                       cbind(snr=3, M=64, L=2, simu_3_64_2),
                       cbind(snr=7, M=64, L=2, simu_7_64_2)) |>
  mutate(MSE3=NA, MSE4=NA) |> 
  select(snr, M, L, MSE1, MSE2, MSE3, MSE4, AMSE)

L4 <- rbind.data.frame(cbind(snr=3, M=32, L=4, simu_3_32_4),
                       cbind(snr=7, M=32, L=4, simu_7_32_4),
                       cbind(snr=3, M=64, L=4, simu_3_64_4),
                       cbind(snr=7, M=64, L=4, simu_7_64_4))

result <- rbind(L2, L4) |>
  rename(MSE_1=MSE1, MSE_2=MSE2, MSE_3=MSE3, MSE_4=MSE4) |>
  mutate(across(c(snr, M, L), factor))



# Gráficos do relatório FAPESP --------------------------------------------

# Figura 13: Boxplot do AMSE para SNR = 3 e 7 considerando M e L.
{
  pdf('figuras/simu2snr3boxplot.pdf', height=4.3, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.2, las=1)
  result |> filter(snr==3) |>
    boxplot(AMSE ~ M + L, data=_,col=rep(c('blue', 'red'), each=2),
            xlab='M', names=c(32, 64, 32, 64))
    legend('topleft', bty='n', col=c('blue', 'red'), lwd=2, cex=1.3,
           legend=c(expression(L == 2), expression(L == 4)))
  dev.off()
  
  pdf('figuras/simu2snr7boxplot.pdf', height=4.3, width=6.75)
  par(mar=c(4.1, 4.1, 1.5, 1), cex=1.2, las=1)
  result |> filter(snr==7) |>
    boxplot(AMSE ~ M + L, data=_,col=rep(c('blue', 'red'), each=2),
            xlab='M', names=c(32, 64, 32, 64))
    legend('topleft', bty='n', col=c('blue', 'red'), lwd=2, cex=1.3,
           legend=c(expression(L == 2), expression(L == 4)))
  dev.off()
}



# Tabelas -----------------------------------------------------------------

# AMSE
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(AMSE), AMSE=mean(AMSE)) |> 
  mutate(res=paste0(round(AMSE, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -AMSE) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE 1
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_1), MSE_1=mean(MSE_1)) |> 
  mutate(res=paste0(round(MSE_1, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_1) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE 2
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_2), MSE_2=mean(MSE_2)) |> 
  mutate(res=paste0(round(MSE_2, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_2) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE 3
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_3), MSE_3=mean(MSE_3)) |> 
  mutate(res=paste0(round(MSE_3, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_3) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE 4
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_4), MSE_4=mean(MSE_4)) |> 
  mutate(res=paste0(round(MSE_4, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_4) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()
