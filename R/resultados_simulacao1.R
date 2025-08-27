# Gráficos e tabelas da simulação com os cenários:
#   - M: 32, 512, 2048 (quantidade de pontos por curva);
#   - snr: 1, 3 e 5;
#   - L:
#       2 (bumps e doppler);
#       4 (funções de DJ);
#       5 (funções de DJ e logit).

# packges
require(ggplot2)
require(dplyr)
require(tidyr)
require(xtable)

# carregando simulação
load('R/simulacao1_FAPESP.RData')

# Os objetos seguem a sintaxe: simu_SNR_M_L

L2 <- rbind(cbind(snr=1, M=32,   L=2, simu_1_32_2),
            cbind(snr=3, M=32,   L=2, simu_3_32_2),
            cbind(snr=5, M=32,   L=2, simu_5_32_2),
            cbind(snr=1, M=512,  L=2, simu_1_512_2),
            cbind(snr=3, M=512,  L=2, simu_3_512_2),
            cbind(snr=5, M=512,  L=2, simu_5_512_2),
            cbind(snr=1, M=2048, L=2, simu_1_2048_2),
            cbind(snr=3, M=2048, L=2, simu_3_2048_2),
            cbind(snr=5, M=2048, L=2, simu_5_2048_2)) |> 
  mutate(MSE_3=NA, MSE_4=NA, MSE_5=NA) |> 
  select(snr, M, L, MSE_1, MSE_2, MSE_3, MSE_4, MSE_5, AMSE)

L4 <- rbind(cbind(snr=1, M=32,   L=4, simu_1_32_4),
            cbind(snr=3, M=32,   L=4, simu_3_32_4),
            cbind(snr=5, M=32,   L=4, simu_5_32_4),
            cbind(snr=1, M=512,  L=4, simu_1_512_4),
            cbind(snr=3, M=512,  L=4, simu_3_512_4),
            cbind(snr=5, M=512,  L=4, simu_5_512_4),
            cbind(snr=1, M=2048, L=4, simu_1_2048_4),
            cbind(snr=3, M=2048, L=4, simu_3_2048_4),
            cbind(snr=5, M=2048, L=4, simu_5_2048_4)) |> 
  mutate(MSE_5=NA) |> 
  select(snr, M, L, MSE_1, MSE_2, MSE_3, MSE_4, MSE_5, AMSE)

L5 <- rbind(cbind(snr=1, M=32,   L=5, simu_1_32_5),
            cbind(snr=3, M=32,   L=5, simu_3_32_5),
            cbind(snr=5, M=32,   L=5, simu_5_32_5),
            cbind(snr=1, M=512,  L=5, simu_1_512_5),
            cbind(snr=3, M=512,  L=5, simu_3_512_5),
            cbind(snr=5, M=512,  L=5, simu_5_512_5),
            cbind(snr=1, M=2048, L=5, simu_1_2048_5),
            cbind(snr=3, M=2048, L=5, simu_3_2048_5),
            cbind(snr=5, M=2048, L=5, simu_5_2048_5))

result <- rbind(L2, L4, L5) |> 
  mutate(across(c(snr, M, L), factor))



# Gráficos do relatório FAPESP --------------------------------------------

# Figura 12: AMSE para SNR, L e M
{
  pdf('figuras/AMSEParam.pdf', height=2.4, width=6)
  par(mfrow=c(1,3), mar=c(4.1, 4.1, 1.5, 1), cex=.7, las=1)
  boxplot(AMSE ~ snr, data=result, xlab='SNR')
  boxplot(AMSE ~ L, data=result)
  boxplot(AMSE ~ M, data=result)
  dev.off()
}



# Tabelas -----------------------------------------------------------------
# L2 |>
#   group_by(snr, M) |>
#   reframe(mean(AMSE))

# AMSE
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(AMSE), AMSE=mean(AMSE)) |> 
  mutate(res=paste0(round(AMSE, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -AMSE) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE_1
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_1), MSE_1=mean(MSE_1)) |> 
  mutate(res=paste0(round(MSE_1, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_1) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE_2
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_2), MSE_2=mean(MSE_2)) |> 
  mutate(res=paste0(round(MSE_2, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_2) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE_3
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_3), MSE_3=mean(MSE_3)) |> 
  mutate(res=paste0(round(MSE_3, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_3) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()

# MSE_4
result |>
  group_by(snr, M, L) |> 
  reframe(sd=sd(MSE_4), MSE_4=mean(MSE_4)) |> 
  mutate(res=paste0(round(MSE_4, 3), ' (', round(sd,3), ')')) |> 
  select(-sd, -MSE_4) |> 
  pivot_wider(names_from=L, values_from=res) |> 
  xtable()
