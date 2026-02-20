library(wham)
library(whamMSE)
library(dplyr)
library(ggplot2)


om_ecov  <- readRDS("C:/Users/swulfing/Documents/GitHub/UMassDartmouth-MSE/SW/om_ecov.rds")

yrs <- list(list(tail(om_ecov$years, 3)))
yrs2 <- list(tail(om_ecov$years, 3))

om_3yr <- project_wham(
  om_ecov,
  proj.opts = list(n.yrs = 50, 
                   #use.last.F = TRUE, 
                   #use.avg.F = FALSE, 
                   use.FXSPR = TRUE,
                   use.FMSY = FALSE, 
                   cont.ecov = TRUE, 
                   use.last.ecov = FALSE, 
                   percentFXSPR = 75,
                   #proj.F = rep(-2.302585, 50), 
                   proj_R_opt = 3, 
                   avg.yrs.R = yrs2, 
                   proj_NAA_opt = 2,
                   avg.yrs.NAA = yrs),
  n.newton = 3,
  do.sdrep = TRUE,
  MakeADFun.silent = FALSE,
  save.sdrep = FALSE,
  check.version = TRUE,
  TMB.bias.correct = FALSE,
  TMB.jointPrecision = FALSE
)

yrs <- list(list(tail(om_ecov$years, 5)))
yrs2 <- list(tail(om_ecov$years, 5))

om_5yr <- project_wham(
  om_ecov,
  proj.opts = list(n.yrs = 50, 
                   #use.last.F = TRUE, 
                   #use.avg.F = FALSE, 
                   use.FXSPR = TRUE,
                   use.FMSY = FALSE, 
                   cont.ecov = TRUE, 
                   use.last.ecov = FALSE, 
                   percentFXSPR = 75,
                   #proj.F = rep(-2.302585, 50), 
                   proj_R_opt = 3, 
                   avg.yrs.R = yrs2, 
                   proj_NAA_opt = 2,
                   avg.yrs.NAA = yrs),
  n.newton = 3,
  do.sdrep = TRUE,
  MakeADFun.silent = FALSE,
  save.sdrep = FALSE,
  check.version = TRUE,
  TMB.bias.correct = FALSE,
  TMB.jointPrecision = FALSE
)


plot_wham_output(
  om_3yr,
  dir.main = file.path('C:/Users/swulfing/Documents/GitHub/UMassDartmouth-MSE/SW/11.04/SensitivityTests/om_3yr')
)

plot_wham_output(
  om_5yr,
  dir.main = file.path('C:/Users/swulfing/Documents/GitHub/UMassDartmouth-MSE/SW/11.04/SensitivityTests/om_5yr')
)


ggplot() +
  # geom_line(aes(x = om_3yr$years_full, y = om_3yr[["rep"]][["SSB"]]), color = 'blue') +
  # geom_line(aes(x = om_5yr$years_full, y = om_5yr[["rep"]][["SSB"]]), color = 'red') +
  geom_line(aes(x = c(om_ecov$years,2023:2072), y = c(om_ecov[["rep"]][["SSB"]], rep(NA, 50))), color = 'black') 


ggplot() +
  geom_line(aes(x = om_3yr$years_full, y = om_3yr[["rep"]][["SSB"]]), color = 'blue') +
  geom_line(aes(x = om_5yr$years_full, y = om_5yr[["rep"]][["SSB"]]), color = 'red') +
  geom_line(aes(x = om_ecov$years, y = om_ecov[["rep"]][["SSB"]]), color = 'black') +
  xlim(2023,2072)

ggplot() +
  geom_line(aes(x = om_3yr$years_full, y = om_3yr[["rep"]][["NAA"]][1,1,,1]), color = 'blue') +
  geom_line(aes(x = om_5yr$years_full, y = om_5yr[["rep"]][["NAA"]][1,1,,1]), color = 'red') +
  geom_line(aes(x = om_ecov$years, y = om_ecov[["rep"]][["NAA"]][1,1,,1]), color = 'black') +
  xlim(2023,2072)




























