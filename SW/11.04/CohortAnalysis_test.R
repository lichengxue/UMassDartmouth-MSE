library(tidyr)
library(dplyr)
library(ggplot2)

cohortsN <- data.frame(readRDS('C:/Users/swulfing/Downloads/Cohort_dataNAARE.rds'))
cohortsP <- data.frame(readRDS('C:/Users/swulfing/Downloads/Cohort_dataPROJ.rds'))
cohorts5 <- data.frame(readRDS('C:/Users/swulfing/Downloads/Cohort_data5YR.rds'))

for(i in 1:nrow(cohortsN)){
  cohortsN$BH[i] <- (exp(cohortsN$SR_Alpha[i]) * cohortsN$SSB[i])/(exp(cohortsN$SR_Beta[i]) + cohortsN$SSB[i])
  cohortsN$TempEffect[i] <- cohortsN$Ecov_beta_R_OM[i] * cohortsN$Temp_OM[i]
  cohortsN$Epsilon[i] <- cohortsN$OM_NAA_devs[i]
}

cohorts_NAARE <- cohortsN %>%
  mutate(R_t = log(BH) + TempEffect + Epsilon) %>%
  mutate(Ecov_t = log(BH) + TempEffect) %>%
  group_by(EM_name, seed_no) %>%
  summarize(var_Ecovt = var(Ecov_t),
            var_Rt = var(R_t)) %>%
  mutate(contribution = var_Ecovt/var_Rt)



for(i in 1:nrow(cohortsP)){
  cohortsP$BH[i] <- (exp(cohortsP$SR_Alpha[i]) * cohortsP$SSB[i])/(exp(cohortsP$SR_Beta[i]) + cohortsP$SSB[i])
  cohortsP$TempEffect[i] <- cohortsP$Ecov_beta_R_OM[i] * cohortsP$Temp_OM[i]
  cohortsP$Epsilon[i] <- cohortsP$OM_NAA_devs[i]
}

cohorts_PROJ <- cohortsP %>%
  mutate(R_t = log(BH) + TempEffect + Epsilon) %>%
  mutate(Ecov_t = log(BH) + TempEffect) %>%
  group_by(EM_name, seed_no) %>%
  summarize(var_Ecovt = var(Ecov_t),
            var_Rt = var(R_t)) %>%
  mutate(contribution = var_Ecovt/var_Rt)


for(i in 1:nrow(cohorts5)){
  cohorts5$BH[i] <- (exp(cohorts5$SR_Alpha[i]) * cohorts5$SSB[i])/(exp(cohorts5$SR_Beta[i]) + cohorts5$SSB[i])
  cohorts5$TempEffect[i] <- cohorts5$Ecov_beta_R_OM[i] * cohorts5$Temp_OM[i]
  cohorts5$Epsilon[i] <- cohorts5$OM_NAA_devs[i]
}

cohorts_5YR <- cohorts5 %>%
  mutate(R_t = log(BH) + TempEffect + Epsilon) %>%
  mutate(Ecov_t = log(BH) + TempEffect) %>%
  group_by(EM_name, seed_no) %>%
  summarize(var_Ecovt = var(Ecov_t),
            var_Rt = var(R_t)) %>%
  mutate(contribution = var_Ecovt/var_Rt)
  


cohorts_NAARE$Experiment <- 'NAARE'
cohorts_PROJ$Experiment <- 'PROJs'
cohorts_5YR$Experiment <- '5Yr'

cohort_combine <- rbind(cohorts_NAARE, cohorts_PROJ, cohorts_5YR)


ggplot(cohort_combine, aes(x = EM_name, y = contribution, fill = Experiment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# cohorts1 <- cohorts %>%
#   select(Year, EM_name, seed_no, CAA1, CAA2, CAA3, CAA4, CAA5, CAA6) %>%
#   rowwise() %>%
#   filter(Year >= 2000 & Year <= 2005)
# 
# cohorts1 <- pivot_longer(cohorts1, names_to = 'Age', values_to = 'CAA', cols = starts_with('CAA'))
# 
# cohorts1$Age <- as.numeric(gsub("[^0-9]", "", cohorts1$Age))
# 
# cohorts1 <- cohorts1 %>%
#   rowwise() %>%
#   mutate(Cohort = (Year + 1) - Age) %>%
#   mutate(CAA = CAA / 10^6)
# 
# cohorts1$Cohort <- as.factor(cohorts1$Cohort)
# 
# bartest <- cohorts1 %>%
#   filter(Year == Cohort)
# 
# 
# ggplot(cohorts1, aes(fill = Cohort, x = Year, y = CAA)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   geom_col(
#     data = bartest,
#     aes(x = Year, y = CAA, group = Cohort), 
#     position = "stack", 
#     color = "black",       # Outline color
#     linewidth = 0.2,       # Outline thickness
#     fill = NA              # Keeps the original fill visible
#   ) +
#   scale_y_continuous(labels = scales::comma) + 
#   ylab('Catch At Age (Millions)') 
