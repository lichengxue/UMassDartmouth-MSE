Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["MAA"]]
Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["NAA"]]

Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["waa_catch"]]


# Possibly all combined
test <- Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["waa_ssb"]]



## Make a for loop to go through mods and seeds

Mods_proj[[1]][["Mod1"]][["om"]][["years"]] # MAKE OUTSIDE LOOPS HERE

years <- Mods_proj[[1]][["Mod1"]][["om"]][["years"]] #CHANGE TO ITERATE MODELS

ages <- c()
ssb <- c()

for(k in 1:6){ # 6 ages
  ages <- append(ages,rep(k, length(years)))
  ssb <- append(ssb,Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["waa_ssb"]][,,k])
  
}

cohorts <- data.frame(
  EM = rep('TESTNAME', length(years)),
  Seed = rep(1, length(years)), # CHANGE TO SEED INDEX
  Year = years,
  Age = ages,
  WAA_SSB = ssb
  
)

cohorts <- cohorts %>%
  rowwise() %>%
  mutate(Cohort = (Year + 1) - Age)

cohorts$Cohort <- as.factor(cohorts$Cohort)


ggplot(cohorts, aes(fill = Cohort, x = Year, y = WAA_SSB)) +
  geom_bar(position = 'stack', stat = 'identity')








###################################

Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["MAA"]]
Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["NAA"]]

Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["waa_catch"]]


# Possibly all combined
test <- Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["MAA"]] * Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["NAA"]]#Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["waa_ssb"]]



## Make a for loop to go through mods and seeds

Mods_proj[[1]][["Mod1"]][["om"]][["years"]] # MAKE OUTSIDE LOOPS HERE

years <- Mods_proj[[1]][["Mod1"]][["om"]][["years"]] #CHANGE TO ITERATE MODELS

ages <- c()
ssb <- c()

for(k in 1:6){ # 6 ages
  ages <- append(ages,rep(k, length(years)))
  ssb <- append(ssb, test[1,1,,k])#append(ssb,Mods_proj[[1]][["Mod1"]][["om"]][["rep"]][["waa_ssb"]][,,k])
  
}

cohorts <- data.frame(
  EM = rep('TESTNAME', length(years)),
  Seed = rep(1, length(years)), # CHANGE TO SEED INDEX
  Year = years,
  Age = ages,
  WAA_SSB = ssb
  
)

cohorts <- cohorts %>%
  rowwise() %>%
  mutate(Cohort = (Year + 1) - Age) %>%
  filter(Year >= 2023)

cohorts$Cohort <- as.factor(cohorts$Cohort)


ggplot(cohorts, aes(fill = Cohort, x = Year, y = WAA_SSB)) +
  geom_bar(position = 'stack', stat = 'identity')
