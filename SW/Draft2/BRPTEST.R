args <- commandArgs(trailingOnly = TRUE)

outdir <- args[which(args == "--outdir") + 1]


library(ggridges)
library(wham)
library(whamMSE)
library(dplyr)
library(ggplot2)
library(ggtern)

# Pull and compile simulations across projections

model_nums <- 1:8 # CHANGE HERE Number of Projections
nsim <- c(0:99) # CHANGE HERE number of simulations/seed
ModNames <- c("Mod_AR1Ecov", "Mod_NoEcov", "Mod_HistAvgEcov", "Mod_RecWindEcov", "Mod_RecTrendEcov", "Mod_TermYear", "Mod_BadProj", "Mod_GoodProj") # CHANGE NAME HERE
assess.years <- c(2022, 2025, 2028, 2031, 2034, 2037, 2040, 2043, 2046, 2049, 2052, 2055, 2058, 2061, 2064, 2067, 2070) #PUT THIS BACK WHEN YOU'RE PUTTING IN CLUSTER# CHANGE ASSESSMENT YEARS HERE
MSE_Length <- 3 # CHANGE HERE FOR 5 YR TESTS

modelLocations <- 'Outputs'

mods <- readRDS('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/Projections/Outputs/MSEmods_50yr.rds')
mods[sapply(mods, is.null)] <- NULL

pad_na <- function(x, target_len) {
  length(x) <- target_len
  return(x)
}


BRP_data <- data.frame(Year = c(), # which EM comes up wiht most accurate BRP, too optimistic/too pessimistic eventual over or under exploitation. % model runs em exceeds catch goal
                       EM_name = c(),
                       Seed_no = c(),
                       
                       OM_Fbar = c(),
                       OM_SSB = c(),
                       OM_log_FMSY_static = c(),
                       OM_log_FMSY = c(),
                       OM_log_SSB_FXSPR_static = c(),
                       OM_log_SSB_FXSPR = c(),
                       OM_log_Fbar_XSPR = c(),
                       OM_log_Fbar_XSPR_static = c(),
                       OM_log_FXSPR_static = c(),
                       
                       EM_Fbar = c(),
                       EM_SSB = c(),
                       EM_log_FMSY_static = c(),
                       EM_log_FMSY = c(),
                       EM_log_SSB_FXSPR_static = c(),
                       EM_log_SSB_FXSPR = c(),
                       EM_log_Fbar_XSPR = c(),
                       EM_log_Fbar_XSPR_static = c(),
                       EM_log_FXSPR_static = c()
)



for(r in 1:length(mods)){ # r is seed No, m is model
  if(!is.null(mods[[r]])){
    for(m in 1:length(mods[[r]])){
      if(is.null(mods[[r]][[m]])){
        mods[[r]] <- NA
      }
      
      model <- mods[[r]][[m]]
      input <- model$om$input
      input$random <- NULL
      
      if (!is.null(input$data$FXSPR_init)) {
        input$data$FXSPR_init[] <- 2
      }
      if (!is.null(input$data$FXSPR_static_init)) {
        input$data$FXSPR_static_init[] <- 2
      }
      
      om_brp <- try(
        fit_wham(input,
                 n.newton         = 10,
                 do.fit           = FALSE,
                 do.brps          = TRUE,
                 MakeADFun.silent = TRUE),
        silent = TRUE
      )
      em_size <- length(model[["em_list"]])
      
      year                    <- om_brp[["years"]]
      eM_name                 <- rep(ModNames[m], length(om_brp[["years"]]))
      seed_no                 <- rep(r, length(om_brp[["years"]]))
      
      # OM Info
      oM_Fbar                 <- om_brp[["rep"]][["Fbar"]][,1]
      oM_SSB                  <- om_brp[["rep"]][["SSB"]]
      oM_log_FMSY_static      <- rep(om_brp[["rep"]][["log_FMSY_static"]], length(om_brp[["years"]]))
      oM_log_FMSY             <- om_brp[["rep"]][["log_FMSY"]]
      oM_log_SSB_FXSPR_static <- rep(om_brp[["rep"]][["log_SSB_FXSPR_static"]][1], length(om_brp[["years"]]))
      oM_log_SSB_FXSPR        <- om_brp[["rep"]][["log_SSB_FXSPR"]][,1]
      oM_log_Fbar_XSPR        <- om_brp[["rep"]][["log_Fbar_XSPR"]][,1]
      oM_log_Fbar_XSPR_static <- rep(om_brp[["rep"]][["log_Fbar_XSPR_static"]][1], length(om_brp[["years"]]))
      oM_log_FXSPR_static     <- rep(om_brp[["rep"]][["log_FXSPR_static"]], length(om_brp[["years"]]))
      
      # EM Comparisons
      eM_Fbar                 <- model[["em_list"]][[em_size]][["Fbar"]][,1]
      eM_SSB                  <- model[["em_list"]][[em_size]][["SSB"]]
      eM_log_FMSY             <- model[["em_list"]][[em_size]][["log_FMSY"]]
      eM_log_SSB_FXSPR        <- model[["em_list"]][[em_size]][["log_SSB_FXSPR"]][,1]
      eM_log_Fbar_XSPR        <- model[["em_list"]][[em_size]][["log_Fbar_XSPR"]][,1]
      
      eM_log_FMSY_static      <- rep(NA, 50)
      eM_log_SSB_FXSPR_static <- rep(NA, 50)
      eM_log_FXSPR_static     <- rep(NA, 50)
      eM_log_Fbar_XSPR_static <- rep(NA, 50)
      
      for(j in 1:length(model[["em_list"]])){
        eM_log_FMSY_static      <- append(eM_log_FMSY_static, rep(model[["em_list"]][[j]][["log_FMSY_static"]], MSE_Length))
        eM_log_SSB_FXSPR_static <- append(eM_log_SSB_FXSPR_static, rep(model[["em_list"]][[j]][["log_SSB_FXSPR_static"]][1], MSE_Length))
        eM_log_FXSPR_static     <- append(eM_log_FXSPR_static, rep(model[["em_list"]][[j]][["log_SPR_FXSPR_static"]][1], MSE_Length))
        eM_log_Fbar_XSPR_static <- append(eM_log_Fbar_XSPR_static, rep(model[["em_list"]][[j]][["log_Fbar_XSPR_static"]][1], MSE_Length))
      }
      
      # Compute max length across all vectors going into data_List
      max_len <- max(
        length(year), length(eM_name), length(seed_no),
        length(oM_Fbar), length(oM_SSB), length(oM_log_FMSY_static),
        length(oM_log_FMSY), length(oM_log_SSB_FXSPR_static), length(oM_log_SSB_FXSPR),
        length(oM_log_Fbar_XSPR), length(oM_log_Fbar_XSPR_static), length(oM_log_FXSPR_static),
        length(eM_Fbar), length(eM_SSB), length(eM_log_FMSY_static),
        length(eM_log_FMSY), length(eM_log_SSB_FXSPR_static), length(eM_log_SSB_FXSPR),
        length(eM_log_Fbar_XSPR), length(eM_log_Fbar_XSPR_static), length(eM_log_FXSPR_static)
      )
      
      data_List <- data.frame(
        Year                    = pad_na(year,                    max_len),
        EM_name                 = pad_na(eM_name,                 max_len),
        Seed_no                 = pad_na(seed_no,                 max_len),
        OM_Fbar                 = pad_na(oM_Fbar,                 max_len),
        OM_SSB                  = pad_na(oM_SSB,                  max_len),
        OM_log_FMSY_static      = pad_na(oM_log_FMSY_static,      max_len),
        OM_log_FMSY             = pad_na(oM_log_FMSY,             max_len),
        OM_log_SSB_FXSPR_static = pad_na(oM_log_SSB_FXSPR_static, max_len),
        OM_log_SSB_FXSPR        = pad_na(oM_log_SSB_FXSPR,        max_len),
        OM_log_Fbar_XSPR        = pad_na(oM_log_Fbar_XSPR,        max_len),
        OM_log_Fbar_XSPR_static = pad_na(oM_log_Fbar_XSPR_static, max_len),
        OM_log_FXSPR_static     = pad_na(oM_log_FXSPR_static,     max_len),
        EM_Fbar                 = pad_na(eM_Fbar,                 max_len),
        EM_SSB                  = pad_na(eM_SSB,                  max_len),
        EM_log_FMSY_static      = pad_na(eM_log_FMSY_static,      max_len),
        EM_log_FMSY             = pad_na(eM_log_FMSY,             max_len),
        EM_log_SSB_FXSPR_static = pad_na(eM_log_SSB_FXSPR_static, max_len),
        EM_log_SSB_FXSPR        = pad_na(eM_log_SSB_FXSPR,        max_len),
        EM_log_Fbar_XSPR        = pad_na(eM_log_Fbar_XSPR,        max_len),
        EM_log_Fbar_XSPR_static = pad_na(eM_log_Fbar_XSPR_static, max_len),
        EM_log_FXSPR_static     = pad_na(eM_log_FXSPR_static,     max_len)
      )
      
      BRP_data <- rbind(BRP_data, data_List)
    }
  }
}

# BRP_data <- BRP_data %>%
#   rowwise() %>%
#   mutate(MSYOverfishing_OM = OM_Fbar/exp(OM_log_FMSY_static),
#          MSYOverfishing_EM = EM_Fbar/exp(EM_log_FMSY_static),
#          
#          Overfished_OM = OM_SSB/exp(OM_log_SSB_FXSPR_static),
#          Overfished_EM = EM_SSB/exp(EM_log_SSB_FXSPR_static),
#          
#          Overfishing_OM = OM_Fbar/exp(OM_log_Fbar_XSPR_static),
#          Overfishing_EM = EM_Fbar/exp(EM_log_Fbar_XSPR_static))



#saveRDS(BRP_data,paste0(modelLocations, '/BRP_data.rds'))
outfile <- file.path(outdir, "BRPTEST_data.rds")
saveRDS(BRP_data, file = outfile)








































