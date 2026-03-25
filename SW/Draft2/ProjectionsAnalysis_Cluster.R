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


mods <- lapply(nsim, function(r) {
  tryCatch({
    mod_list <- lapply(model_nums, function(m) {
      simLocations <- paste0('Outputs/block', m)
      newr <- r + ((m -1) * 100)
      file_path <- file.path(sprintf(paste0('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/Projections/',simLocations,"/block_%d_sim_%d_output.rds"), m, newr))
      readRDS(file_path)
      #print(file_path)
    })

    names(mod_list) <- paste0("Mod", model_nums)

    return(mod_list)
  }, error=function(e){})
})


outfile <- file.path(outdir, "MSEmods_50yr.rds")
saveRDS(mods, file = outfile)

#saveRDS(mods, file = paste0(modelLocations,"/MSEmods_50yr.rds")) # CHANGE HERE filename

# Throw away whole seed if one of the models did not converge

# for (i in 1:length(mods)){
#   for(j in 1:length(mods[[i]])){
#     if(is.null(mods[[i]][[j]])){
#       mods[[i]] <- NA
#     }
#   }
#
# }

# mods <- readRDS('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/Projections/Outputs/MSEmods_50yr.rds')

# Create Plot Output
# test <- mods
mods[sapply(mods, is.null)] <- NULL

plot_mse_output(mods,
                main_dir = modelLocations,
                #output_dir = "Plots", # TO DO: MAKE SUBFOLDER IN DIRECTORY
                output_format = c("png"), # or html but this is giving me an error
                width = 10, height = 7, dpi = 300,
                col.opt = "D",
                new_model_names = ModNames,
                base.model = "Mod_AR1Ecov",
                start.years = 51,
                use.n.years.first = 1,
                use.n.years.last = 10)


# Skipping spider and boxplot fixes FOR NOW

#


Proj_ests <- data.frame(Year = c(),
                        EM_name = c(),
                        seed_no = c(),
                        Temp_OM = c(),
                        Temp_OM_Obs = c(),
                        Temp_EM = c(),
                        Temp_EM_Obs = c(),
                        SSB_OM = c(), # SSB_real
                        SSB_EM = c(),
                        FBar = c(), #f bar OM
                        Pred_catch = c(), #Predited catch om
                        Catch_advice = c(),  # catch Advice from EM
                        #THESE STILL AREN'T RUNNING. WHY
                        # log_FMSY_stat = c(), #IN LOG SPACE. TRANSFORM
                        # log_SSB_FXSPR_stat = c(), #IN LOG SPACE. TRANSFORM
                        # log_FBar_SXPR_stat = c() #IN LOG SPACE. TRANSFORM
                        recPars_a = c(), #IN LOG SPACE. TRANSFORM ? # EM Stock Recruit Relationship 1 # Intercept
                        recPars_b = c(), #IN LOG SPACE. TRANSFORM ? # EM Stock Recruit Relationship 1
                        Ecov_beta_R_OM= c(), #logitrho for ecov param
                        Ecov_beta_R_EM = c(),
                        OM_RecPar_A = c(),#IN LOG SPACE. TRANSFORM # Real stock recruit relationship
                        OM_RecPar_B = c(),#IN LOG SPACE. TRANSFORM # Real stock recruit relationship
                        EM_q = c(),
                        OM_q = c(),
                        OM_NAA_oldest = c(),
                        OM_NAA_youngest = c(),
                        EM_NAA_oldest = c(),
                        EM_NAA_youngest = c(),

                        OM_NAAsigma_oldest = c(),
                        OM_NAAsigma_youngest = c(),
                        EM_NAAsigma_oldest = c(),
                        EM_NAAsigma_youngest = c(),
                        OM_NAA_devs = c(),
                        EM_NAA_devs = c()
)


pad_na <- function(x, target_len) {
  length(x) <- target_len
  return(x)
}

for(i in 1:length(mods)){ # i = seed, j = model
  if(!is.null(mods[[i]])){
  for(k in 1:length(ModNames)){

    years <- min(mods[[i]][[paste0("Mod",k)]][["em_full"]][[1]][["years"]]):(max(mods[[i]][[paste0("Mod",k)]][["em_full"]][[1]][["years"]]) + 3)
    em_name <- rep(ModNames[k], length(years))
    Seed <- rep(i, length(years))
    om_ecovx <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["Ecov_x"]][1:101,1]
    om_ecovObs <- mods[[i]][[paste0("Mod",k)]][["om"]][["input"]][["data"]][["Ecov_obs"]][1:101,1]
    em_ecovx <- mods[[i]][[paste0("Mod",k)]][["em_full"]][[1]][["rep"]][["Ecov_x"]][1:101,1]
    em_ecovObs <- mods[[i]][[paste0("Mod",k)]][["em_input"]][[1]][["data"]][["Ecov_obs"]][1:101,1]

    SSB_om <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["SSB"]]
    SSB_em <- mods[[i]][[paste0("Mod",k)]][["em_full"]][[1]][["rep"]][["SSB"]]
    FBar_om <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["Fbar"]][,1]
    pred_catch_om <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["pred_catch"]]
    cat_advice_em <- c(rep(NA, 51), unlist(mods[[i]][[paste0("Mod",k)]][["catch_advice"]]))

    Ecov_beta_R_OM <- rep(mods[[i]][[paste0("Mod",k)]][["om"]][["par"]][["Ecov_beta_R"]], length(years))

    om_recpar_a <- rep(mods[[i]][[paste0("Mod",k)]][["om"]][['parList']][['mean_rec_pars']][,1], length(years))
    om_recpar_b <- rep(mods[[i]][[paste0("Mod",k)]][["om"]][['parList']][['mean_rec_pars']][,2], length(years))

    em_q <- rep(mods[[i]][[paste0("Mod",k)]][["em_full"]][[1]][["parList"]][["logit_q"]][1], length(years))
    om_q <- rep(mods[[i]][[paste0("Mod",k)]][["om"]][['parList']][['logit_q']][1], length(years))

    om_naa_oldest <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["NAA"]][,,,6]
    om_naasigma_oldest <- rep(mods[[i]][[paste0("Mod",k)]][["om"]][["parList"]][["log_NAA_sigma"]][,,6], length(years))
    om_naa_youngest <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["NAA"]][,,,1]
    om_naasigma_youngest <- rep(mods[[i]][[paste0("Mod",k)]][["om"]][["parList"]][["log_NAA_sigma"]][,,1], length(years))

    em_naa_oldest <- mods[[i]][[paste0("Mod",k)]][["em_list"]][[length(mods[[i]][[paste0("Mod", k)]][["em_list"]])]][["NAA"]][,,,6]
    em_naa_youngest <- mods[[i]][[paste0("Mod",k)]][["em_list"]][[length(mods[[i]][[paste0("Mod", k)]][["em_list"]])]][["NAA"]][,,,1]

    om_naa_devs <- mods[[i]][[paste0("Mod",k)]][["om"]][["rep"]][["NAA_devs"]][1,1,,1]
    em_naa_devs <- mods[[i]][[paste0("Mod",k)]][["em_list"]][[length(mods[[i]][[paste0("Mod", k)]][["em_list"]])]][["NAA_devs"]][1,1,,1]

    # RECORDING INFO THAT IS REPORTED PER ASSESSMENT YEAR
    df_perAssessment <- data.frame(years_trace = years,
                                   recPars_a = NA,
                                   recPars_b = NA,
                                   Ecov_beta_R_EM = NA,
                                   em_naasigma_oldest = NA,
                                   em_naasigma_youngest = NA)
    for(q in 1:length(years)){
      if(years[q] %in% assess.years){
        assess.no <- which(assess.years == years[q])
        df_perAssessment$recPars_a[q:(q+2)] <- mods[[i]][[paste0("Mod",k)]][['par.est']][[assess.no]]$mean_rec_pars[,1]
        df_perAssessment$recPars_b[q:(q+2)] <- mods[[i]][[paste0("Mod",k)]][['par.est']][[assess.no]]$mean_rec_pars[,2]
        df_perAssessment$Ecov_beta_R_EM[q:(q+2)] <- mods[[i]][[paste0("Mod",k)]][['par.est']][[assess.no]][["Ecov_beta_R"]]
        df_perAssessment$em_naasigma_oldest[q:(q+2)] <- mods[[i]][[paste0("Mod",k)]][['par.est']][[assess.no]][["log_NAA_sigma"]][,,6]
        df_perAssessment$em_naasigma_youngest[q:(q+2)] <- mods[[i]][[paste0("Mod",k)]][['par.est']][[assess.no]][["log_NAA_sigma"]][,,1]
      }
    }

    # Compute max length across all vectors going into data_List
    max_len <- max(
      length(years), length(em_name), length(Seed),
      length(om_ecovx), length(om_ecovObs), length(em_ecovx), length(em_ecovObs),
      length(SSB_om), length(SSB_em), length(FBar_om), length(pred_catch_om), length(cat_advice_em),
      length(Ecov_beta_R_OM), length(om_recpar_a), length(om_recpar_b),
      length(em_q), length(om_q),
      length(om_naa_oldest), length(om_naa_youngest), length(em_naa_oldest), length(em_naa_youngest),
      length(om_naasigma_oldest), length(om_naasigma_youngest),
      length(om_naa_devs), length(em_naa_devs),
      nrow(df_perAssessment)
    )

    data_List <- data.frame(
      Year                  = pad_na(years,                            max_len),
      EM_name               = pad_na(em_name,                          max_len),
      seed_no               = pad_na(Seed,                             max_len),
      Temp_OM               = pad_na(om_ecovx,                         max_len),
      Temp_OM_Obs           = pad_na(om_ecovObs,                       max_len),
      Temp_EM               = pad_na(em_ecovx,                         max_len),
      Temp_EM_Obs           = pad_na(em_ecovObs,                       max_len),
      SSB                   = pad_na(SSB_om,                           max_len),
      SSB_EM                = pad_na(SSB_em,                           max_len),
      FBar                  = pad_na(FBar_om,                          max_len),
      Pred_catch            = pad_na(pred_catch_om,                    max_len),
      Catch_advice          = pad_na(cat_advice_em,                    max_len),
      EM_recPars_a          = pad_na(df_perAssessment$recPars_a,       max_len),
      EM_recPars_b          = pad_na(df_perAssessment$recPars_b,       max_len),
      Ecov_beta_R_OM        = pad_na(Ecov_beta_R_OM,                   max_len),
      Ecov_beta_R_EM        = pad_na(df_perAssessment$Ecov_beta_R_EM,  max_len),
      OM_RecPar_A           = pad_na(om_recpar_a,                      max_len),
      OM_RecPar_B           = pad_na(om_recpar_b,                      max_len),
      EM_q                  = pad_na(em_q,                             max_len),
      OM_q                  = pad_na(om_q,                             max_len),
      OM_NAA_oldest         = pad_na(om_naa_oldest,                    max_len),
      OM_NAA_youngest       = pad_na(om_naa_youngest,                  max_len),
      EM_NAA_oldest         = pad_na(em_naa_oldest,                    max_len),
      EM_NAA_youngest       = pad_na(em_naa_youngest,                  max_len),
      OM_NAAsigma_oldest    = pad_na(om_naasigma_oldest,               max_len),
      OM_NAAsigma_youngest  = pad_na(om_naasigma_youngest,             max_len),
      EM_NAAsigma_oldest    = pad_na(df_perAssessment$em_naasigma_oldest,   max_len),
      EM_NAAsigma_youngest  = pad_na(df_perAssessment$em_naasigma_youngest, max_len),
      OM_NAA_devs           = pad_na(om_naa_devs,                      max_len),
      EM_NAA_devs           = pad_na(em_naa_devs,                      max_len)
    )

    Proj_ests <- rbind(Proj_ests, data_List)
  }
}
}

Proj_ests <- Proj_ests[!is.na(Proj_ests$Year), ]


outfile <- file.path(outdir, "Projections_data.rds")
saveRDS(Proj_ests, file = outfile)

# #saveRDS(Proj_ests,paste0(modelLocations, '/Projections_data.rds'))
#
#
# # Calculating and pulling out BRP data
#
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

BRP_data <- BRP_data %>%
  rowwise() %>%
  mutate(MSYOverfishing_OM = OM_Fbar/exp(OM_log_FMSY_static),
         MSYOverfishing_EM = EM_Fbar/exp(EM_log_FMSY_static),

         Overfished_OM = OM_SSB/exp(OM_log_SSB_FXSPR_static),
         Overfished_EM = EM_SSB/exp(EM_log_SSB_FXSPR_static),

         Overfishing_OM = OM_Fbar/exp(OM_log_Fbar_XSPR_static),
         Overfishing_EM = EM_Fbar/exp(EM_log_Fbar_XSPR_static))



#saveRDS(BRP_data,paste0(modelLocations, '/BRP_data.rds'))
outfile <- file.path(outdir, "BRP_data.rds")
saveRDS(BRP_data, file = outfile)










