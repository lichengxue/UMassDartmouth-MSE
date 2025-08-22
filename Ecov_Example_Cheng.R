library(wham)
library(whamMSE)
library(dplyr)

#' A flexible function to set up a WHAM OM input object
#'
#' @param base_om The initial fitted WHAM model object (your 'om_ecov').
#' @param hist_dat The historical data input (your 'gb_dat').
#' @param ecov_dat The specific environmental covariate data list for the scenario.
#' @return A fully configured WHAM input list for one scenario.

setup_om_input <- function(base_om, hist_dat, ecov_dat) {
  if (!requireNamespace("wham", quietly = TRUE)) {
    stop("The 'wham' package is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("whamMSE", quietly = TRUE)) {
    stop("The 'whamMSE' package is required. Please install it.", call. = FALSE)
  }
  # -- Define Model Dimensions and Timeframe --
  n_stocks <- 1
  n_regions <- 1
  n_ages <- 6
  year_start <- 1973
  year_end   <- 2022
  MSE_years  <- 30
  n_proj_years <- year_end - year_start + 1 + MSE_years
  
  # -- Extend Biological Data --
  user_maturity <- array(NA, dim = c(n_stocks, n_proj_years, n_ages))
  user_maturity[, 1:50, ] <- base_om$input$data$mature
  for (i in 51:n_proj_years) user_maturity[, i, ] <- base_om$input$data$mature[, 50, , drop = FALSE]
  
  user_waa <- list()
  user_waa$waa <- array(NA, dim = c(5, n_proj_years, n_ages))
  user_waa$waa[, 1:50, ] <- base_om$input$data$waa
  for (i in 51:n_proj_years) user_waa$waa[, i, ] <- base_om$input$data$waa[, 50, , drop = FALSE]
  user_waa$waa_pointer_fleets <- base_om$input$data$waa_pointer_fleets
  user_waa$waa_pointer_indices <- base_om$input$data$waa_pointer_indices
  user_waa$waa_pointer_totcatch <- base_om$input$data$waa_pointer_ssb
  user_waa$waa_pointer_ssb <- base_om$input$data$waa_pointer_ssb
  user_waa$waa_pointer_M <- base_om$input$data$waa_pointer_M
  
  fracyr_spawn <- hist_dat[[1]]$dat$fracyr_spawn
  
  # -- Define Catch and Survey Index Specs --
  catch_info <- list(catch_cv = 0.05, catch_Neff = 50, use_agg_catch = 1, use_catch_paa = 1)
  index_info <- list(
    index_cv = rep(0.5, 3), index_Neff = rep(25, 3), fracyr_indices = c(0.25, 0.75, 0.25),
    q = c(2.110e-4, 2.247e-4, 2.683e-4), use_indices = rep(1, 3), use_index_paa = rep(1, 3),
    units_indices = rep(2, 3), units_index_paa = rep(2, 3)
  )
  
  # -- Generate Basic Model Structure --
  info <- whamMSE::generate_basic_info(
    n_stocks = 1L, n_regions = 1L, n_indices = 3L, n_fleets = 1L, n_seasons = 1L,
    base.years = year_start:year_end, n_feedback_years = MSE_years, n_ages = 6,
    catch_info = catch_info, index_info = index_info, user_waa = user_waa,
    user_maturity = user_maturity, fracyr_spawn = fracyr_spawn
  )
  
  basic_info <- info$basic_info
  catch_info_use <- info$catch_info
  index_info_use <- info$index_info
  F_info <- info$F
  F_info$F[1:50, ] <- base_om$rep$Fbar[1:50, 1, drop = FALSE]
  
  # -- Define Model Processes --
  sel3 <- list(
    model = c("age-specific", "logistic", "logistic", "logistic"), re = c("ar1_y", "none", "none", "none"),
    initial_pars = list(c(0.017,0.249,0.751,1,1,1), c(2.305,0.327), c(1.611,0.483), c(2.135,0.217)),
    fix_pars = list(c(4:6), NULL, NULL, NULL)
  )
  M <- list(model = "constant", initial_means = array(c(0.57, 0.33, 0.26, 0.23, 0.22, 0.22), dim = c(n_stocks, n_regions, n_ages)))
  sigma_vals <- array(0.5, dim = c(n_stocks, n_regions, n_ages)); sigma_vals[,,1] <- 0.655; sigma_vals[,,2:6] <- 0.720
  NAA_re <- list(
    N1_model = rep("age-specific-fe", n_stocks), sigma = rep("rec+1", n_stocks), cor = rep("iid", n_stocks),
    recruit_model = 3, recruit_pars = list(c(exp(8.082), 1.018e-4)), sigma_vals = sigma_vals
  )
  
  # -- Prepare the WHAM input object, now using the specific ecov_dat --
  om_input <- wham::prepare_wham_input(
    basic_info = basic_info, selectivity = sel3, M = M, NAA_re = NAA_re, ecov = ecov_dat,
    catch_info = catch_info_use, index_info = index_info_use, F = F_info,
    age_comp = "logistic-normal-pool0"
  )
  
  # -- Perform ALL subsequent modifications on this specific object --
  waa_info <- info$par_inputs$user_waa
  om_input <- whamMSE::update_waa(om_input, waa_info = waa_info)
  
  om_input$par$Ecov_process_pars <- base_om$parList$Ecov_process_pars
  om_input$par$Ecov_beta_R <- base_om$parList$Ecov_beta_R
  om_input$par$catch_paa_pars <- base_om$parList$catch_paa_pars
  om_input$par$index_paa_pars <- base_om$parList$index_paa_pars
  om_input$par$sel_repars <- base_om$parList$sel_repars
  for (i in 1:n_regions) om_input$par$log_N1[i, i, ] <- base_om$parList$log_NAA[i, i, 1, ]
  
  input_hist <- wham::prepare_wham_input(hist_dat)
  om_input$data$agg_index_sigma[1:50, ] <- input_hist$data$agg_index_sigma
  om_input$data$use_indices[1:50, ] <- input_hist$data$use_indices
  om_input$data$use_index_paa[1:50, ] <- input_hist$data$use_index_paa
  
  generate_remove_years <- function(mat) {
    lz <- lapply(seq_len(ncol(mat)), function(j) which(mat[, j] == 0)); max_rows <- max(sapply(lz, length))
    if (max_rows == -Inf) return(NULL)
    out <- matrix(NA, nrow = max_rows, ncol = length(lz))
    for (j in seq_along(lz)) out[1:length(lz[[j]]), j] <- lz[[j]]
    out
  }
  remove_agg_years1 <- generate_remove_years(input_hist$data$use_indices)
  remove_paa_years1 <- generate_remove_years(input_hist$data$use_index_paa)
  
  for (i in 51:n_proj_years) om_input$data$agg_index_sigma[i, ] <- input_hist$data$agg_index_sigma[50, ]
  
  om_input <- whamMSE::update_input_index_info(
    om_input, agg_index_sigma = om_input$data$agg_index_sigma, index_Neff = om_input$data$index_Neff,
    remove_agg = TRUE, remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
    remove_paa = TRUE, remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
  )
  
  return(om_input)
}

# -- 1. Initial Data Loading and Setup --
gb_dat <- read_asap3_dat("ASAPfiles_5.14Pull/GBK.DAT")
om_ecov <- readRDS("om_ecov.rds") # This is your base operating model
env.dat_me <- read.csv("CI_indices.csv") %>% filter(Year > 1971)

# -- 2. Create the Two Environmental Scenarios --
# Base environmental data structure
ecov_om <- list(
  label = "bt_temp",
  mean = as.matrix(c(env.dat_me$bt_temp, rep(NA, 31))), # Future values will be filled
  logsigma = "est_1",
  year = 1972:2052,
  use_obs = matrix(1, ncol = 1, nrow = 81),
  process_model = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear")
)

# Scenario 1: Constant Future (mean of historical)
ecov_om_C <- ecov_om
ecov_om_C$mean[51:81,] <- mean(env.dat_me$bt_temp)

# Scenario 2: Projected Future Trend
om_trend <- summary(lm(bt_temp ~ Year, data = env.dat_me))
projection_om <- data.frame(year_proj = c(2022:2052))
error_sd <- om_trend[["sigma"]] # Correct way to get residual standard error
projection_om$SD <- rnorm(n = nrow(projection_om), mean = 0, sd = error_sd)
projection_om <- projection_om %>%
  rowwise() %>%
  mutate(Temp_proj = om_trend[["coefficients"]][1,1] + (om_trend[["coefficients"]][2,1] * year_proj) + SD)

ecov_om_P <- ecov_om
ecov_om_P$mean[51:81,] <- projection_om$Temp_proj

# -- 3. Build the Full Input Object for Each Scenario --
# (Assuming the setup_om_input function from above is in your environment)

# Build the input object for the CONSTANT scenario
om_input_C <- setup_om_input(base_om = om_ecov, hist_dat = gb_dat, ecov_dat = ecov_om_C)

# Build the input object for the PROJECTED scenario
om_input_P <- setup_om_input(base_om = om_ecov, hist_dat = gb_dat, ecov_dat = ecov_om_P)

# -- 4. Prepare and Run the Models --
# You now have two distinct, complete input objects: om_input_C and om_input_P

# Prepare the Constant OM
om_C <- fit_wham(om_input_C, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

# Prepare the Projected OM
om_P <- fit_wham(om_input_P, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

om_C$parList$Ecov_process_pars
om_P$parList$Ecov_process_pars


plot(om_C$input$data$Ecov_obs,type = "l", col="red") 
lines(om_C$rep$Ecov_x,type = "l", col="blue") # Not surprising becuase random effects on Ecov was set ON in om_input_C 

plot(om_P$input$data$Ecov_obs,type = "l", col="red") 
lines(om_P$rep$Ecov_x,type = "l", col="blue") # Not surprising becuase random effects on Ecov was set ON in om_input_P 

# input <- unfitted.om$input
Ecov_re <- readRDS("Ecov_re.RDS") # this MUST come from a fitted OM (estimated random effects of Ecov from a fitted.mod) 

#' Finalize a WHAM OM input for simulation with fixed environmental conditions
#'
#' @param input A fully prepared WHAM input object (e.g., om_input_P).
#' @param scenario_type A character string, either 'P' for projected or 'C' for constant.
#' @param base_om The original fitted model from which to get historical Ecov effects.
#' @param ecov_lmean The long-term mean of the environmental covariate.
#' @return A final, unfitted OM object ready for simulation.

finalize_om_for_simulation <- function(input, scenario_type, Ecov_re, ecov_lmean) {
  
  # Extract historical environmental random effects from the base model
  Ecov_re_hist <- Ecov_re
  input$par$Ecov_re[1:51,] <- Ecov_re_hist
  
  # Set future environmental random effects based on the scenario
  if (scenario_type == "P") {
    # For the projected trend, the future "random effect" is the difference
    # between the projected value and the long-term process mean.
    input$par$Ecov_re[52:81,] <- input$data$Ecov_obs[52:81,] - ecov_lmean
  } else { # Assumes "C" for constant scenario
    # For a constant future, the true value should equal the mean.
    # Since Ecov_true = Ecov_mean + Ecov_re, the future random effect must be 0.
    input$par$Ecov_re[52:81,] <- 0
  }
  
  # Tell WHAM not to simulate these values, but to use the ones we just set
  input$data$do_simulate_Ecov_re <- 0
  
  # Remove Ecov_re from the list of random effects TMB will estimate
  if ("Ecov_re" %in% input$random) {
    input$random <- input$random[input$random != "Ecov_re"]
  }
  
  # Build the final TMB object without fitting
  unfitted_om <- fit_wham(input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
  
  return(unfitted_om)
}

# Define the long-term mean of your environmental covariate.
# Your code used 7.341484. Let's define it explicitly.
long_term_bt_mean <- 7.341484 

# Finalize the CONSTANT scenario OM
unfitted_om_C <- finalize_om_for_simulation(
  input = om_input_C, 
  scenario_type = "C", 
  Ecov_re = Ecov_re, 
  ecov_lmean = long_term_bt_mean
)

# Finalize the PROJECTED scenario OM
unfitted_om_P <- finalize_om_for_simulation(
  input = om_input_P, 
  scenario_type = "P", 
  Ecov_re = Ecov_re, 
  ecov_lmean = long_term_bt_mean
)

# You now have two distinct OMs ready for your MSE loops:
# unfitted_om_C and unfitted_om_P

hcr <- list(hcr.type = 1, hcr.opts = list(use_FXSPR = TRUE, percentFXSPR = 75))

# Define MSE schedule
assess.interval <- 3
base.years <- year_start:year_end
assess.years <- seq(tail(base.years,1), tail(om_ecov$years,1)-assess.interval, by = assess.interval)

# Constant 
i = 1
random1 = random[-2] # We want to turn off the random effects on Ecov
om_with_data <- update_om_fn(unfitted_om_C, seed = 123+i, random = random1)

plot(om_with_data$input$data$Ecov_obs, type = "l")
lines(om_with_data$rep$Ecov_x, type = "l", col = "red") # Not surprising with observation error

# Projected Temp 
i = 1
random1 = random[-2] # We want to turn off the random effects on Ecov
om_with_data <- update_om_fn(unfitted_om_P, seed = 123+i, random = random1)

plot(om_with_data$input$data$Ecov_obs, type = "l")
lines(om_with_data$rep$Ecov_x, type = "l", col = "red") # Not surprising with observation error

# No Ecov Assumed in the EM
mod1 <- loop_through_fn(
  om = om_with_data,
  em_info = info,
  random = random1,
  M_em = M,
  sel_em = sel3,
  NAA_re_em = NAA_re,
  age_comp_em = "logistic-normal-pool0",
  em.opt = list(separate.em = FALSE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
  update_index_info = list(
    agg_index_sigma = input_Ecov$data$agg_index_sigma,
    index_Neff = input_Ecov$data$index_Neff,
    remove_agg = TRUE, remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
    remove_paa = TRUE, remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
  ),
  update_catch_info = list(
    agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
    catch_Neff = input_Ecov$data$catch_Neff
  ),
  assess_years = assess.years,
  assess_interval = assess.interval,
  base_years = base.years,
  year.use = 50,
  add.years = TRUE,
  seed = 123 + i,
  hcr = hcr,
  save.sdrep = FALSE,
  save.last.em = TRUE,
  FXSPR_init = 0.5
)


plot(mod[["om"]][["rep"]][["Ecov_x"]], type = "l")
lines(mod[["om"]][["input"]][["data"]][["Ecov_obs"]], type = "l", col = "red")
plot_wham_output(mod$em_full[[1]],out.type = "html")
mod$om$input$par$Ecov_process_pars
mod$em_full[[1]]$parList$Ecov_process_pars

# Ecov Assumed in the EM! (See below as an example for Ecov_C/Ecov_P)
ecov_me <- list(
  label = "bt_temp",
  mean = as.matrix(c(env.dat_me$bt_temp,rep(mean(env.dat_me$bt_temp),31))),
  # logsigma = as.matrix(rep(log(0.4733709), 50)),
  logsigma = "est_1",
  year = 1972:2052,
  use_obs = matrix(1, ncol = 1, nrow = 81),
  process_model = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear")
)

mod2 <- loop_through_fn(
  om = om_with_data,
  em_info = info,
  random = random1,
  ecov_em = ecov_me,
  M_em = M,
  sel_em = sel3,
  NAA_re_em = NAA_re,
  age_comp_em = "logistic-normal-pool0",
  em.opt = list(separate.em = FALSE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
  update_index_info = list(
    agg_index_sigma = input_Ecov$data$agg_index_sigma,
    index_Neff = input_Ecov$data$index_Neff,
    remove_agg = TRUE, remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
    remove_paa = TRUE, remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
  ),
  update_catch_info = list(
    agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
    catch_Neff = input_Ecov$data$catch_Neff
  ),
  assess_years = assess.years,
  assess_interval = assess.interval,
  base_years = base.years,
  year.use = 50,
  add.years = TRUE,
  seed = 123 + i,
  hcr = hcr,
  save.sdrep = FALSE,
  save.last.em = TRUE,
  FXSPR_init = 0.5
)

# plot_wham_output(mod2$em_full[[1]],out.type = "html")

plot(mod2$em_full[[1]]$input$data$Ecov_obs, type = "l")
lines(mod2$em_full[[1]]$rep$Ecov_x, type = "l", col = "red")
lines(mod2$om$rep$Ecov_x, type = "l", col = "blue")

data <- update_om_fn(om_ecov, seed = 1234, random = random1)
plot(data$rep$Ecov_x)
plot(data$input$data$Ecov_obs,col = "red", type = "l")

