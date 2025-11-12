## ============================================================
## PHASE A. Fit a historical OM with Ecov (no projection years)
## ============================================================

# (Optional) install dev versions
# devtools::install_github("lichengxue/whamMSE", ref = "UMassD-MSE", dependencies = FALSE)
# devtools::install_github("timjmiller/wham",     ref = "lab",         dependencies = FALSE)

setwd("~/Desktop/UMassD-Cheng")  # <- ADJUST if needed

suppressPackageStartupMessages({
  library(wham)
  library(whamMSE)
  library(dplyr)
  library(ggplot2)
})

## ------------------------------
## A1) Historical data & inputs
## ------------------------------
gb_dat <- read_asap3_dat("GBK.DAT")
input  <- prepare_wham_input(gb_dat)

# Define historical time span (assessment years)
year_start <- 1973L
year_end   <- 2022L
hist_len   <- year_end - year_start + 1L  # = 50

# MSE_years = 0 here because we are fitting the historical OM only
MSE_years_histfit <- 0L

## Maturity-at-age and Weight-at-age (historical only)
user_maturity <- array(NA_real_, dim = c(1, hist_len, 6))
user_maturity[, 1:hist_len, ] <- input$data$mature  # (fitted ASAP -> WHAM input)

user_waa <- list()
user_waa$waa <- array(NA_real_, dim = c(5, hist_len, 6))
user_waa$waa[, 1:hist_len, ] <- input$data$waa
user_waa$waa_pointer_fleets   <- input$data$waa_pointer_fleets
user_waa$waa_pointer_indices  <- input$data$waa_pointer_indices
user_waa$waa_pointer_totcatch <- input$data$waa_pointer_ssb
user_waa$waa_pointer_ssb      <- input$data$waa_pointer_ssb
user_waa$waa_pointer_M        <- input$data$waa_pointer_M

# Spawning fraction by year (ASAP → WHAM)
fracyr_spawn <- gb_dat[[1]]$dat$fracyr_spawn

## Catch & index observation model (historical)
catch_info <- list(
  catch_cv      = 0.05,
  catch_Neff    = 50,
  use_agg_catch = 1,
  use_catch_paa = 1
)

index_info <- list(
  index_cv        = rep(0.5, 3),
  index_Neff      = rep(25, 3),
  fracyr_indices  = c(0.25, 0.75, 0.25),
  q               = rep(0.2, 3),     # placeholder q for a clean historical fit
  use_indices     = rep(1, 3),
  use_index_paa   = rep(1, 3),
  units_indices   = rep(2, 3),
  units_index_paa = rep(2, 3)
)

## Basic info for historical OM (no feedback/projection years here)
info_hist <- whamMSE::generate_basic_info(
  n_stocks = 1L, n_regions = 1L, n_indices = 3L,
  n_fleets = 1L, n_seasons = 1L,
  base.years = year_start:year_end,
  n_feedback_years = MSE_years_histfit,
  n_ages = 6L,
  catch_info = catch_info, index_info = index_info,
  user_waa = user_waa, user_maturity = user_maturity,
  fracyr_spawn = fracyr_spawn
)

basic_info_hist     <- info_hist$basic_info
catch_info_use_hist <- info_hist$catch_info
index_info_use_hist <- info_hist$index_info
F_info_hist         <- info_hist$F

## Fill observed catch & index series from the ASAP-derived input (historical)
catch_info_use_hist$agg_catch[1:hist_len, ]       <- input$data$agg_catch
catch_info_use_hist$catch_paa[, 1:hist_len, ]     <- input$data$catch_paa
index_info_use_hist$agg_indices[1:hist_len, ]     <- input$data$agg_indices
index_info_use_hist$index_paa[, 1:hist_len, ]     <- input$data$index_paa
index_info_use_hist$use_indices[1:hist_len, ]     <- input$data$use_indices
index_info_use_hist$use_index_paa[1:hist_len, ]   <- input$data$use_index_paa

## Selectivity (1 fleet + 3 indices)
sel_hist <- list(
  model = c("age-specific", "logistic", "logistic", "logistic"),
  re    = c("ar1_y", "none", "none", "none"),
  initial_pars = list(
    c(0.1, 0.25, 0.5, 1, 1, 1),
    c(2, 0.3), c(2, 0.3), c(2, 0.3)
  ),
  fix_pars = list(c(4:6), NULL, NULL, NULL)
)

## Natural mortality, recruitment, and NAA (historical OM)
n_stocks  <- as.integer(basic_info_hist["n_stocks"])
n_regions <- as.integer(basic_info_hist["n_regions"])
n_ages    <- as.integer(basic_info_hist["n_ages"])

M_hist <- list(
  model = "constant",
  initial_means = array(c(0.57, 0.33, 0.26, 0.23, 0.22, 0.22),
                        dim = c(n_stocks, n_regions, n_ages))
)

sigma_vals_hist <- array(0.5, dim = c(n_stocks, n_regions, n_ages))
NAA_re_hist <- list(
  N1_model = rep("equilibrium", n_stocks),  # equilibrium start for the historical fit
  sigma    = rep("rec+1", n_stocks),
  cor      = rep("iid", n_stocks),
  recruit_model = 3,
  sigma_vals    = sigma_vals_hist
)

## Historical Ecov series: CSV provides 1972..2021 (lag lead + hist years)
env.dat_me <- read.csv("CI_indices.csv") %>% filter(Year > 1971)
stopifnot(all(c("Year", "bt_temp") %in% names(env.dat_me)))

ecov_me_hist <- list(
  label           = "bt_temp",
  mean            = as.matrix(env.dat_me$bt_temp),
  logsigma      = as.matrix(rep(log(0.2), length(env.dat_me$bt_temp))),
  year            = env.dat_me$Year,  # expected 1972..2021 (50 pts)
  use_obs         = matrix(1L, ncol = 1, nrow = nrow(env.dat_me)),
  process_model   = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear")
)

## Prepare and fit the historical OM
input_Ecov_hist <- prepare_wham_input(
  basic_info  = basic_info_hist,
  selectivity = sel_hist,
  M           = M_hist,
  NAA_re      = NAA_re_hist,
  ecov        = ecov_me_hist,
  catch_info  = catch_info_use_hist,
  index_info  = index_info_use_hist,
  F           = F_info_hist,
  age_comp    = "logistic-normal-pool0"
)

# Ensure WAA pointers are honored
input_Ecov_hist <- update_waa(input_Ecov_hist, waa_info = info_hist$par_inputs$user_waa)

# Use ASAP-derived index usage/sigmas (already copied via index_info_use_hist)
# but we can explicitly set them again if desired:
input_e <- prepare_wham_input(gb_dat)
input_Ecov_hist$data$agg_index_sigma[1:hist_len, ] <- input_e$data$agg_index_sigma
input_Ecov_hist$data$use_indices[1:hist_len, ]     <- input_e$data$use_indices
input_Ecov_hist$data$use_index_paa[1:hist_len, ]   <- input_e$data$use_index_paa

# Remove years where indices weren’t used
generate_remove_years <- function(mat){
  lz <- lapply(seq_len(ncol(mat)), function(j) which(mat[, j] == 0))
  if (all(lengths(lz) == 0)) return(NULL)
  max_rows <- max(sapply(lz, length))
  out <- matrix(NA_integer_, nrow = max_rows, ncol = length(lz))
  for (j in seq_along(lz)) if (length(lz[[j]])) out[seq_along(lz[[j]]), j] <- lz[[j]]
  out
}
remove_agg_years1 <- generate_remove_years(input_e$data$use_indices)
remove_paa_years1 <- generate_remove_years(input_e$data$use_index_paa)

input_Ecov_hist <- update_input_index_info(
  input_Ecov_hist,
  agg_index_sigma = input_Ecov_hist$data$agg_index_sigma,
  index_Neff      = input_Ecov_hist$data$index_Neff,
  remove_agg = TRUE, remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
  remove_paa = TRUE, remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
)

# Fit historical OM
om_ecov <- fit_wham(input_Ecov_hist, do.fit = TRUE, do.brps = TRUE, 
                    do.retro = FALSE, do.osa = FALSE, MakeADFun.silent = TRUE)
check_convergence(om_ecov)
saveRDS(om_ecov,                "om_ecov.RDS")
saveRDS(om_ecov$parList$Ecov_re, "Ecov_re.RDS")

## ============================================================
## PHASE B. Project Ecov, compile unfitted OM, simulate, and run EMs
## ============================================================

## ------------------------------
## B1) Reload essentials & set horizon
## ------------------------------
gb_dat     <- read_asap3_dat("GBK.DAT")
input_hist <- prepare_wham_input(gb_dat)
om_ecov    <- readRDS("om_ecov.RDS")
Ecov_re    <- readRDS("Ecov_re.RDS")
env.dat_me <- read.csv("CI_indices.csv") %>% filter(Year > 1971)

stopifnot(all(c("Year","bt_temp") %in% names(env.dat_me)))

# Dimensions for OM/EM and MSE horizon
n_stocks   <- 1L
n_regions  <- 1L
n_ages     <- 6L
n_indices  <- 3L
n_fleets   <- 1L
n_seasons  <- 1L

year_start <- 1973L
year_end   <- 2022L
hist_len   <- year_end - year_start + 1L       # = 50

MSE_years       <- 30L                          # feedback/projection years
assess.interval <- 3L                           # assess every 3 years
stopifnot(MSE_years %% assess.interval == 0L)   # enforce divisibility

n_model_years <- hist_len + MSE_years           # = 80 total modeled years
n_proj_years  <- n_model_years                  # arrays span all modeled years

base.years   <- year_start:year_end
assess.years <- seq(tail(base.years, 1), year_end + MSE_years - assess.interval, by = assess.interval)

## ------------------------------
## B2) Biology over historical + projection years
## ------------------------------
user_maturity <- array(NA_real_, dim = c(n_stocks, n_proj_years, n_ages))
user_maturity[, 1:hist_len, ] <- om_ecov$input$data$mature[, 1:hist_len, ]
if (hist_len < n_proj_years) {
  for (iy in (hist_len + 1):n_proj_years) {
    user_maturity[, iy, ] <- om_ecov$input$data$mature[, hist_len, , drop = FALSE]
  }
}

user_waa <- list()
user_waa$waa <- array(NA_real_, dim = c(5, n_proj_years, n_ages))
user_waa$waa[, 1:hist_len, ] <- om_ecov$input$data$waa[, 1:hist_len, ]
if (hist_len < n_proj_years) {
  for (iy in (hist_len + 1):n_proj_years) {
    user_waa$waa[, iy, ] <- om_ecov$input$data$waa[, hist_len, , drop = FALSE]
  }
}
user_waa$waa_pointer_fleets   <- om_ecov$input$data$waa_pointer_fleets
user_waa$waa_pointer_indices  <- om_ecov$input$data$waa_pointer_indices
user_waa$waa_pointer_totcatch <- om_ecov$input$data$waa_pointer_totcatch
user_waa$waa_pointer_ssb      <- om_ecov$input$data$waa_pointer_ssb
user_waa$waa_pointer_M        <- om_ecov$input$data$waa_pointer_M

fracyr_spawn <- gb_dat[[1]]$dat$fracyr_spawn

## Observation model (same as Phase A, but we’ll use different q values typical to your later runs)
catch_info <- list(
  catch_cv       = 0.05,
  catch_Neff     = 50,
  use_agg_catch  = 1,
  use_catch_paa  = 1
)
index_info <- list(
  index_cv        = rep(0.5, n_indices),
  index_Neff      = rep(25,  n_indices),
  fracyr_indices  = c(0.25, 0.75, 0.25),
  q               = c(2.110e-4, 2.247e-4, 2.683e-4),
  use_indices     = rep(1, n_indices),
  use_index_paa   = rep(1, n_indices),
  units_indices   = rep(2, n_indices),
  units_index_paa = rep(2, n_indices)
)

## Basic info for the (projecting) OM/EM build
info <- whamMSE::generate_basic_info(
  n_stocks = n_stocks, n_regions = n_regions, n_indices = n_indices,
  n_fleets = n_fleets, n_seasons = n_seasons,
  base.years = base.years,
  n_feedback_years = MSE_years,
  n_ages = n_ages,
  catch_info = catch_info, index_info = index_info,
  user_waa = user_waa, user_maturity = user_maturity,
  fracyr_spawn = fracyr_spawn
)

basic_info     <- info$basic_info
catch_info_use <- info$catch_info
index_info_use <- info$index_info
F_info         <- info$F

# Initialize historical Fbar (leave projection years to be controlled by HCR)
nF_hist <- min(nrow(F_info$F), hist_len)
F_info$F[1:nF_hist, ] <- om_ecov$rep$Fbar[1:nF_hist, 1, drop = FALSE]

## Selectivity & NAA for projecting OM
sel3 <- list(
  model = c("age-specific", "logistic", "logistic", "logistic"),
  re    = c("ar1_y", "none", "none", "none"),
  initial_pars = list(
    c(0.017, 0.249, 0.751, 1, 1, 1),
    c(2.305, 0.327),
    c(1.611, 0.483),
    c(2.135, 0.217)
  ),
  fix_pars = list(c(4:6), NULL, NULL, NULL)
)

M <- list(
  model = "constant",
  initial_means = array(c(0.57, 0.33, 0.26, 0.23, 0.22, 0.22),
                        dim = c(n_stocks, n_regions, n_ages))
)

sigma_vals <- array(0.5, dim = c(n_stocks, n_regions, n_ages))
sigma_vals[,,1]   <- 0.655
sigma_vals[,,2:6] <- 0.720

NAA_re <- list(
  N1_model      = rep("age-specific-fe", n_stocks),
  sigma         = rep("rec+1", n_stocks),
  cor           = rep("iid", n_stocks),
  recruit_model = 3,
  recruit_pars  = list(c(exp(8.082), 1.018e-4)),
  sigma_vals    = sigma_vals
)

## ------------------------------
## B3) Projected Ecov mean path (1972..2052)
## ------------------------------
# Ecov has a 1-year lag; to model 80 years (1973..2052), we need 81 Ecov rows: 1972..2052
ecov_years_full <- (year_start - 1L):(year_end + MSE_years)  # 1972..2052 (81)
stopifnot(length(ecov_years_full) == (n_model_years + 1L))

ecov_om <- list(
  label           = "bt_temp",
  mean            = matrix(c(env.dat_me$bt_temp, rep(NA_real_, MSE_years + 1L)), ncol = 1),
  logsigma        = "est_1",
  year            = ecov_years_full,
  use_obs         = matrix(1L, ncol = 1, nrow = length(ecov_years_full)),
  process_model   = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear", ncol = 1)
)

# Linear projection from 1972..2021 → fill 2022..2052
lm_fit   <- lm(bt_temp ~ Year, data = env.dat_me)
lm_sum   <- summary(lm_fit)
error_sd <- lm_sum[["sigma"]]

proj_df <- data.frame(year_proj = ecov_years_full)
set.seed(42)
proj_df$eps <- rnorm(nrow(proj_df), 0, error_sd)
proj_df$Temp_proj <- lm_sum$coefficients[1,1] + lm_sum$coefficients[2,1]*proj_df$year_proj + proj_df$eps

fill_idx <- which(ecov_years_full >= 2022L)
ecov_om_P <- ecov_om
ecov_om_P$mean[fill_idx, 1] <- proj_df$Temp_proj[fill_idx]

## ------------------------------
## B4) Build projecting OM input (unfitted)
## ------------------------------
om_input_P <- prepare_wham_input(
  basic_info  = basic_info,
  selectivity = sel3,
  M           = M,
  NAA_re      = NAA_re,
  ecov        = ecov_om_P,
  catch_info  = catch_info_use,
  index_info  = index_info_use,
  F           = F_info,
  age_comp    = "logistic-normal-pool0"
)

# Apply WAA for full horizon
om_input_P <- whamMSE::update_waa(om_input_P, waa_info = info$par_inputs$user_waa)

## 8.2 Copy parameter vectors from the fitted OM (keeps dynamics close to base fit)
om_input_P$par$Ecov_process_pars <- om_ecov$parList$Ecov_process_pars
om_input_P$par$Ecov_beta_R       <- om_ecov$parList$Ecov_beta_R
om_input_P$par$catch_paa_pars    <- om_ecov$parList$catch_paa_pars
om_input_P$par$index_paa_pars    <- om_ecov$parList$index_paa_pars
om_input_P$par$sel_repars        <- om_ecov$parList$sel_repars

# Initialize N1 from the historical OM
for (i in 1:n_regions) {
  om_input_P$par$log_N1[i, i, ] <- log(om_ecov$rep$N1)[i, i, ]
}

# Extend index usage/sigma into projections by holding last historical row
om_input_P$data$agg_index_sigma[1:hist_len, ] <- input_hist$data$agg_index_sigma
om_input_P$data$use_indices[1:hist_len, ]     <- input_hist$data$use_indices
om_input_P$data$use_index_paa[1:hist_len, ]   <- input_hist$data$use_index_paa
if (hist_len < n_proj_years) {
  for (iy in (hist_len + 1):n_proj_years) {
    om_input_P$data$agg_index_sigma[iy, ] <- input_hist$data$agg_index_sigma[hist_len, ]
    om_input_P$data$use_indices[iy, ]     <- input_hist$data$use_indices[hist_len, ]
    om_input_P$data$use_index_paa[iy, ]   <- input_hist$data$use_index_paa[hist_len, ]
  }
}

# Drop years where an index wasn’t used (historical pattern)
loc_zeros_agg <- lapply(seq_len(ncol(input_hist$data$use_indices)), function(j)
  which(input_hist$data$use_indices[, j] == 0))
loc_zeros_paa <- lapply(seq_len(ncol(input_hist$data$use_index_paa)), function(j)
  which(input_hist$data$use_index_paa[, j] == 0))

max_agg <- max(sapply(loc_zeros_agg, length))
max_paa <- max(sapply(loc_zeros_paa, length))
remove_agg_years1 <- if (max_agg > 0) {
  out <- matrix(NA_integer_, nrow = max_agg, ncol = length(loc_zeros_agg))
  for (j in seq_along(loc_zeros_agg)) if (length(loc_zeros_agg[[j]]))
    out[seq_along(loc_zeros_agg[[j]]), j] <- loc_zeros_agg[[j]]
  out
} else NULL
remove_paa_years1 <- if (max_paa > 0) {
  out <- matrix(NA_integer_, nrow = max_paa, ncol = length(loc_zeros_paa))
  for (j in seq_along(loc_zeros_paa)) if (length(loc_zeros_paa[[j]]))
    out[seq_along(loc_zeros_paa[[j]]), j] <- loc_zeros_paa[[j]]
  out
} else NULL

om_input_P <- whamMSE::update_input_index_info(
  om_input_P,
  agg_index_sigma = om_input_P$data$agg_index_sigma,
  index_Neff      = om_input_P$data$index_Neff,
  remove_agg = TRUE, remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
  remove_paa = TRUE, remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
)

## ------------------------------
## B5) Ecov residuals: historical + AR(1) projection
## ------------------------------
hist_ecov_len  <- hist_len + 1L          # includes lag-lead 1972
total_ecov_len <- n_model_years + 1L     # 81 = 1972..2052
stopifnot(length(ecov_years_full) == total_ecov_len)

# Ensure slot exists
if (is.null(om_input_P$par$Ecov_re)) {
  om_input_P$par$Ecov_re <- matrix(0, nrow = total_ecov_len, ncol = 1)
} else {
  om_input_P$par$Ecov_re <- matrix(om_input_P$par$Ecov_re, ncol = 1)
  if (nrow(om_input_P$par$Ecov_re) != total_ecov_len)
    om_input_P$par$Ecov_re <- matrix(0, nrow = total_ecov_len, ncol = 1)
}

# Historical residuals (1972..2022)
Ecov_re_vec <- as.vector(Ecov_re)
stopifnot(length(Ecov_re_vec) == hist_ecov_len)
om_input_P$par$Ecov_re[1:hist_ecov_len, 1] <- Ecov_re_vec

hist_ecov = om_input_P$par$Ecov_process_pars[1,1] + Ecov_re_vec
proj_ecov = proj_df$Temp_proj[fill_idx[-1]] 

om_input_P$par$Ecov_process_pars[1,1] = mean(c(hist_ecov,proj_ecov)) #  7.705472

om_input_P$par$Ecov_re = matrix(c(hist_ecov,proj_ecov) - mean(c(hist_ecov,proj_ecov)), ncol =  1)

# Make sure we do NOT re-simulate Ecov_re
if (!is.null(om_input_P$random) && "Ecov_re" %in% om_input_P$random) {
  om_input_P$random <- om_input_P$random[om_input_P$random != "Ecov_re"]
}
# If supported in your build, force no simulation of Ecov RE:
om_input_P$data$do_simulate_Ecov_re <- 0L

## ------------------------------
## B6) Compile unfitted OM and simulate a concrete dataset
## ------------------------------
unfitted_om_P <- fit_wham(
  om_input_P,
  do.fit = FALSE,
  do.brps = FALSE,
  MakeADFun.silent = TRUE
)

random_P <- unfitted_om_P$input$random
if (!is.null(random_P)) random_P <- setdiff(random_P, "Ecov_re")

set.seed(12345)
om_with_data_P <- update_om_fn(
  unfitted_om_P,
  seed   = 12345,
  random = random_P
)

# Quick checks (optional)
plot(om_with_data_P$input$data$Ecov_obs, type="l", col="red")
lines(om_with_data_P$rep$Ecov_x, col="blue")

## ------------------------------
## B7) Build EM-side Ecov scenarios
## ------------------------------
years_all <- ecov_years_full       # 1972..2052
lag_ecov  <- 1L
n_total   <- n_model_years + 1L    # includes lag-lead year
stopifnot(length(om_with_data_P$input$data$Ecov_obs) == n_total)

# Template: EM sees the OM-simulated observed Ecov
ecov_me <- list(
  label           = "bt_temp",
  mean            = as.matrix(om_with_data_P$input$data$Ecov_obs),
  logsigma        = "est_1",
  year            = years_all,
  use_obs         = matrix(1L, ncol = 1, nrow = n_total),
  process_model   = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear")
)

ecov_AR1 <- ecov_me

# Historical-average future (deterministic flat mean for t >= 2023)
hist_rows    <- (1 + lag_ecov):(hist_len + lag_ecov)        # rows for 1973..2022
hist_avg     <- mean(ecov_me$mean[hist_rows, 1], na.rm = TRUE)
proj_rows    <- (hist_len + lag_ecov + 1):n_total           # rows for 2023..end
ecov_HistAvg <- ecov_me
ecov_HistAvg$mean[proj_rows, 1] <- hist_avg

# No Ecov → Rec link (Ecov state still estimated)
ecov_NoLink <- ecov_AR1
ecov_NoLink$recruitment_how[] <- "none"

# Terminal-hold (future mean = last historical year’s mean)
ecov_TermHold <- ecov_HistAvg
id_last_hist  <- which(years_all == year_end)               # row for 2022
ecov_TermHold$mean[proj_rows, 1] <- ecov_TermHold$mean[id_last_hist, 1]

scenarios <- list(AR1 = ecov_AR1, HistAvg = ecov_HistAvg, NoLink = ecov_NoLink, TermHold = ecov_TermHold)

## ------------------------------
## B8) Run EMs under different Ecov scenarios
## ------------------------------
hcr <- list(hcr.type = 1, hcr.opts = list(use_FXSPR = TRUE, percentFXSPR = 75))

update_index_info_list <- list(
  agg_index_sigma   = om_input_P$data$agg_index_sigma,
  index_Neff        = om_input_P$data$index_Neff,
  remove_agg        = TRUE,  remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
  remove_paa        = TRUE,  remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
)
update_catch_info_list <- list(
  agg_catch_sigma = om_input_P$data$agg_catch_sigma,
  catch_Neff      = om_input_P$data$catch_Neff
)

M_em        <- M
sel_em      <- sel3
NAA_re_em   <- NAA_re
NAA_re_em$N1_model <- "equilibrium"   # IMPORTANT: EM starts from equilibrium (common for EMs)
age_comp_em <- "logistic-normal-pool0"

# Helper to run one EM
run_em <- function(ecov_obj, ecov_em_opts, seed = 12345L, label = "mod"){
  loop_through_fn(
    om = om_with_data_P,
    em_info = info,
    random  = random_P,
    M_em = M_em,
    sel_em = sel_em,
    NAA_re_em = NAA_re_em,
    ecov_em = ecov_obj,
    ecov_em_opts = ecov_em_opts,
    age_comp_em = age_comp_em,
    em.opt = list(separate.em = FALSE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
    update_index_info = update_index_info_list,
    update_catch_info = update_catch_info_list,
    assess_years = assess.years,
    assess_interval = assess.interval,
    base_years = base.years,
    year.use = hist_len,
    add.years = TRUE,
    seed = seed,
    hcr = hcr,
    save.sdrep = FALSE,
    save.last.em = TRUE,
    FXSPR_init = 0.5,
    do.brps = TRUE
  )
}

# mod1: EM knows true AR1 Ecov (best case / baseline)
set.seed(12345)
mod1 <- run_em(scenarios$AR1, 
               ecov_em_opts = NULL,
               seed = 12345L, 
               label = "mod1")
saveRDS(mod1, "mod1.RDS")

# mod2: EM uses historical-average future
proj_idx <- which(years_all >= (year_end + 1L))  # years >= 2023
set.seed(12345)
mod2 <- run_em(
  scenarios$HistAvg,
  ecov_em_opts = list(use_ecov_em = TRUE, period = proj_idx, lag = 1),
  seed = 12345L, label = "mod2"
)
saveRDS(mod2, "mod2.RDS")

# mod3: EM estimates Ecov but removes Ecov→Recruitment link
set.seed(12345)
mod3 <- run_em(
  scenarios$NoLink,
  ecov_em_opts = list(use_ecov_em = TRUE, period = proj_idx, lag = 1),
  seed = 12345L, label = "mod3"
)
saveRDS(mod3, "mod3.RDS")

# mod4: EM uses terminal-hold future Ecov
set.seed(12345)
mod4 <- run_em(
  scenarios$TermHold,
  ecov_em_opts = list(use_ecov_em = TRUE, period = proj_idx, lag = 1),
  seed = 12345L, label = "mod4"
)
saveRDS(mod4, "mod4.RDS")

## ------------------------------
## B9) Quick comparisons (feedback years only: 51..80)
## ------------------------------
par(mfrow = c(1,1))
yr_idx <- (hist_len + 1):n_model_years  # 51..80

plot(mod1$om$rep$SSB[yr_idx, ], type = "l", lwd = 2, main = "SSB (Feedback Years)")
lines(mod2$om$rep$SSB[yr_idx, ], col = "red",   lwd = 2)
lines(mod3$om$rep$SSB[yr_idx, ], col = "blue",  lwd = 2)
lines(mod4$om$rep$SSB[yr_idx, ], col = "green", lwd = 2)
legend("topleft",
       legend = c("AR1", "Average", "No_Ecov", "Last_Ecov"),
       col    = c("black", "red", "blue", "green"),
       lty = 1, lwd = 2, bty = "n", inset = 0.02)

plot(mod1$om$rep$pred_catch[yr_idx, ], type = "l", lwd = 2, main = "Catch (Feedback Years)")
lines(mod2$om$rep$pred_catch[yr_idx, ], col = "red",   lwd = 2)
lines(mod3$om$rep$pred_catch[yr_idx, ], col = "blue",  lwd = 2)
lines(mod4$om$rep$pred_catch[yr_idx, ], col = "green", lwd = 2)
legend("topleft",
       legend = c("AR1", "Average", "No_Ecov", "Last_Ecov"),
       col    = c("black", "red", "blue", "green"),
       lty = 1, lwd = 2, bty = "n", inset = 0.02)

plot(mod1$om$rep$Fbar[yr_idx, 1], type = "l", lwd = 2, main = "Fbar (Feedback Years)")
lines(mod2$om$rep$Fbar[yr_idx, 1], col = "red",   lwd = 2)
lines(mod3$om$rep$Fbar[yr_idx, 1], col = "blue",  lwd = 2)
lines(mod4$om$rep$Fbar[yr_idx, 1], col = "green", lwd = 2)
legend("topleft",
       legend = c("AR1", "Average", "No_Ecov", "Last_Ecov"),
       col    = c("black", "red", "blue", "green"),
       lty = 1, lwd = 2, bty = "n", inset = 0.02)

mod1$om$rep$Ecov_x
mod2$om$rep$Ecov_x

sum(as.numeric(mod1$converge_list)) # you should have 20! 
sum(as.numeric(mod2$converge_list)) # you should have 20!
sum(as.numeric(mod3$converge_list)) # you should have 20! 
sum(as.numeric(mod4$converge_list)) # you should have 20! 
