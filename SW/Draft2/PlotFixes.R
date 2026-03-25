args <- commandArgs(trailingOnly = TRUE)

outdir <- args[which(args == "--outdir") + 1]

library(ggridges)
library(wham)
library(whamMSE)
library(dplyr)
library(ggplot2)
library(ggtern)



source('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/SpiderPlotFix.R')

filepath <- '/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/5Yr/SensitivityTest/Outputs'


# PROJECTIONS
mods_cleaned <- readRDS(paste0(filepath,'/MSEmods_50yr.rds'))
mods_cleaned[sapply(mods_cleaned, is.null)] <- NULL

mods <- mods_cleaned
main_dir <- filepath
main.dir <- main_dir
output_dir <- "FigFix"
sub.dir <- output_dir
output_format <- c("png") # or html or png
width <- 10
height <- 7
dpi <- 300
col.opt <- "D"
new_model_names <- c("AR(1)", "No Ecov", "Historical Avg", "Recent Window", "Recent Trend", "Terminal Year", "Bad Projection","Good Projection")
#new_model_names <- c("AR(1)", "No Ecov", "Low", "Medium Low", "Medium", "High")
base.model <- "Mod_AR1Ecov"
start.years <- 51
use.n.years.first <- 1
use.n.years.last <- 10
is.nsim <- if (!is.list(mods[[1]][[1]][[1]])) FALSE else TRUE
method <- "median"

#whamMSE:::plot_model_performance_radar(mods_cleaned, is.nsim, '.', output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)


tryCatch(
FIXEDplot_model_performance_radar(mods_cleaned, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)
)
tryCatch(
FIXEDplot_model_performance_radar2(mods_cleaned, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)
)
tryCatch(
FIXEDplot_model_performance_radar3(mods_cleaned, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)
)





mods <- mods_cleaned
mods <- mods_cleaned
main_dir <- filepath
main.dir <- main_dir
output_dir <- "FigFix"
sub.dir <- output_dir
output_format <- c("png") # or html or png
width <- 10
height <- 7
dpi <- 300
col.opt <- "D"
base.model <- "Mod_AR1Ecov"
start.years <- 51
use.n.years.first <- 1
use.n.years.last <- 10
is.nsim <- if (!is.list(mods[[1]][[1]][[1]])) FALSE else TRUE
method <- "median"

source('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/FIXEDplot_ssb_performance.R')
tryCatch(
FIXEDplot_ssb_performance(mods, is.nsim, main.dir, sub.dir, var = "SSB", width, 
                          height, dpi, col.opt, method, outlier.opt = NA, 
                          plot.style = "median_iqr", show.whisker = TRUE, use.n.years = 10, 
                          new_model_names)
)

source('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/FIXEDplot_fbar_performance.R')
tryCatch(
FIXEDplot_fbar_performance(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, method, outlier.opt = NA, plot.style = "median_iqr", show.whisker = TRUE, f.ymin = NULL, f.ymax = NULL, use.n.years = 10, new_model_names)
)

source('/work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/FIXEDplot_catch_performance.R')
tryCatch(
FIXEDplot_catch_performance(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt = NA, plot.style = "median_iqr", show.whisker = TRUE, use.n.years = 10, new_model_names)
)
