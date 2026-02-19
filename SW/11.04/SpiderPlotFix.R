FIXEDplot_model_performance_radar <- function (mods, is.nsim, main.dir, sub.dir, width = 10, height = 10, 
          dpi = 300, col.opt = "D", method = NULL, use.n.years.first = 5, 
          use.n.years.last = 5, start.years = 1, new_model_names = NULL) 
{
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    results <- list()
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        fbar_ts <- rep$Fbar[, n_fleets + n_regions + 
                              1]
        if (is.null(method)) 
          method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + 
                                                               use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
          tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, 
                                           use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
        }
        else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
          tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
        }
      }
      for (v in c("Catch_first", "SSB_first", "Catch_last", 
                  "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        }
        else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]]))/range_val
        }
      }
      for (v in c("Fbar_first", "Fbar_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        }
        else {
          norm_f <- (tmp[[v]] - min(tmp[[v]]))/range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      results[[r]] <- tmp
    }
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>% group_by(Model) %>% summarise(across(-Realization, 
                                                                       median), .groups = "drop")
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model, 
                                    levels = paste0("Model", seq_along(new_model_names)), 
                                    labels = new_model_names)
    }
    plot_df <- as.data.frame(scores_median)
    rownames(plot_df) <- plot_df$Model
    plot_df$Model <- NULL
    plot_df <- as.data.frame(t(plot_df))
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                     plot_df)
  }
  else {
    n_models <- length(mods)
    results <- list()
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
      if (is.null(method)) 
        method = "median"
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
        tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
      }
      else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + 
                                                       use.n.years.first - 1)])
        tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
      }
    }
    for (v in c("Catch_first", "SSB_first", "Catch_last", 
                "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      }
      else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]]))/range_val
      }
    }
    for (v in c("Fbar_first", "Fbar_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      }
      else {
        norm_f <- (tmp[[v]] - min(tmp[[v]]))/range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    results[[1]] <- tmp
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>% group_by(Model) %>% summarise(across(-Realization, 
                                                                       median), .groups = "drop")
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model, 
                                    levels = paste0("Model", seq_along(new_model_names)), 
                                    labels = new_model_names)
    }
    plot_df <- as.data.frame(scores_median)
    rownames(plot_df) <- plot_df$Model
    plot_df$Model <- NULL
    plot_df <- as.data.frame(t(plot_df))
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                     plot_df)
  }
  if (ncol(plot_df) < 3) {
    message("Radar chart needs at least 3 models. Showing barplot instead.")
    return(invisible(NULL))
  }
  colors <- viridisLite::viridis(n = ncol(plot_df), option = col.opt) ####CHANGED COLORS
  new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
  if (!file.exists(new_sub_dir)) {
    dir.create(new_sub_dir)
  }
  my_legend_labels <- new_model_names #c(expression(Catch[ST]), expression(SSB[ST]), 
                        # expression(F[ST]), expression(Catch[LT]), expression(SSB[LT]), 
                        # expression(F[LT])) ##### HERE ########
  output_file <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot", 
                           "model_performance_radar.png")
  png(filename = output_file, width = width, height = height, 
      units = "in", res = dpi)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(1, 1, 1, 1))
  ##### HERE ########
  plot_df <- plot_df[-c(1:2),]
  plot_df <- data.frame(t(plot_df))
  plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
          plot_df)
  ##### HERE ########
  radarchart(plot_df, axistype = 4, pcol = colors, plwd = 3, 
             plty = 1:length(colors), cglcol = "grey80", cglty = 1, 
             axislabcol = "grey30", vlcex = 1.2)
  par(mar = c(1, 1, 1, 1))
  plot.new()
  legend("center", legend = my_legend_labels, col = colors, 
         lty = 1:length(colors), lwd = 3, cex = 0.9, y.intersp = 1.5)
  dev.off()
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(1, 1, 2, 1))
  radarchart(plot_df, axistype = 4, pcol = colors, plwd = 3, 
             plty = 1:length(colors), cglcol = "grey80", cglty = 1, 
             axislabcol = "grey30", vlcex = 1.2)
  par(mar = c(1, 1, 2, 1))
  plot.new()
  legend("center", legend = my_legend_labels, col = colors, 
         lty = 1:length(colors), lwd = 3, cex = 0.9, y.intersp = 2, 
         bty = "n")
}

FIXEDplot_model_performance_radar2 <- function (mods, is.nsim, main.dir, sub.dir, width = 10, height = 10, 
          dpi = 300, col.opt = "D", method = NULL, use.n.years.first = 5, 
          use.n.years.last = 5, start.years = 1, new_model_names = NULL) 
{
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  library(viridisLite)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs)/sum(values[-length(values)])
    return(aacv)
  }
  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    results <- list()
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
        ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
        fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
        if (is.null(method)) 
          method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + 
                                                               use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, 
                                           use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }
        else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }
      }
      for (v in c("Catch_first", "SSB_first", "Catch_last", 
                  "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        }
        else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]]))/range_val
        }
      }
      for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        }
        else {
          norm_f <- (tmp[[v]] - min(tmp[[v]]))/range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      results[[r]] <- tmp
    }
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>% group_by(Model) %>% summarise(across(-Realization, 
                                                                       median), .groups = "drop")
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model, 
                                    levels = paste0("Model", seq_along(new_model_names)), 
                                    labels = new_model_names)
    }
    plot_df <- as.data.frame(scores_median)
    rownames(plot_df) <- plot_df$Model
    plot_df$Model <- NULL
    plot_df <- as.data.frame(t(plot_df))
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                     plot_df)
  }
  else {
    n_models <- length(mods)
    results <- list()
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
      ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
      fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
      if (is.null(method)) 
        method = "median"
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
      else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + 
                                                       use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
    }
    for (v in c("Catch_first", "SSB_first", "Catch_last", 
                "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      }
      else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]]))/range_val
      }
    }
    for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      }
      else {
        norm_f <- (tmp[[v]] - min(tmp[[v]]))/range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    results[[1]] <- tmp
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>% group_by(Model) %>% summarise(across(-Realization, 
                                                                       median), .groups = "drop")
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model, 
                                    levels = paste0("Model", seq_along(new_model_names)), 
                                    labels = new_model_names)
    }
    plot_df <- as.data.frame(scores_median)
    rownames(plot_df) <- plot_df$Model
    plot_df$Model <- NULL
    plot_df <- as.data.frame(t(plot_df))
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                     plot_df)
  }
  if (ncol(plot_df) < 3) {
    message("Radar chart needs at least 3 models. Showing barplot instead.")
    return(invisible(NULL))
  }
  colors <- viridisLite::viridis(n = ncol(plot_df), option = col.opt) #FIXED HERE FROM NROW(PLOTDF)-2
  new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
  if (!file.exists(new_sub_dir)) {
    dir.create(new_sub_dir)
  }
  rownames(plot_df)
  my_legend_labels <- new_model_names#c(expression(Catch[ST]), expression(SSB[ST]), 
                        # expression(Catch[LT]), expression(SSB[LT]), expression(Catch[AAV]), 
                        # expression(SSB[AAV]), expression(F[AAV]))
  output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", 
                                     "model_performance_radar2.png"))
  png(filename = output_file, width = width, height = height, 
      units = "in", res = dpi)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(1, 1, 1, 1))
  ##### HERE ########
  plot_df <- plot_df[-c(1:2),]
  plot_df <- data.frame(t(plot_df))
  plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                   plot_df)
  ##### HERE ########
  radarchart(plot_df, axistype = 4, pcol = colors, plwd = 3, 
             plty = 1:length(colors), cglcol = "grey80", cglty = 1, 
             axislabcol = "grey30", vlcex = 1.2)
  par(mar = c(1, 1, 1, 1))
  plot.new()
  legend("center", legend = my_legend_labels, col = colors, 
         lty = 1:length(colors), lwd = 3, cex = 0.9, y.intersp = 1.5)
  dev.off()
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(1, 1, 2, 1))
  radarchart(plot_df, axistype = 4, pcol = colors, plwd = 3, 
             plty = 1:length(colors), cglcol = "grey80", cglty = 1, 
             axislabcol = "grey30", vlcex = 1.2)
  par(mar = c(1, 1, 2, 1))
  plot.new()
  legend("center", legend = my_legend_labels, col = colors, 
         lty = 1:length(colors), lwd = 3, cex = 0.9, y.intersp = 2, 
         bty = "n")
}

FIXEDplot_model_performance_radar3 <- function (mods, is.nsim, main.dir, sub.dir, width = 10, height = 10, 
          dpi = 300, col.opt = "D", method = NULL, use.n.years.first = 5, 
          use.n.years.last = 5, start.years = 1, new_model_names = NULL) 
{
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  library(viridisLite)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs)/sum(values[-length(values)])
    return(aacv)
  }
  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    results <- list()
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
        ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
        fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
        if (is.null(method)) 
          method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + 
                                                               use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, 
                                           use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }
        else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }
        if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
          message("Biological Reference Point has not been calculated internally!")
          return(invisible(NULL))
        }
        tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
        tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
        temp <- tmp1/tmp2
        temp1 <- temp[start.years:(start.years + use.n.years.first - 
                                     1)]
        temp2 <- tail(temp, use.n.years.last)
        tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
        tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
      }
      for (v in c("Catch_first", "SSB_first", "Catch_last", 
                  "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        }
        else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]]))/range_val
        }
      }
      for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv", 
                  "prob_first", "prob_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        }
        else {
          norm_f <- (tmp[[v]] - min(tmp[[v]]))/range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      results[[r]] <- tmp
    }
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>% group_by(Model) %>% summarise(across(-Realization, 
                                                                       median), .groups = "drop")
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model, 
                                    levels = paste0("Model", seq_along(new_model_names)), 
                                    labels = new_model_names)
    }
    plot_df <- as.data.frame(scores_median)
    rownames(plot_df) <- plot_df$Model
    plot_df$Model <- NULL
    plot_df <- as.data.frame(t(plot_df))
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                     plot_df)
  }
  else {
    n_models <- length(mods)
    results <- list()
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
      ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
      fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
      if (is.null(method)) 
        method = "median"
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + 
                                                             use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + 
                                                         use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
      else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + 
                                                           use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + 
                                                       use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
    }
    if (is.null(rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
    tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
    tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
    temp <- tmp1/tmp2
    temp1 <- temp[start.years:(start.years + use.n.years.first - 
                                 1)]
    temp2 <- tail(temp, use.n.years.last)
    tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
    tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
    for (v in c("Catch_first", "SSB_first", "Catch_last", 
                "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      }
      else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]]))/range_val
      }
    }
    for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv", "prob_first", 
                "prob_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      }
      else {
        norm_f <- (tmp[[v]] - min(tmp[[v]]))/range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    results[[1]] <- tmp
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>% group_by(Model) %>% summarise(across(-Realization, 
                                                                       median), .groups = "drop")
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model, 
                                    levels = paste0("Model", seq_along(new_model_names)), 
                                    labels = new_model_names)
    }
    plot_df <- as.data.frame(scores_median)
    rownames(plot_df) <- plot_df$Model
    plot_df$Model <- NULL
    plot_df <- as.data.frame(t(plot_df))
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                     plot_df)
  }
  if (ncol(plot_df) < 3) {
    message("Radar chart needs at least 3 models. Showing barplot instead.")
    return(invisible(NULL))
  }
  colors <- viridisLite::viridis(n = ncol(plot_df), option = col.opt) #CHANGED HERE FROM Nnrow(plot_df) - 2
  new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
  if (!file.exists(new_sub_dir)) {
    dir.create(new_sub_dir)
  }
  rownames(plot_df)
  my_legend_labels <- new_model_names#c(expression(Catch[ST]), expression(SSB[ST]), 
                        # expression(Catch[LT]), expression(SSB[LT]), expression(Catch[AAV]), 
                        # expression(SSB[AAV]), expression(F[AAV]), expression(P(F[ST] > 
                        #                                                          F[MSY])), expression(P(F[LT] > F[MSY])))
  output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", 
                                     "model_performance_radar3.png"))
  png(filename = output_file, width = width, height = height, 
      units = "in", res = dpi)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(1, 1, 1, 1))
  ##### HERE ########
  plot_df <- plot_df[-c(1:2),]
  plot_df <- data.frame(t(plot_df))
  plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), 
                   plot_df)
  ##### HERE ########
  radarchart(plot_df, axistype = 4, pcol = colors, plwd = 3, 
             plty = 1:length(colors), cglcol = "grey80", cglty = 1, 
             axislabcol = "grey30", vlcex = 1.2)
  par(mar = c(1, 1, 1, 1))
  plot.new()
  legend("center", legend = my_legend_labels, col = colors, 
         lty = 1:length(colors), lwd = 3, cex = 0.9, y.intersp = 1.5)
  dev.off()
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(1, 1, 2, 1))
  radarchart(plot_df, axistype = 4, pcol = colors, plwd = 3, 
             plty = 1:length(colors), cglcol = "grey80", cglty = 1, 
             axislabcol = "grey30", vlcex = 1.2)
  par(mar = c(1, 1, 2, 1))
  plot.new()
  legend("center", legend = my_legend_labels, col = colors, 
         lty = 1:length(colors), lwd = 3, cex = 0.9, y.intersp = 2, 
         bty = "n")
}
