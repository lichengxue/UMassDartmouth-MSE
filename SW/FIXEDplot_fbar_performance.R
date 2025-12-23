FIXEDplot_fbar_performance <- function (mods, is.nsim, main.dir, sub.dir, var = "Fbar", width = 10, 
          height = 7, dpi = 300, col.opt = "D", method = NULL, outlier.opt = NA, 
          plot.style = "median_iqr", show.whisker = TRUE, f.ymin = NULL, 
          f.ymax = NULL, use.n.years = NULL, new_model_names = NULL, 
          base.model = NULL) 
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  make_plot_data <- function(index_range, label_prefix) {
    if (!is.nsim) {
      Years <- mods[[1]]$om$years
      n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
      n_regions <- mods[[1]]$om$input$data$n_regions[1]
      lapply(seq_along(mods), function(i) {
        tmp <- mods[[i]]$om$rep$Fbar[, index_range, drop = FALSE]
        tmp <- as.data.frame(tmp)
        names(tmp) <- paste0(label_prefix, seq_along(index_range))
        tmp$Model <- paste0("Model", i)
        tmp$Year <- Years
        tmp$Realization <- 1
        tmp <- tail(tmp, use.n.years)
      }) %>% bind_rows()
    }
    else {
      Years <- mods[[1]][[1]]$om$years
      lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          tmp <- mods[[r]][[m]]$om$rep$Fbar[, index_range, 
                                            drop = FALSE]
          tmp <- as.data.frame(tmp)
          names(tmp) <- paste0(label_prefix, seq_along(index_range))
          tmp$Model <- paste0("Model", m)
          tmp$Year <- Years
          tmp$Realization <- r
          tmp <- tail(tmp, use.n.years)
        }) %>% bind_rows()
      }) %>% bind_rows()
    }
  }
  if (!is.nsim) {
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
  }
  else {
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
  }
  res_fleet <- make_plot_data(1:n_fleets, "Fleet_")
  res_region <- make_plot_data((n_fleets + 1):(n_fleets + n_regions), 
                               "Region_")
  res_global <- make_plot_data(n_fleets + n_regions + 1, "Global")
  plot_and_save <- function(res, title, ylab_text, filename) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(res$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      res$Model <- factor(res$Model, levels = paste0("Model", 
                                                     seq_along(new_model_names)), labels = new_model_names)
      if (!is.null(base.model)) {
        if (!(base.model %in% new_model_names)) {
          warning("base.model does not match any of the new_model_names.")
        }
      }
    }
    res_long <- pivot_longer(res, cols = starts_with(c("Fleet_", 
                                                       "Region_", "Global")), names_to = "Label", values_to = "Fbar")
    if (!is.null(base.model)) {
      base_df <- res_long %>% filter(Model == base.model) %>% 
        rename(base_val = Fbar) %>% select(Realization, 
                                           Year, Label, base_val)
      res_long <- left_join(res_long, base_df, by = c("Realization", 
                                                      "Year", "Label")) %>% mutate(Fbar = Fbar/base_val - 
                                                                                     1)
    }
    if (!is.null(base.model)) {
      if (!is.null(f.ymin)) 
        y1 = f.ymin
      else y1 = -1
      if (!is.null(f.ymax)) 
        y2 = f.ymax
      else y2 = 2
    }
    else {
      if (!is.null(f.ymin)) 
        y1 = f.ymin
      else y1 = 0
      if (!is.null(f.ymax)) 
        y2 = f.ymax
      else y2 = 2
    }
    if (!is.null(method)) {
      res_long <- res_long %>% group_by(Model, Realization, 
                                        Label) %>% summarise(`:=`(!!var, if (method == 
                                                                             "mean") 
                                          mean(!!sym(var))
                                          else median(!!sym(var))), .groups = "drop")
    }
    if (plot.style == "boxplot") {
      p <- ggplot(res_long, aes(x = Model, y = Fbar, color = Model)) + 
        geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) + 
        coord_cartesian(ylim = c(y1, y2)) + facet_grid(Label ~ 
                                                         ., scales = "free") + scale_fill_viridis_d(option = col.opt) + 
        ggtitle(paste0(ifelse(is.null(base.model), title, 
                              paste0("Relative ", title, " vs ", base.model)), 
                       ": Last ", use.n.years, " Years")) + ylab(ifelse(is.null(base.model), 
                                                                        ylab_text, "Relative Difference in Fbar")) + 
        xlab("Model") + theme_bw()
    }
    else if (plot.style == "median_iqr") {
      res_summary <- res_long %>% group_by(Model, Label) %>% 
        summarise(q1 = quantile(!!sym(var), 0.25), med = median(!!sym(var)), 
                  q3 = quantile(!!sym(var), 0.75), iqr = q3 - 
                    q1, .groups = "drop") %>% mutate(x = as.numeric(factor(Model)), 
                                                     ymin = if (show.whisker) 
                                                       q1 - 1.5 * iqr
                                                     else NA_real_, ymax = if (show.whisker) 
                                                       q3 + 1.5 * iqr
                                                     else NA_real_)
      res_limits <- res_long %>% group_by(Model, Label) %>% 
        summarise(min_val = min(!!sym(var), na.rm = TRUE), 
                  max_val = max(!!sym(var), na.rm = TRUE), .groups = "drop")
      res_summary <- left_join(res_summary, res_limits, 
                               by = c("Model", "Label")) %>% mutate(ymin = pmax(ymin, 
                                                                                min_val), ymax = pmin(ymax, max_val))
      p <- ggplot(res_summary, aes(x = x, color = Model)) + 
        {
          if (show.whisker) 
            geom_segment(aes(x = x, xend = x, y = ymin, 
                             yend = q1), color = 'black')
        } + {
          if (show.whisker) 
            geom_segment(aes(x = x, xend = x, y = q3, yend = ymax), color = 'black')
        } + geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, 
                          ymin = q1, ymax = q3, fill = Model), color = 'black', 
                      linewidth = 0.8) + geom_segment(aes(x = x - 0.3, 
                                                          xend = x + 0.3, y = med, yend = med), color = 'black', 
                                                      linewidth = 0.8) + scale_x_continuous(breaks = res_summary$x, 
                                                                                            labels = res_summary$Model) + facet_grid(Label ~ 
                                                                                                                                       ., scales = "free") + scale_fill_viridis_d(option = col.opt) + 
        ggtitle(paste0(ifelse(is.null(base.model), title, 
                              paste0("Relative ", title, " vs ", base.model)), 
                       ": Last ", use.n.years, " Years")) + ylab(ifelse(is.null(base.model), 
                                                                        ylab_text, "Relative Difference in Fbar")) + 
        xlab("Model") + theme_bw()
    }
    else {
      stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
    }
    plot_name <- paste0(filename, ifelse(is.null(base.model), 
                                         "", "_Relative"), ".PNG")
    new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
    if (!file.exists(new_sub_dir)) {
      dir.create(new_sub_dir)
    }
    ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, 
           height = height, dpi = dpi)
    return(p)
  }
  p_fleet <- plot_and_save(res_fleet, "Fbar by Fleet", "Fbar", 
                           paste0(var, "_fleet_last_", use.n.years, "_years"))
  p_region <- plot_and_save(res_region, "Fbar by Region", "Fbar", 
                            paste0(var, "_region_last_", use.n.years, "_years"))
  p_global <- plot_and_save(res_global, "Global Fbar", "Fbar", 
                            paste0(var, "_global_last_", use.n.years, "_years"))
  return(list(fleet = p_fleet, region = p_region, global = p_global))
}