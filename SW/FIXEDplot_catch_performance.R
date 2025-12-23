FIXEDplot_catch_performance <- function (mods, is.nsim, main.dir, sub.dir, var = "Catch", width = 10, 
          height = 7, dpi = 300, col.opt = "D", method = NULL, outlier.opt = NA, 
          plot.style = "median_iqr", show.whisker = TRUE, use.n.years = NULL, 
          new_model_names = NULL, base.model = NULL) 
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  res <- NULL
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    res <- lapply(seq_along(mods), function(i) {
      data.frame(Catch = mods[[i]]$om$rep$pred_catch, Model = paste0("Model", 
                                                                     i), Year = Years, Realization = 1) %>% tail(use.n.years)
    }) %>% bind_rows()
  }
  else {
    Years <- mods[[1]][[1]]$om$years
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        data.frame(Catch = mods[[r]][[m]]$om$rep$pred_catch, 
                   Model = paste0("Model", m), Year = Years, Realization = r) %>% 
          tail(use.n.years)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  res <- res %>% rowwise() %>% mutate(Catch_Global = sum(c_across(starts_with("Catch")), 
                                                         na.rm = TRUE)) %>% ungroup()
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
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", 
                      values_to = var)
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>% rename(base_val = !!sym(var)) %>% 
      select(Realization, Year, Label, base_val)
    res <- left_join(res, base_df, by = c("Realization", 
                                          "Year", "Label")) %>% mutate(`:=`(!!var, (!!sym(var))/base_val - 
                                                                              1))
  }
  if (!is.null(method)) {
    res <- res %>% group_by(Model, Realization, Label) %>% 
      summarise(`:=`(!!var, if (method == "mean") 
        mean(!!sym(var))
        else median(!!sym(var))), .groups = "drop")
  }
  if (plot.style == "boxplot") {
    p <- ggplot(res, aes(x = Model, y = Catch, color = Model)) + 
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) + 
      facet_grid(Label ~ ., scales = "free") + scale_fill_viridis_d(option = col.opt) + 
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", 
                                                             var, " vs ", base.model)), ": Last ", use.n.years, 
                     " Years")) + ylab(ifelse(is.null(base.model), 
                                              var, paste0("Relative ", var, " Difference"))) + 
      xlab("Model") + theme_bw()
  }
  else if (plot.style == "median_iqr") {
    res_summary <- res %>% group_by(Model, Label) %>% summarise(q1 = quantile(!!sym(var), 
                                                                              0.25), med = median(!!sym(var)), q3 = quantile(!!sym(var), 
                                                                                                                             0.75), iqr = q3 - q1, .groups = "drop") %>% mutate(x = as.numeric(factor(Model)), 
                                                                                                                                                                                ymin = if (show.whisker) 
                                                                                                                                                                                  q1 - 1.5 * iqr
                                                                                                                                                                                else NA_real_, ymax = if (show.whisker) 
                                                                                                                                                                                  q3 + 1.5 * iqr
                                                                                                                                                                                else NA_real_)
    res_limits <- res %>% group_by(Model, Label) %>% summarise(min_val = min(!!sym(var), 
                                                                             na.rm = TRUE), max_val = max(!!sym(var), na.rm = TRUE), 
                                                               .groups = "drop")
    res_summary <- left_join(res_summary, res_limits, by = c("Model", 
                                                             "Label")) %>% mutate(ymin = pmax(ymin, min_val), 
                                                                                  ymax = pmin(ymax, max_val))
    p <- ggplot(res_summary, aes(x = x, color = Model)) + 
      {
        if (show.whisker) 
          geom_segment(aes(x = x, xend = x, y = ymin, 
                           yend = q1), color = 'black')
      } + {
        if (show.whisker) 
          geom_segment(aes(x = x, xend = x, y = q3, yend = ymax), color = 'black')
      } + geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, 
                        ymax = q3, fill = Model), color = 'black', linewidth = 0.8) + 
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, 
                       yend = med), color = 'black', linewidth = 0.8) + 
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) + 
      facet_grid(Label ~ ., scales = "free") + scale_fill_viridis_d(option = col.opt) + 
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", 
                                                             var, " vs ", base.model)), ": Last ", use.n.years, 
                     " Years")) + ylab(ifelse(is.null(base.model), 
                                              var, paste0("Relative ", var, " Difference"))) + 
      xlab("Model") + theme_bw()
  }
  else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  plot_name <- paste0(var, ifelse(is.null(base.model), "", 
                                  "_Relative"), "_last_", use.n.years, "_years.PNG")
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  if (!file.exists(new_sub_dir)) {
    dir.create(new_sub_dir)
  }
  ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, 
         height = height, dpi = dpi)
  return(p)
}