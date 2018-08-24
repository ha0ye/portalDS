#### setup ----

library(tidyverse)
library(lubridate)
library(rEDM)
library(parallel)

# devtools::install_github("weecology/portalr")
# devtools::install_github("weecology/LDATS")

source(here::here("R/dynamic_stability_functions.R"))
source(here::here("R/plotting_functions.R"))

#### plotting code ----

make_time_series_plots <- function()
{
    block <- readRDS(here::here("data/portal_block_50.RDS"))
    
    ## plot time series ----
    sp_list <- c("DM", "DS", "PP")
    palette <- viridis(NCOL(block)-1, option = "plasma")[match(sp_list, sort(colnames(select(block, -censusdate))))]
    portal_time_series <- 
        select(block, c("censusdate", sp_list)) %>%
        plot_time_series() + 
        scale_x_date(breaks = seq(from = as.Date("1980-01-01"), to = as.Date("2015-01-01"), by = "5 years"), 
                     date_labels = "%Y", expand = c(0.01, 0)) + 
        scale_color_manual(values = palette) + 
        theme(panel.background = element_rect(fill = "#AAAABB", color = NA), 
              panel.grid.major = element_line(color = "grey30", size = 1), 
              panel.grid.minor = element_line(color = "grey30", size = 1), 
              legend.key = element_rect(fill = "#AAAABB"))
    
    ggsave(here::here("figures/portal_time_series.pdf"),
           portal_time_series, width = 12, height = 6)
    ggsave(here::here("figures/portal_time_series_highlight.pdf"),
           add_regime_shift_highlight(portal_time_series), 
           width = 12, height = 6)
    
    return(portal_time_series)
}

portal_time_series <- make_time_series_plots()

make_eigenvalue_plots <- function()
{
    results_file <- here::here("output/portal_ds_results_50.RDS")
    results <- readRDS(results_file)
    
    dom_ev_plot <- results$eigenvalues %>%
        plot_eigenvalues() + 
        scale_y_continuous(limits = c(0.90, 1.03))
    ggsave(here::here("figures/portal_eigenvalues_50_dom.pdf"), 
           dom_ev_plot, width = 12, height = 6)
    
    two_ev_plot <- results$eigenvalues %>%
        plot_eigenvalues(num_values = 2) + 
        scale_y_continuous(limits = c(0.90, 1.03)) 
    ggsave(here::here("figures/portal_eigenvalues_50.pdf"), 
           two_ev_plot, width = 12, height = 6)
    
    two_ev_plot_complex <- results$eigenvalues %>%
        plot_eigenvalues(num_values = 2, highlight_complex = TRUE) + 
        scale_y_continuous(limits = c(0.90, 1.03)) 
    ggsave(here::here("figures/portal_eigenvalues_50_complex.pdf"), 
           two_ev_plot_complex, width = 12, height = 6)
    
    lower_date <- as.Date("1984-09-29")
    upper_date <- as.Date("1996-06-13")
    two_ev_plot_a <- results$eigenvalues %>%
        plot_eigenvalues(num_values = 2, highlight_complex = TRUE) %>% 
        add_regime_shift_highlight(lower_date, upper_date) + 
        scale_y_continuous(limits = c(0.90, 1.03))
    ggsave(here::here("figures/portal_eigenvalues_50_a.pdf"), 
           two_ev_plot_a, width = 12, height = 6)
    
    ggsave(here::here("figures/portal_eigenvalues_50_aa.pdf"), 
           cowplot::plot_grid(two_ev_plot_a + 
                                  scale_x_date(limits = as.Date(c("1979-09-22", "2015-03-26")), 
                                               breaks = seq(from = as.Date("1980-01-01"), 
                                                            to = as.Date("2015-01-01"), 
                                                            by = "5 years"), 
                                               date_labels = "%Y", expand = c(0.01, 0)), 
                              portal_time_series %>% 
                                  add_regime_shift_highlight(lower_date, upper_date) + 
                                  scale_x_date(limits = as.Date(c("1979-09-22", "2015-03-26")), 
                                               breaks = seq(from = as.Date("1980-01-01"), 
                                                            to = as.Date("2015-01-01"), 
                                                            by = "5 years"), 
                                               date_labels = "%Y", expand = c(0.01, 0)), 
                              nrow = 2, align = "v", axis = "lr"), 
           width = 12, height = 8)
    
    lower_date <- as.Date(c("1999-01-12", "2003-09-17", "2009-08-13"))
    upper_date <- as.Date(c("2000-01-12", "2004-09-17", "2011-01-15"))
    two_ev_plot_b <- results$eigenvalues %>%
        plot_eigenvalues(num_values = 2, highlight_complex = TRUE) %>%
        add_regime_shift_highlight(lower_date, upper_date) + 
        scale_y_continuous(limits = c(0.90, 1.03)) 
    ggsave(here::here("figures/portal_eigenvalues_50_b.pdf"), 
           two_ev_plot_b, width = 12, height = 6)
    
    ggsave(here::here("figures/portal_eigenvalues_50_bb.pdf"), 
           cowplot::plot_grid(two_ev_plot_b + 
                                  scale_x_date(limits = as.Date(c("1979-09-22", "2015-03-26")), 
                                               breaks = seq(from = as.Date("1980-01-01"), 
                                                            to = as.Date("2015-01-01"), 
                                                            by = "5 years"), 
                                               date_labels = "%Y", expand = c(0.01, 0)), 
                              portal_time_series %>%
                                  add_regime_shift_highlight(lower_date, upper_date) + 
                                  scale_x_date(limits = as.Date(c("1979-09-22", "2015-03-26")), 
                                               breaks = seq(from = as.Date("1980-01-01"), 
                                                            to = as.Date("2015-01-01"), 
                                                            by = "5 years"), 
                                               date_labels = "%Y", expand = c(0.01, 0)), 
                              nrow = 2, align = "v", axis = "lr"), 
           width = 12, height = 8)
}

make_eigenvalue_plots()

make_eigenvector_plots <- function()
{
    results <- readRDS("output/portal_ds_results_50.RDS")
    
    eigenvectors <- results$eigenvectors
    
    # extract vars
    non_null_idx <- first(which(!sapply(eigenvectors, is.null)))
    var_names <- rownames(eigenvectors[[non_null_idx]])
    var_idx <- grep("_0", var_names)
    var_names <- gsub("_0", "", var_names[var_idx])
    
    vector_scale <- function(v)
    {
        sum_sq <- sum(abs(v)^2)
        v / sqrt(sum_sq)
    }
    
    # make data.frame of eigenvector components
    v_df <- map_dfr(seq(eigenvectors), function(i) {
        v <- eigenvectors[[i]]
        if (is.null(v))
            return()
        out <- reshape2::melt(v[var_idx, ])
        out$t <- i
        return(out)
    }) %>% 
        rename(variable = Var1, rank = Var2) %>%
        mutate(censusdate = as.Date(names(eigenvectors)[t]), 
               variable = as.factor(var_names[variable]), 
               value = abs(Re(value))) %>%
        group_by(t, rank) %>%
        mutate(value = vector_scale(value)) %>%
        ungroup()
    
    lower_date <- as.Date(c("1999-01-12", "2003-09-17", "2009-08-13"))
    upper_date <- as.Date(c("2000-01-12", "2004-09-17", "2011-01-15"))
    
    v_df %>%
        filter(censusdate >= lower_date[1], censusdate <= upper_date[1], 
               rank == 1) %>%
        ggplot(aes(x = censusdate, y = value, color = variable)) + 
        scale_color_viridis(discrete = TRUE, option = "plasma") + 
        scale_x_date(expand = c(0.01, 0)) + 
        geom_line(size = NA) + 
        labs(x = "", y = "magnitude", color = "Species") + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.background = element_rect(fill = "#AAAABB", color = NA), 
              panel.grid.major = element_line(color = "grey30", size = 1), 
              panel.grid.minor = element_line(color = "grey30", size = 1), 
              legend.key = element_rect(fill = "#AAAABB")) + 
        guides(color = guide_legend(override.aes = list(size = 1)))
    ggsave(here::here("figures/portal_eigenvector_null.pdf"),
           width = 12, height = 6)
    
    v_df %>%
        filter(censusdate >= lower_date[1], censusdate <= upper_date[1], 
               rank == 1) %>%
        ggplot(aes(x = censusdate, y = value, color = variable)) + 
        scale_color_viridis(discrete = TRUE, option = "plasma") + 
        scale_x_date(expand = c(0.01, 0)) + 
        geom_line(size = 1) + 
        labs(x = "", y = "magnitude", color = "Species") + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.background = element_rect(fill = "#AAAABB", color = NA), 
              panel.grid.major = element_line(color = "grey30", size = 1), 
              panel.grid.minor = element_line(color = "grey30", size = 1), 
              legend.key = element_rect(fill = "#AAAABB"))
    ggsave(here::here("figures/portal_eigenvector_a.pdf"),
           width = 12, height = 6)
    
    v_df %>%
        filter(censusdate >= lower_date[2], censusdate <= upper_date[2], 
               rank == 1) %>%
        ggplot(aes(x = censusdate, y = value, color = variable)) + 
        scale_color_viridis(discrete = TRUE, option = "plasma") + 
        scale_x_date(expand = c(0.01, 0)) + 
        geom_line(size = 1) + 
        labs(x = "", y = "magnitude", color = "Species") + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.background = element_rect(fill = "#AAAABB", color = NA), 
              panel.grid.major = element_line(color = "grey30", size = 1), 
              panel.grid.minor = element_line(color = "grey30", size = 1), 
              legend.key = element_rect(fill = "#AAAABB"))
    ggsave(here::here("figures/portal_eigenvector_b.pdf"),
           width = 12, height = 6)
    
    v_df %>%
        filter(censusdate >= lower_date[3], censusdate <= upper_date[3], 
               rank == 1) %>%
        ggplot(aes(x = censusdate, y = value, color = variable)) + 
        scale_color_viridis(discrete = TRUE, option = "plasma") + 
        scale_x_date(expand = c(0.01, 0)) + 
        geom_line(size = 1) + 
        labs(x = "", y = "magnitude", color = "Species") + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.background = element_rect(fill = "#AAAABB", color = NA), 
              panel.grid.major = element_line(color = "grey30", size = 1), 
              panel.grid.minor = element_line(color = "grey30", size = 1), 
              legend.key = element_rect(fill = "#AAAABB"))
    ggsave(here::here("figures/portal_eigenvector_c.pdf"),
           width = 12, height = 6)
}

make_eigenvector_plots()