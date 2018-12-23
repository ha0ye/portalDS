library(portalDS)
library(dplyr)

#### plotting code ----

make_time_series_plots <- function(block = readRDS(here::here("data/portal_block_50.RDS")),
                                   sp_list)
{
    ## plot time series ----
    sp_list <- c("DM", "DS", "PP")
    palette <- viridis(NCOL(block) - 1, option = "plasma")[
        match(sp_list, sort(colnames(select(block, -censusdate))))]
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

    # empty plot
    empty_plot <- results$eigenvectors %>%
        plot_eigenvectors(line_size = NA) +
        labs(x = "", y = "magnitude", color = "Species")
    ggsave(here::here("figures/portal_eigenvector_null.pdf"),
           empty_plot, width = 12, height = 6)

    # base plot of eigenvectors
    base_plot <- results$eigenvectors %>%
        plot_eigenvectors() +
        labs(x = "", y = "magnitude", color = "Species")

    ggsave(here::here("figures/portal_eigenvector_figure.pdf"),
           add_regime_shift_highlight(base_plot) +
               scale_x_date(limits = as.Date(c("1997-01-01", "2015-01-01")),
                            expand = c(0.01, 0)),
           width = 32.5/3, height = 12.5/3)

    lower_date <- as.Date(c("1999-01-12", "2003-09-17", "2009-08-13"))
    upper_date <- as.Date(c("2000-01-12", "2004-09-17", "2011-01-15"))
    ggsave(here::here("figures/portal_eigenvector_a.pdf"),
           base_plot + scale_x_date(limits = c(lower_date[1], upper_date[1]),
                                    expand = c(0.01, 0)),
           width = 12, height = 6)

    ggsave(here::here("figures/portal_eigenvector_b.pdf"),
           base_plot + scale_x_date(limits = c(lower_date[2], upper_date[2]),
                                    expand = c(0.01, 0)),
           width = 12, height = 6)

    ggsave(here::here("figures/portal_eigenvector_c.pdf"),
           base_plot + scale_x_date(limits = c(lower_date[3], upper_date[3]),
                                    expand = c(0.01, 0)),
           width = 12, height = 6)
}

make_eigenvector_plots()

make_combined_network <- function(plot_file = NULL)
{
    portal_network <- here::here("output/portal_ds_results.RDS") %>%
        readRDS() %>%
        .$ccm_links %>%
        plot_network()

    portal_network_50 <- here::here("output/portal_ds_results_50.RDS") %>%
        readRDS() %>%
        .$ccm_links %>%
        plot_network(palette = portal_network$palette,
                     existing_graph = portal_network$graph)

    portal_network_33 <- here::here("output/portal_ds_results_33.RDS") %>%
        readRDS() %>%
        .$ccm_links %>%
        plot_network(palette = portal_network$palette,
                     existing_graph = portal_network$graph)

    combined_network_plot <- plot_grid(portal_network$plot,
                                       portal_network_50$plot,
                                       portal_network_33$plot,
                                       labels = c("all species",
                                                  "species present >= 50%",
                                                  "species present >= 33%"),
                                       align = "v", axis = "l",
                                       ncol = 1)

    if (!is.null(plot_file))
    {
        ggsave(plot_file, combined_network_plot, width = 6, height = 12)
    }
    return(combined_network_plot)
}

