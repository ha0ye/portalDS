library(portalDS)
library(ggplot2)
library(dplyr)
library(viridis)

block <- readRDS(here::here("data/portal_block_50.RDS"))

portal_time_series <- block %>%
  plot_time_series(base_size = 20) +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2015-01-01")),
    breaks = seq(as.Date("1980-01-01"), as.Date("2015-01-01"), by = "5 years"),
    date_labels = "%Y", expand = c(0.01, 0)
  ) +
  scale_color_viridis_d(option = "plasma") +
  guides(color = FALSE) +
  theme(
    panel.background = element_rect(fill = "#AAAABB", color = NA),
    panel.grid.major = element_line(color = "grey30", size = 1),
    panel.grid.minor = element_line(color = "grey30", size = 1),
    legend.key = element_rect(fill = "#AAAABB")
  )

ggsave(here::here("figures/portal_time_series.pdf"),
  portal_time_series,
  width = 12, height = 6
)

portal_time_series_regimes <- portal_time_series %>%
  add_regime_shift_highlight(fill = "black")
ggsave(here::here("figures/portal_time_series_regimes.pdf"),
  portal_time_series_regimes,
  width = 12, height = 6
)

results_file <- here::here("output/portal_ds_results_50.RDS")
results <- readRDS(results_file)
dom_ev_plot <- results$eigenvalues %>%
  plot_eigenvalues() +
  scale_y_continuous(limits = c(0.92, 1.03)) +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2015-01-01")),
    breaks = seq(as.Date("1980-01-01"), as.Date("2015-01-01"), by = "5 years"),
    date_labels = "%Y", expand = c(0.01, 0)
  )
dom_ev_plot <- add_regime_shift_highlight(dom_ev_plot)

ggsave(here::here("figures/portal_eigenvalues_50_dom.pdf"),
  dom_ev_plot,
  width = 12, height = 6
)

blank_ev_plot <- results$eigenvalues %>%
  plot_eigenvalues() +
  scale_y_continuous(limits = c(0.92, 1.03)) +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2015-01-01")),
    breaks = seq(as.Date("1980-01-01"), as.Date("2015-01-01"), by = "5 years"),
    date_labels = "%Y", expand = c(0.01, 0)
  ) + scale_color_manual(values = NA)
ggsave(here::here("figures/portal_eigenvalues_50_blank.pdf"),
  blank_ev_plot,
  width = 12, height = 6
)

eigenvectors <- results$eigenvectors

# empty plot
empty_plot <- results$eigenvectors %>%
  plot_eigenvectors(line_size = NA) +
  labs(x = "", y = "eigenvector composition", color = "Species") +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2015-01-01")),
    breaks = seq(as.Date("1980-01-01"), as.Date("2015-01-01"), by = "5 years"),
    date_labels = "%Y", expand = c(0.01, 0)
  )
ggsave(here::here("figures/portal_eigenvector_null.pdf"),
  empty_plot,
  width = 12, height = 6
)

# base plot of eigenvectors
base_plot <- results$eigenvectors %>%
  plot_eigenvectors() +
  scale_x_date(
    limits = as.Date(c("1980-01-01", "2015-01-01")),
    breaks = seq(as.Date("1980-01-01"), as.Date("2015-01-01"), by = "5 years"),
    date_labels = "%Y", expand = c(0.01, 0)
  ) +
  labs(x = "", y = "eigenvector composition", color = "Species")

ggsave(here::here("figures/portal_eigenvector_figure.pdf"),
  base_plot,
  width = 12, height = 6
)

ggsave(here::here("figures/portal_eigenvector_figure_regimes.pdf"),
  add_regime_shift_highlight(base_plot, fill = "black"),
  width = 12, height = 6
)


ggsave(here::here("figures/portal_eigenvector_figure_wide.pdf"),
  add_regime_shift_highlight(base_plot, fill = "black"),
  width = 24, height = 6
)
