extract_matrix_values <- function(values_list, id_var = "censusdate")
{
    purrr::map_dfr(values_list, .id = id_var, 
                   function(vals) {
                       if (any(is.na(vals))) {
                           return(data.frame())
                       }
                       data.frame(
                           value = vals, 
                           rank = seq(-vals)
                       )
                   }) %>%
        dplyr::mutate_at(vars(id_var), as.Date)
}

plot_matrix_values <- function(values_dist, id_var = "censusdate", 
                               y_label = "dynamic stability \n(higher is more unstable)", 
                               line_size = 1, base_size = 16)
{
    ggplot(values_dist, 
           aes(
               x = !!sym(id_var), y = .data$value,
               color = as.factor(.data$rank), group = rev(.data$rank)
           )) +
        geom_line(size = line_size) +
        scale_color_viridis_d(option = "inferno") +
        geom_hline(yintercept = 1.0, size = 1, linetype = 2) +
        labs(x = NULL, y = y_label, color = "rank") +
        theme_bw(
            base_size = base_size, base_family = "Helvetica",
            base_line_size = 1
        ) +
        theme(panel.grid.minor = element_line(size = 1)) +
        guides(color = FALSE)
}