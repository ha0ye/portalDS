extract_matrix_values <- function(values_list, id_var = "censusdate")
{
    purrr::map_dfr(values_list, 
                   .id = id_var, 
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

# normalize vectors so that length = 1
vector_scale <- function(v) {
    sum_sq <- sum(abs(v)^2)
    v / sqrt(sum_sq)
}

extract_matrix_vectors <- function(vectors_list, id_var = "censusdate", 
                                   rescale = TRUE, 
                                   row_idx = NULL, 
                                   col_idx = NULL)
{
    non_null_idx <- dplyr::first(which(!vapply(vectors_list, is.null, FALSE)))
    if (is.null(row_idx))
    {
        row_idx <- seq_len(NROW(vectors_list[[non_null_idx]]))
    }
    if (is.null(col_idx))
    {
        col_idx <- seq_len(NCOL(vectors_list[[non_null_idx]]))
    }
    
    out <- purrr::map_dfr(eigenvectors, 
                           .id = id_var, 
                           function(v) {
                               if (anyNA(v) || is.null(v)) {
                                   return()
                               }
                               reshape2::melt(v[row_idx, col_idx, drop = FALSE], as.is = TRUE)
                           }) %>% 
        dplyr::mutate_at(vars(id_var), as.Date) %>%
        dplyr::rename(variable = .data$Var1, rank = .data$Var2) %>%
        dplyr::mutate(value = abs(Re(.data$value)))
    
    if (rescale)
    {
        out <- out %>% 
            dplyr::group_by_at(c(id_var, "rank")) %>%
            dplyr::mutate(value = vector_scale(.data$value)) %>%
            dplyr::ungroup()
    }
    
    return(out)
}