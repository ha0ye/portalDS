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

make_matrix_value_plot <- function(values_dist, id_var = "censusdate", 
                                   y_label = "dynamic stability \n(higher is more unstable)", 
                                   line_size = 1, base_size = 16)
{
    ggplot(values_dist, 
           aes(x = !!sym(id_var), y = .data$value,
               color = as.factor(.data$rank), group = rev(.data$rank)
           )) +
        geom_line(size = line_size) +
        scale_color_viridis_d(option = "inferno") +
        geom_hline(yintercept = 1.0, size = 1, linetype = 2) +
        labs(x = NULL, y = y_label, color = "rank") +
        my_theme(base_size = base_size) +
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
    non_null_idx <- dplyr::first(which(!vapply(vectors_list, anyNA, FALSE)))
    if (is.null(row_idx))
    {
        row_idx <- seq_len(NROW(vectors_list[[non_null_idx]]))
    }
    if (is.null(col_idx))
    {
        col_idx <- seq_len(NCOL(vectors_list[[non_null_idx]]))
    }
    
    out <- purrr::map_dfr(vectors_list, 
                          .id = id_var, 
                          function(v) {
                              if (anyNA(v) || is.null(v)) {
                                  return()
                              }
                              reshape2::melt(v[row_idx, col_idx, drop = FALSE], as.is = TRUE)
                          }) %>% 
        dplyr::mutate_at(vars(id_var), as.Date) %>%
        dplyr::rename(variable = .data$Var1, rank = .data$Var2) %>%
        dplyr::mutate(value = Re(.data$value))
    
    if (rescale)
    {
        out <- out %>% 
            dplyr::group_by_at(c(id_var, "rank")) %>%
            dplyr::mutate(value = vector_scale(.data$value)) %>%
            dplyr::ungroup()
    }
    
    return(out)
}

# compute IPR = Inverse Participation Ratio
#   for each eigenvector
#     normalize so that sum([v_i]^2) = 1
#     IPR = sum([v_i]^4)
#     ranges from 1/N (N = length of eigenvector) to 1
add_IPR <- function(v_df, comp_name = "eigenvector", id_var = "censusdate")
{
    ipr_df <- v_df %>%
        dplyr::group_by_at(c(id_var, "rank")) %>%
        dplyr::summarize(value = sum(abs(.data$value)^4)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(variable = "IPR")
    
    v_df$component <- comp_name
    ipr_df$component <- "IPR"
    
    dat <- dplyr::bind_rows(v_df, ipr_df)
    dat$variable <- as.factor(dat$variable)
    dat$variable <- forcats::fct_relevel(dat$variable, "IPR", after = Inf)
    return(dat)
}

make_matrix_vector_plot <- function(v_df, 
                                comp_name = "eigenvector", 
                                num_values = 1, 
                                id_var = "censusdate", 
                                add_IPR = FALSE, 
                                palette_option = "plasma", 
                                line_size = 1, base_size = 16)
{
    if (add_IPR)
    {
        v_df <-  add_IPR(v_df, comp_name = comp_name, id_var = id_var)
        row_facet_groups <- rlang::quos(component, rank)
    } else {
        row_facet_groups <- rlang::quos(rank)
    }
    
    my_plot <- v_df %>%
        ggplot(aes(x = !!sym(id_var), y = .data$value, color = .data$variable)) +
        facet_grid(rows = vars(!!!row_facet_groups), scales = "free", switch = "y") +
        scale_x_date(expand = c(0.01, 0)) +
        scale_y_continuous(limits = c(-1, 1)) +
        scale_color_viridis_d(option = palette_option) + 
        geom_line(size = line_size) +
        labs(x = "censusdate", y = "value", color = "variable") +
        my_theme(base_size = base_size, 
                 panel.background = element_rect(fill = "#AAAABB", color = NA),
                 panel.grid.major = element_line(color = "grey30", size = 1),
                 panel.grid.minor = element_line(color = "grey30", size = 1),
                 legend.key = element_rect(fill = "#AAAABB")
        ) +
        guides(color = guide_legend(override.aes = list(size = 1)))
    
    if (num_values == 1) {
        my_plot <- my_plot + theme(
            strip.background = element_blank(),
            strip.text.y = element_blank()
        )
    }
    return(my_plot)
}

#' @title add_regime_shift_highlight
#' @description add transparent bars to highlight specific time spans
#' @param my_plot the original ggplot object
#' @param lower_date a vector of the beginnings of the time spans
#' @param upper_date a vector of the ends of the time spans
#' @param alpha the transparency of the bars to add to the plot
#' @param fill the fill color of the bars to add to the plot
#'
#' @return A ggplot object with the bars added
#'
#' @export
add_regime_shift_highlight <- function(my_plot,
                                       ## using dates from updated analysis code (weecology/LDA-kratplots)
                                       lower_date = as.Date(c("1983-11-12", "1990-01-06", "1998-12-22", "2009-05-23")),
                                       upper_date = as.Date(c("1985-03-16", "1992-04-04", "1999-11-06", "2011-01-05")),
                                       alpha = 0.5, fill = "grey30") {
    ## using dates from Christensen et al. 2018
    # lower_date <- as.Date(c("1983-12-01", "1988-10-01", "1998-09-01", "2009-06-01"))
    # upper_date <- as.Date(c("1984-07-01", "1996-01-01", "1999-12-01", "2010-09-01"))
    
    my_plot + geom_rect(
        data = data.frame(
            xmin = lower_date, xmax = upper_date,
            ymin = -Inf, ymax = Inf
        ),
        mapping = aes(
            xmin = .data$xmin, xmax = .data$xmax,
            ymin = .data$ymin, ymax = .data$ymax
        ),
        alpha = alpha, inherit.aes = FALSE, fill = fill
    )
}

my_theme <- function(base_size, ...)
{
    theme_bw(base_size = base_size, base_family = "Helvetica",
             base_line_size = 1) +
        theme(panel.grid.minor = element_line(size = 1)) + 
        theme(...)
}