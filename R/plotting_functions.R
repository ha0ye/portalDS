#' @importFrom rlang .data

#' @title plot_network
#' @description Visualize the network of interactions, created from running CCM
#'   on community time series as part of the dynamic stability analysis
#' @param ccm_links a data.frame containing the inferred interactions from CCM,
#'   it should have a `lib_column` and `target_column` to specify the causal
#'   links (where directed edges are from `target_column` to `lib_column`)
#' @param palette a data.frame with the colors for each vertex; if NULL, one
#'   will be created using [viridis::viridis()]
#' @param palette_option the color palette to use (see [viridis::viridis()] for
#'   more info)
#' @param existing_graph a data.frame specifying the layout of nodes (e.g. from a
#'   previous call to plot_network); if NULL, one will be created
#' @inheritParams ggraph::ggraph
#'
#' @return a list with three elements:
#' \tabular{ll}{
#'   \code{plot} \tab the ggraph object to plot or save\cr
#'   \code{palette} \tab the palette (for a future `plot_network` call)\cr
#'   \code{graph} \tab the graph (for a future `plot_network` call)\cr
#' }
#'
#' @export
plot_network <- function(ccm_links,
                         layout = "circle",
                         palette = NULL,
                         palette_option = "plasma",
                         existing_graph = NULL) {
    my_graph <- ccm_links %>%
        dplyr::filter(.data$lib_column != .data$target_column) %>%
        dplyr::arrange(.data$target_column) %>%
        dplyr::select(c("target_column", "lib_column")) %>%
        igraph::graph_from_data_frame(vertices = levels(ccm_links$target_column))
    
    if (is.null(palette)) {
        vertices <- igraph::V(my_graph)
        palette <- viridis::viridis(length(vertices), option = palette_option)
        names(palette) <- igraph::as_ids(vertices)
    }
    
    my_graph <- create_layout(my_graph, layout = layout)
    
    if (!is.null(existing_graph)) {
        idx <- match(my_graph$name, existing_graph$name)
        my_graph$x <- existing_graph$x[idx]
        my_graph$y <- existing_graph$y[idx]
    }
    
    my_plot <- ggraph(my_graph) +
        geom_edge_link(
            edge_width = 0.5, start_cap = circle(0.3, "inches"),
            end_cap = circle(0.3, "inches"),
            arrow = arrow(angle = 20, type = "closed")
        ) +
        geom_node_circle(aes(r = 0.08, fill = .data$name)) +
        theme_graph(
            foreground = "black", fg_text_colour = "white",
            background = "transparent"
        ) +
        coord_fixed() +
        theme(
            text = element_text(family = "Helvetica"),
            panel.border = element_rect(color = NA, fill = NA)
        ) +
        scale_fill_manual(values = palette) +
        guides(fill = guide_legend(title = "Species"))
    
    return(list(
        plot = my_plot,
        palette = palette,
        graph = my_graph
    ))
}

#' @title plot_smap_coeffs
#' @description Visualize the smap-coefficients from running the S-map model on
#'   the community time series as part of the dynamic stability analysis
#' @param smap_matrices a list of matrices for the s-map coefficients:
#'   the number of elements in the list corresponds to the time points of the
#'   s-map model, and each element is a matrix of the s-map coefficients at
#'   that time step
#' @param base_size the base font size
#' @param plot_file the filepath to where to save the plot; if `NULL` (the
#'   default), then the plot is not saved to a file
#' @param width width of the saved plot
#' @param height height of the saved plot
#'
#' @return A ggplot object of the smap coefficients plot
#'
#' @export
plot_smap_coeffs <- function(smap_matrices, base_size = 16,
                             plot_file = NULL, width = 6, height = NULL)
{
    # make data.frame of smap coefficients
    smap_coeff_df <- purrr::map_dfr(seq(smap_matrices), function(i) {
        m <- smap_matrices[[i]]
        if (is.null(dim(m))) {
            return()
        }
        row_idx <- grep("_0", rownames(m))
        out <- reshape2::melt(m[row_idx, ])
        out$t <- i
        return(out)
    }) %>%
        dplyr::rename(target = .data$Var1, predictor = .data$Var2)
    
    # identify coefficients that matter
    to_keep <- smap_coeff_df %>%
        dplyr::group_by(.data$target, .data$predictor) %>%
        dplyr::filter(max(abs(.data$value)) > 0) %>%
        dplyr::mutate(coeff_name = paste0(.data$target, .data$predictor))
    smap_coeff_df <- smap_coeff_df %>%
        dplyr::mutate(coeff_name = paste0(.data$target, .data$predictor)) %>%
        dplyr::filter(.data$coeff_name %in% to_keep$coeff_name)
    
    # convert time index into dates
    smap_coeff_df$censusdate <- as.Date(names(smap_matrices)[smap_coeff_df$t])
    
    # time series plot
    ts_plot <- ggplot(smap_coeff_df,
                      aes(x = .data$censusdate, y = abs(.data$value), 
                          color = .data$predictor)) +
        facet_grid(target ~ ., scales = "free_y", switch = "y") +
        geom_hline(yintercept = 1, size = 1, linetype = 2) +
        scale_color_viridis_d(option = "E") +
        scale_x_date(
            breaks = seq(
                from = as.Date("1985-01-01"),
                to = as.Date("2015-01-01"),
                by = "5 years"
            ),
            date_labels = "%Y", expand = c(0.01, 0)
        ) +
        geom_line() +
        labs(x = "censusdate", y = "abs(value)", color = "predictor") +
        my_theme(base_size = base_size, 
                 panel.grid.minor = element_line(size = 0.5)) + 
        guides(color = FALSE, fill = FALSE)
    
    # density plot
    density_plot <- ggplot(smap_coeff_df,
                           aes(x = abs(.data$value), color = .data$predictor)) +
        facet_grid(target ~ ., switch = "y") +
        geom_vline(xintercept = 1, size = 1, linetype = 2) +
        scale_color_viridis_d(option = "E") +
        geom_density(fill = NA, weight = 0.5) +
        coord_flip() +
        labs(x = "abs(value)", y = "density", color = "predictor") +
        my_theme(base_size = base_size, 
                 panel.grid.minor = element_line(size = 0.5)) + 
        guides(color = FALSE, fill = FALSE)
    
    combined_plot <- cowplot::plot_grid(ts_plot, density_plot,
                                        nrow = 1,
                                        rel_widths = c(3, 1)
    )
    if (is.null(height)) {
        height <- nlevels(smap_coeff_df$target)
    }
    
    if (!is.null(plot_file)) {  cowplot::ggsave(plot_file, combined_plot, width = width, height = height)  }
    return(combined_plot)
}

#' @title plot_eigenvalues
#' @aliases plot_svd_values
#' @description [plot_eigenvalues()] visualizes the dominant eigenvalue(s) from 
#'   running the S-map model on the community time series
#' @param eigenvalues a list of vectors for the eigenvalues:
#'   the number of elements in the list corresponds to the time points of the
#'   s-map model, and each element is a vector of the eigenvalues, computed
#'   from the matrix of the s-map coefficients at that time step
#' @param num_values the number of eigenvalues to plot
#' @param id_var when constructing the long-format tibble,what should be the 
#'   variable name containing the time index
#' @param highlight_complex whether to also draw points to indicate when the
#'   dominant eigenvalue is complex
#' @param line_size the line width for the plot
#' @inheritParams plot_smap_coeffs
#'
#' @return A ggplot object of the plot
#'
#' @export
plot_eigenvalues <- function(eigenvalues, num_values = 1,
                             id_var = "censusdate",
                             highlight_complex = FALSE, 
                             line_size = 1, base_size = 16,
                             plot_file = NULL, width = 6, height = NULL)
{
    eigenvalue_dist <- extract_matrix_values(eigenvalues, id_var = id_var) %>%
        dplyr::filter(.data$rank <= num_values) %>%
        dplyr::mutate(value = abs(.data$value))
    
    my_plot <- make_matrix_value_plot(eigenvalue_dist, 
                                      id_var = id_var, 
                                      y_label = "dynamic stability \n(higher is more unstable)", 
                                      line_size = line_size, 
                                      base_size = base_size)
    
    if (highlight_complex && num_values >= 2) {
        complex_df <- data.frame(
            censusdate = eigenvalue_dist %>%
                tidyr::spread(.data$rank, .data$value) %>%
                dplyr::filter(.data$`1` < .data$`2` + 0.001) %>%
                dplyr::select(.data$censusdate),
            value = min(eigenvalue_dist$value, na.rm = TRUE),
            rank = 1
        ) %>%
            tidyr::complete(censusdate = eigenvalue_dist$censusdate, fill = list(lambda = NA, rank = 1))
        my_plot <- my_plot +
            geom_point(data = complex_df, color = "red")
    }
    
    if (!is.null(plot_file)) {  cowplot::ggsave(plot_file, my_plot, width = width, height = height)  }
    return(my_plot)
}

#' @rdname plot_eigenvalues
#' @description [plot_svd_values()] visualizes the dominant singular value(s) 
#'   from running the S-map model on the community time series
#' @param singular_values a list of vectors for the singular values:
#'   the number of elements in the list corresponds to the time points of the
#'   s-map model, and each element is a vector of the singular values, computed
#'   from the matrix of the s-map coefficients at that time step
#' @inheritParams plot_eigenvalues
#'
#' @export
plot_svd_values <- function(singular_values, num_values = 1,
                            id_var = "censusdate",
                            line_size = 1,
                            base_size = 16,
                            plot_file = NULL, width = 6, height = NULL)
{
    sigma_dist <- extract_matrix_values(singular_values) %>%
        dplyr::filter(.data$rank <= num_values)
    
    my_plot <- make_matrix_value_plot(sigma_dist, 
                                      id_var = id_var, 
                                      y_label = "local convergence \n(higher is more divergence)", 
                                      line_size = line_size, 
                                      base_size = base_size)
    
    if (!is.null(plot_file)) {  cowplot::ggsave(plot_file, my_plot, width = width, height = height)  }
    return(my_plot)
}

#' @title Plot time-varying vector components
#' @aliases plot_svd_vectors
#' @description [plot_eigenvectors()] visualizes the dominant eigenvector(s) 
#'   from running the S-map model on the community time series
#' @param eigenvectors a list of matrices for the eigenvectors:
#'   the number of elements in the list corresponds to the time points of the
#'   s-map model, and each element is a matrix, where the columns are the
#'   eigenvectors, in descending order according to the eigenvalues
#' @param num_values the number of eigenvectors to plot
#' @param add_IPR whether to also plot the Inverse Participation Ratio, a
#'   numerical quantity that measures how evenly the different components
#'   contribute to the eigenvector
#' @inheritParams plot_eigenvalues
#' @inheritParams plot_network
#' @inheritParams plot_smap_coeffs
#'
#' @return A ggplot object of the plot
#'
#' @export
plot_eigenvectors <- function(eigenvectors, num_values = 1,
                              id_var = "censusdate", 
                              add_IPR = FALSE, 
                              palette_option = "plasma",
                              line_size = 1, base_size = 16,
                              plot_file = NULL, width = 6, height = NULL)
{
    non_null_idx <- dplyr::first(which(!vapply(eigenvectors, anyNA, FALSE)))
    var_names <- rownames(eigenvectors[[non_null_idx]])
    var_idx <- grep("_0", var_names)
    
    v_df <- extract_matrix_vectors(eigenvectors, 
                                   id_var = id_var, 
                                   rescale = TRUE, 
                                   row_idx = var_idx, 
                                   col_idx = seq_len(num_values)) %>%
        dplyr::mutate(variable = gsub("_0", "", .data$variable))
    
    my_plot <- make_matrix_vector_plot(v_df, 
                                       comp_name = "eigenvector", 
                                       num_values = num_values, 
                                       id_var = id_var, 
                                       add_IPR = add_IPR, 
                                       palette_option = palette_option, 
                                       line_size = line_size, 
                                       base_size = base_size)
    
    if (!is.null(plot_file)) {  cowplot::ggsave(plot_file, my_plot, width = width, height = height)  }
    return(my_plot)
}

#' @rdname plot_eigenvectors
#' @aliases plot_svd_vectors
#' @description [plot_svd_vectors()] visualizes the dominant SVD vector(s) 
#'   from running the S-map model on the community time series
#' @param svd_vectors a list of matrices for the SVD vectors:
#'   the number of elements in the list corresponds to the time points of the
#'   s-map model, and each element is a matrix, where the columns are the
#'   the SVD vectors, in descending order according to the singular values
#'
#' @export
plot_svd_vectors <- function(svd_vectors, num_values = 1,
                             id_var = "censusdate", 
                             add_IPR = FALSE, 
                             palette_option = "plasma",
                             line_size = 1, base_size = 16,
                             plot_file = NULL, width = 6, height = NULL)
{
    v_df <- extract_matrix_vectors(svd_vectors, 
                                   rescale = FALSE, 
                                   col_idx = seq_len(num_values))
    
    my_plot <- make_matrix_vector_plot(v_df, 
                                       comp_name = "svd vector", 
                                       num_values = num_values, 
                                       id_var = id_var, 
                                       add_IPR = add_IPR, 
                                       palette_option = palette_option, 
                                       line_size = line_size, 
                                       base_size = base_size)
    
    if (!is.null(plot_file)) {  cowplot::ggsave(plot_file, my_plot, width = width, height = height)  }
    return(my_plot)
}

#' @title plot_time_series
#' @description plot the time series in `block`, using the appropriate rescaling
#'   of the data
#' @param block A data.frame containing time series for the community. Each
#'   column is a time series of abundances.
#' @param time_column The name of the column in the block that has the time,
#'   which could be a numeric or a date/time type
#' @param scale How to scale the time series:
#'   `unif` -- scale each time series to be on [0, 1]
#'   `norm` -- scale each time series to have mean 0 and variance 1
#'   (anything else) -- no scaling
#' @inheritParams ggplot2::theme_bw
#' @inheritParams plot_network
#' @inheritParams plot_eigenvalues
#' 
#' @return A ggplot object of the time series plot
#'
#' @export
plot_time_series <- function(block,
                             time_column = "censusdate",
                             scale = "unif",
                             palette_option = "plasma",
                             line_size = 1, base_size = 11)
{
    time <- dplyr::pull(block, time_column)
    abundances <- dplyr::select(block, -time_column)
    
    # do scaling and setup y-axis label
    y_label <- "abundance"
    if (is.null(scale) || length(scale) == 0)
    {
    } else if (tolower(scale) == "unif") {
        abundances <- abundances %>%
            dplyr::mutate_all(scales::rescale)
        y_label <- "relative abundance"
    } else if (tolower(scale) == "norm") {
        abundances <- abundances %>%
            dplyr::mutate_all(norm_rescale)
        y_label <- "scaled abundance"
    }
    
    # convert to long format
    to_plot <- abundances %>%
        tibble::add_column(time = time) %>%
        tidyr::gather("species", "abundance", -.data$time)

    # make plot
    to_plot %>%
        ggplot(aes(x = .data$time, y = .data$abundance, 
                   color = .data$species)) +
        geom_line(size = line_size) +
        scale_color_viridis_d(option = palette_option) + 
        labs(x = NULL, y = y_label, color = "species") +
        my_theme(base_size = base_size)
}
