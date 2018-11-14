#' @importFrom magrittr "%>%"
#' @import ggraph
#' @import ggplot2

#### make networks from all ccm results, match color and location to first ----

make_combined_network <- function(plot_file = NULL)
{
    portal_network <- here::here("output/portal_ds_results.RDS") %>%
        readRDS() %>%
        .$ccm_results %>%
        make_network_from_ccm_results()

    portal_network_50 <- here::here("output/portal_ds_results_50.RDS") %>%
        readRDS() %>%
        .$ccm_results %>%
        make_network_from_ccm_results(palette = portal_network$palette,
                                      existing_graph = portal_network$graph)

    portal_network_33 <- here::here("output/portal_ds_results_33.RDS") %>%
        readRDS() %>%
        .$ccm_results %>%
        make_network_from_ccm_results(palette = portal_network$palette,
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

#### generate network figure from ccm_results ----

make_network_from_ccm_results <- function(results_file = here::here("output/portal_ds_results.RDS"),
                                          layout = "circle",
                                          palette = NULL,
                                          palette_option = "plasma",
                                          existing_graph = NULL)
{
    results <- readRDS(results_file)

    my_graph <- results$ccm_links %>%
        dplyr::filter(lib_column != target_column) %>%
        dplyr::arrange(target_column) %>%
        dplyr::select(target_column, lib_column) %>%
        igraph::graph_from_data_frame(vertices = levels(.$target_column))

    if (is.null(palette))
    {
        palette <- viridis::viridis(length(V(my_graph)), option = palette_option)
        names(palette) <- igraph::as_ids(V(my_graph))
    }

    my_graph <- create_layout(my_graph, layout = layout)

    if (!is.null(existing_graph))
    {
        idx <- match(my_graph$name, existing_graph$name)
        my_graph$x <- existing_graph$x[idx]
        my_graph$y <- existing_graph$y[idx]
    }

    my_plot <- ggraph(my_graph) +
        geom_edge_link(edge_width = 0.5, start_cap = circle(0.3, "inches"),
                       end_cap = circle(0.3, "inches"),
                       arrow = arrow(angle = 20, type = "closed")) +
        geom_node_circle(aes(r = 0.08, fill = name)) +
        theme_graph(foreground = "black", fg_text_colour = "white",
                    background = "transparent") +
        coord_fixed() +
        theme(text = element_text(family = "Helvetica"),
              panel.border = element_rect(color = NA, fill = NA)) +
        scale_fill_manual(values = palette) +
        guides(fill = guide_legend(title = "Species"))

    return(list(plot = my_plot,
                palette = palette,
                graph = my_graph))
}

#### function to produce smap coeffs plot ----

plot_smap_coeffs <- function(results_file = here::here("output/portal_ds_results.RDS"),
                             plot_file = NULL, width = 6, height = NULL)
{
    results <- readRDS(results_file)
    smap_matrices <- results$smap_matrices

    # make data.frame of smap coefficients
    smap_coeff_df <- map_dfr(seq(smap_matrices), function(i) {
        m <- smap_matrices[[i]]
        if (is.null(dim(m)))
            return()
        row_idx <- grep("_0", rownames(m))
        out <- reshape2::melt(m[row_idx,])
        out$t <- i
        return(out)
    }) %>%
        rename(target = Var1, predictor = Var2)

    # identify coefficients that matter
    to_keep <- smap_coeff_df %>%
        dplyr::group_by(target, predictor) %>%
        dplyr::summarize(v = max(abs(value))) %>%
        dplyr::filter(v > 0) %>%
        dplyr::mutate(coeff_name = paste0(target, predictor))
    smap_coeff_df <- smap_coeff_df %>%
        dplyr::mutate(coeff_name = paste0(target, predictor)) %>%
        dplyr::filter(coeff_name %in% to_keep$coeff_name)

    # convert time index into dates
    smap_coeff_df$censusdate <- as.Date(names(smap_matrices)[smap_coeff_df$t])

    # time series plot
    ts_plot <- ggplot(smap_coeff_df,
                      aes(x = censusdate, y = abs(value), color = predictor)) +
        facet_grid(target ~ ., scales = "free_y", switch = "y") +
        geom_hline(yintercept = 1, size = 1, linetype = 2) +
        scale_color_viridis_d(option = "E") +
        scale_x_date(breaks = seq(from = as.Date("1985-01-01"),
                                  to = as.Date("2015-01-01"),
                                  by = "5 years"),
                     date_labels = "%Y", expand = c(0.01, 0)) +
        geom_line() +
        theme_bw() +
        guides(color = FALSE, fill = FALSE)

    # density plot
    density_plot <- ggplot(smap_coeff_df,
                           aes(x = abs(value), color = predictor)) +
        facet_grid(target ~ ., switch = "y") +
        geom_vline(xintercept = 1, size = 1, linetype = 2) +
        scale_color_viridis_d(option = "E") +
        geom_density(fill = NA, weight = 0.5) +
        coord_flip() +
        theme_bw() +
        guides(color = FALSE, fill = FALSE)

    combined_plot <- plot_grid(ts_plot, density_plot, nrow = 1,
                                        rel_widths = c(3, 1))
    if (is.null(height))
    {
        height <- nlevels(smap_coeff_df$target)
    }

    # save output
    if (!is.null(plot_file))
    {
         ggsave(plot_file, combined_plot,
               width = width, height = height)
    }
    return(combined_plot)
}

#### function to produce eigenvalues plot ----

#' @export
plot_eigenvalues <- function(eigenvalues, num_values = 1, highlight_complex = FALSE, line_size = 1)
{
    # generate df for plotting
    eigenvalue_dist <- purrr::map_dfr(seq(eigenvalues), function(i) {
        lambda <- eigenvalues[[i]]
        if (any(is.na(lambda)))
            return(data.frame())
        lambda <- sort(abs(lambda), decreasing = TRUE)
        data.frame(lambda = lambda, censusdate = names(eigenvalues)[i], rank = seq(lambda))
    }) %>%
        dplyr::filter(rank <= num_values) %>%
        dplyr::mutate(censusdate = as.Date(censusdate))

    ds_plot <- eigenvalue_dist %>%
        ggplot(aes(x = censusdate, y = lambda, color = as.factor(rank), group = rev(rank))) +
        geom_line(size = line_size) +
        scale_color_manual(values = viridis(7, option = "magma")[c(1, 3, 5, 6, 7)]) +
        geom_hline(yintercept = 1.0, size = 1, linetype = 2) +
        labs(x = NULL, y = "dynamic stability \n(higher is more unstable)", color = "rank") +
        theme_bw(base_size = 20, base_family = "Helvetica",
                 base_line_size = 1) +
        theme(panel.grid.minor = element_line(size = 1)) +
        guides(color = FALSE)

    if (highlight_complex)
    {
        complex_df <- data.frame(censusdate = eigenvalue_dist %>%
                                     dplyr::spread(rank, lambda) %>%
                                     dplyr::filter(`1` < `2` + 0.001) %>%
                                     dplyr::select(censusdate),
                                 lambda = min(eigenvalue_dist$lambda, na.rm = TRUE),
                                 rank = 1) %>%
            tidyr::complete(censusdate = eigenvalue_dist$censusdate, fill = list(lambda = NA, rank = 1))
        ds_plot <- ds_plot +
            geom_point(data = complex_df, color = "red")
    }
    return(ds_plot)
}

#' @export
plot_eigenvectors <- function(eigenvectors, num_values = 1, add_IPR = FALSE, line_size = 1)
{
    # extract vars
    non_null_idx <- dplyr::first(which(!sapply(eigenvectors, is.null)))
    var_names <- rownames(eigenvectors[[non_null_idx]])
    var_idx <- grep("_0", var_names)
    var_names <- gsub("_0", "", var_names[var_idx])

    # normalize eigenvectors so that length = 1
    vector_scale <- function(v)
    {
        sum_sq <- sum(abs(v)^2)
        v / sqrt(sum_sq)
    }

    # make data.frame of eigenvector components
    v_df <- purrr::map_dfr(seq(eigenvectors), function(i) {
        v <- eigenvectors[[i]]
        if (is.null(v))
            return()
        out <- reshape2::melt(v[var_idx, seq(num_values), drop = FALSE])
        out$t <- i
        return(out)
    }) %>%
        dplyr::rename(variable = Var1, rank = Var2) %>%
        dplyr::mutate(censusdate = as.Date(names(eigenvectors)[t]),
               variable = as.factor(var_names[variable]),
               value = abs(Re(value))) %>%
        dplyr::group_by(t, rank) %>%
        dplyr::mutate(value = vector_scale(value)) %>%
        dplyr::ungroup()

    if (add_IPR)
    {
        # compute IPR = Inverse Participation Ratio
        #   for each eigenvector
        #     normalize so that sum([v_i]^2) = 1
        #     IPR = sum([v_i]^4)
        #     ranges from 1/N (N = length of eigenvector) to 1
        ipr_df <- v_df %>%
            dplyr::group_by(t, rank) %>%
            dplyr::summarize(value = sum(abs(value)^4)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(censusdate = as.Date(names(eigenvectors)[t]),
                   variable = "IPR")

        v_df$component <- "eigenvector"
        ipr_df$component <- "IPR"

        dat <- dplyr::bind_rows(v_df, ipr_df)
        dat$variable <- as.factor(dat$variable)
        dat$variable <- forcats::fct_relevel(dat$variable, c(var_names, "IPR"))
        palette <- c(viridis::viridis(length(var_names), option = "plasma"), "black")

        my_plot <- dat %>%
            ggplot(aes(x = censusdate, y = value, color = variable)) +
            facet_grid(component + rank ~ ., switch = "y") +
            scale_y_continuous(limits = c(0, 1)) +
            scale_color_manual(values = palette)
    } else {
        my_plot <- v_df %>%
            ggplot(aes(x = censusdate, y = value, color = variable)) +
            facet_grid(rank ~ ., scales = "free", switch = "y") +
            scale_color_viridis(discrete = TRUE, option = "plasma")
    }
    my_plot <- my_plot +
        scale_x_date(expand = c(0.01, 0)) +
        geom_line(size = line_size) +
        theme_bw(base_size = 20, base_family = "Helvetica",
                 base_line_size = 1) +
        theme(panel.background = element_rect(fill = "#AAAABB", color = NA),
              panel.grid.major = element_line(color = "grey30", size = 1),
              panel.grid.minor = element_line(color = "grey30", size = 1),
              legend.key = element_rect(fill = "#AAAABB")) +
        guides(color = guide_legend(override.aes = list(size = 1)))

    if (num_values == 1)
    {
        my_plot <- my_plot + theme(strip.background = element_blank(),
                                   strip.text.y = element_blank())
    }

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
#'
#' @return A ggplot object of the time series plot
#'
#' @export

plot_time_series <- function(block,
                             time_column = censusdate,
                             scale = "unif",
                             theme = "print",
                             base_size = 11,
                             base_family = "Helvetica",
                             base_line_size = 1)
{
    time_column <- rlang::enquo(time_column)
    time <- dplyr::pull(block, !!time_column)
    block <- dplyr::select(block, -!!time_column)

    if (tolower(scale) == "unif")
    {
        block <- block %>%
            dplyr::mutate_all(scales::rescale)
    } else if (tolower(scale) == "norm") {
        block <- block %>%
            dplyr::mutate_all(norm_rescale)
    }
    block %>%
        dplyr::mutate(time = time) %>%
        tidyr::gather(species, abundance, -time) %>%
        ggplot(aes(x = time, y = abundance, color = species)) +
        geom_line(size = 1) +
        labs(x = NULL, y = "relative abundance") +
        theme_bw(base_size = base_size, base_family = base_family,
                 base_line_size = base_line_size) +
        theme(panel.grid.minor = element_line(size = 1))
}

#' @export
add_regime_shift_highlight <- function(my_plot,
                                       ## using dates from updated analysis code (weecology/LDA-kratplots)
                                       lower_date = as.Date(c("1983-11-12", "1990-01-06", "1998-12-22", "2009-05-23")),
                                       upper_date = as.Date(c("1985-03-16", "1992-04-04", "1999-11-06", "2011-01-05")),
                                       alpha = 0.5, fill = "grey30")
{
    ## using dates from Christensen et al. 2018
    # lower_date <- as.Date(c("1983-12-01", "1988-10-01", "1998-09-01", "2009-06-01"))
    # upper_date <- as.Date(c("1984-07-01", "1996-01-01", "1999-12-01", "2010-09-01"))

    my_plot + geom_rect(data = data.frame(xmin = lower_date, xmax = upper_date,
                                          ymin = -Inf, ymax = Inf),
                        mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                        alpha = alpha, inherit.aes = FALSE, fill = fill)
}
