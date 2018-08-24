library(ggraph)
library(igraph)
library(viridis)
library(cowplot)

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
        filter(lib_column != target_column) %>% 
        arrange(target_column) %>%
        select(target_column, lib_column) %>%
        graph_from_data_frame(vertices = levels(.$target_column))
    
    if (is.null(palette))
    {
        palette <- viridis(length(V(my_graph)), option = palette_option)
        names(palette) <- as_ids(V(my_graph))
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
        group_by(target, predictor) %>%
        summarize(v = max(abs(value))) %>%
        filter(v > 0) %>%
        mutate(coeff_name = paste0(target, predictor))
    smap_coeff_df <- smap_coeff_df %>% 
        mutate(coeff_name = paste0(target, predictor)) %>%
        filter(coeff_name %in% to_keep$coeff_name)
    
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

plot_eigenvalues <- function(results_file = here::here("output/portal_ds_results.RDS"), 
                             num_values = 1, plot_file = NULL, highlight_complex = FALSE, 
                             highlight_shifts = TRUE, width = 8, height = 4.5)
{
    results <- readRDS(results_file)
    eigenvalues <- results$eigenvalues
    
    # generate df for plotting
    eigenvalue_dist <- map_dfr(seq(eigenvalues), function(i) {
        lambda <- eigenvalues[[i]]
        if (any(is.na(lambda)))
            return(data.frame())
        lambda <- sort(abs(lambda), decreasing = TRUE)
        data.frame(lambda = lambda, censusdate = names(eigenvalues)[i], rank = seq(lambda))
    })
    eigenvalue_dist$censusdate <- as.Date(eigenvalue_dist$censusdate)
    
    to_plot <- eigenvalue_dist %>%
        filter(rank <= num_values)
    
    # construct plot
    ds_plot <- to_plot %>% 
        ggplot(aes(x = censusdate, y = lambda, color = as.factor(rank), group = rev(rank))) + 
        scale_color_manual(values = viridis(7, option = "A")[c(1, 3, 4, 5, 6)]) +
        geom_line(size = 1) + 
        geom_hline(yintercept = 1.0, size = 1, linetype = 2) + 
        scale_x_date(breaks = seq(from = as.Date("1985-01-01"), 
                                  to = as.Date("2015-01-01"), 
                                  by = "5 years"), 
                     date_labels = "%Y", expand = c(0.01, 0)) + 
        #scale_y_continuous(limits = c(0.90, 1.04)) + 
        labs(x = NULL, y = "dynamic stability \n(higher is more unstable)", color = "rank") +
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.grid.minor = element_line(size = 1)) + 
        guides(color = FALSE)
    
    if (highlight_complex)
    {
        complex_df <- data.frame(censusdate = to_plot %>% 
                                     spread(rank, lambda) %>%
                                     filter(`1` < `2` + 0.001) %>%
                                     select(censusdate), 
                                 lambda = min(to_plot$lambda, na.rm = TRUE), 
                                 rank = 1) %>%
            complete(censusdate = to_plot$censusdate, fill = list(lambda = NA, rank = 1))
        ds_plot <- ds_plot + 
            geom_point(data = complex_df, color = "red")
    }
    
    if (highlight_shifts)
    {
        ds_plot <- add_regime_shift_highlight(ds_plot)
    }
    
    # save output
    if (!is.null(plot_file)) 
    {
        ggsave(plot_file, ds_plot, 
               width = width, height = height)
    }
    return(ds_plot)
}

plot_eigenvectors <- function(results_file = here::here("output/portal_ds_results.RDS"), 
                              num_values = 1, plot_file = NULL, 
                              highlight_shifts = TRUE, width = 8, height = NULL)
{
    results <- readRDS(results_file)
    eigenvectors <- results$eigenvectors
    
    # extract vars
    non_null_idx <- first(which(!sapply(eigenvectors, is.null)))
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
    v_df <- map_dfr(seq(eigenvectors), function(i) {
        v <- eigenvectors[[i]]
        if (is.null(v))
            return()
        out <- reshape2::melt(v[var_idx, seq(num_values)])
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

    # compute IPR = Inverse Participation Ratio
    #   for each eigenvector
    #     normalize so that sum([v_i]^2) = 1
    #     IPR = sum([v_i]^4)
    #     ranges from 1/N (N = length of eigenvector) to 1
    ipr_df <- v_df %>% 
        group_by(t, rank) %>%
        summarize(value = sum(abs(value)^4)) %>%
        ungroup() %>%
        mutate(censusdate = as.Date(names(eigenvectors)[t]), 
               variable = "IPR")
    
    v_df$component <- "eigenvector"
    ipr_df$component <- "IPR"
    
    ev_plot <- ggplot(v_df, 
           aes(x = censusdate, y = value, color = variable)) + 
        facet_grid(component + rank ~ ., scales = "free", switch = "y") + 
        scale_color_viridis(discrete = TRUE, option = "plasma") + 
        scale_x_date(breaks = seq(from = as.Date("1985-01-01"), 
                                  to = as.Date("2015-01-01"), 
                                  by = "5 years"), 
                     date_labels = "%Y", expand = c(0.01, 0)) + 
        geom_line() + 
        theme_bw()
    ipr_plot <- ggplot(ipr_df, 
           aes(x = censusdate, y = value)) + 
        facet_grid(component + rank ~ ., switch = "y") + 
        scale_y_continuous(limits = c(1/length(var_idx), 1)) + 
        scale_x_date(breaks = seq(from = as.Date("1985-01-01"), 
                                  to = as.Date("2015-01-01"), 
                                  by = "5 years"), 
                     date_labels = "%Y", expand = c(0.01, 0)) + 
        geom_line(color = "black") + 
        theme_bw()
    
    if (highlight_shifts)
    {
        ev_plot <- add_regime_shift_highlight(ev_plot)
        ipr_plot <- add_regime_shift_highlight(ipr_plot)
    }
    
    if (is.null(height))
    {
        height <- num_values * 3
    }
    
    my_plot <- ev_plot
    # my_plot <- cowplot::plot_grid(ev_plot, ipr_plot, 
    #                               align = "v", axis = "lr", 
    #                               ncol = 1, labels = NA)
    
    # save output
    if (!is.null(plot_file))
    {
        ggsave(plot_file, my_plot, 
               width = width, height = height)
    }
    return(my_plot)
}

plot_time_series <- function(block)
{
    block %>%
        mutate_at(vars(-censusdate), 
                  function(x) {(x - min(x, na.rm = TRUE)) / 
                          (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}) %>%
        gather(species, abundance, -censusdate) %>%
        ggplot(aes(x = censusdate, y = abundance, color = species)) + 
        geom_line(size = 1) + 
        labs(x = NULL, y = "relative abundance") + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.grid.minor = element_line(size = 1))
}

add_regime_shift_highlight <- function(my_plot, alpha = 0.5, fill = "grey30")
{
    ## using dates from Christensen et al. 2018
    # lower_date <- as.Date(c("1983-12-01", "1988-10-01", "1998-09-01", "2009-06-01"))
    # upper_date <- as.Date(c("1984-07-01", "1996-01-01", "1999-12-01", "2010-09-01"))
    
    ## using dates from updated analysis code (weecology/LDA-kratplots)
    lower_date <- as.Date(c("1983-11-12", "1990-01-06", "1998-12-22", "2009-05-23"))
    upper_date <- as.Date(c("1985-03-16", "1992-04-04", "1999-11-06", "2011-01-05"))
    
    my_plot + geom_rect(data = data.frame(xmin = lower_date, xmax = upper_date,
                                          ymin = -Inf, ymax = Inf),
                        mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                        alpha = alpha, inherit.aes = FALSE, fill = fill)
}
