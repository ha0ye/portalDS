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
    
    if (is.null(plot_file))
    {
        ggsave(plot_file, combined_network_plot, width = 6, height = 12)
    } else {
        print(combined_network_plot)
    }
    return()
}

#### generate network figure from ccm_results ----

make_network_from_ccm_results <- function(ccm_results, 
                                          layout = "circle", 
                                          palette = NULL, 
                                          palette_option = "plasma", 
                                          existing_graph = NULL)
{
    my_graph <- ccm_results %>% 
        identify_ccm_links() %>%
        arrange(target_column) %>%
        select(target_column, lib_column) %>%
        graph_from_data_frame()
    
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
        geom_edge_fan(edge_width = 0.5) + 
        geom_node_circle(aes(r = 0.05, fill = name)) + 
        theme_graph(foreground = "black", fg_text_colour = "white", 
                    background = "transparent") + 
        coord_fixed() + 
        theme(text = element_text(family = "Helvetica"),
              panel.border = element_rect(color = NA, fill = NA)) + 
        scale_fill_manual(values = palette)
    
    return(list(plot = my_plot, 
                palette = palette, 
                graph = my_graph))
}

#### function to produce smap coeffs plot ----

plot_smap_coeffs <- function(smap_matrices, plot_file = NULL, 
                             width = 6, height = NULL)
{
    # make data.frame of smap coefficients
    smap_coeff_df <- map_df(seq(smap_matrices), function(i) {
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
    if (is.null(plot_file)) 
    {
        print(combined_plot)
    } else {
        ggsave(plot_file, combined_plot, 
               width = width, height = height)
    }
    return()
}

#### function to produce eigenvalues plot ----

plot_eigenvalues <- function(eigenvalues, num_values = 1, plot_file = NULL, 
                             width = 8, height = 4.5, highlight_shifts = FALSE)
{
    # generate df for plotting
    eigenvalue_dist <- map_df(seq(eigenvalues), function(i) {
        lambda <- eigenvalues[[i]]
        if(any(is.na(lambda)))
            return(data.frame())
        lambda <- sort(abs(lambda), decreasing = TRUE)
        data.frame(lambda = lambda, censusdate = names(eigenvalues)[i], rank = seq(lambda))
    })
    eigenvalue_dist$censusdate <- as.Date(eigenvalue_dist$censusdate)
    
    # construct plot
    portal_ds_plot <- eigenvalue_dist %>%
        filter(rank <= num_values) %>% 
        ggplot(aes(x = censusdate, y = lambda, color = as.factor(rank))) + 
        scale_color_viridis_d() +
        geom_line() + 
        geom_hline(yintercept = 1.0, size = 1, linetype = 2) + 
        scale_x_date(breaks = seq(from = as.Date("1985-01-01"), 
                                  to = as.Date("2015-01-01"), 
                                  by = "5 years"), 
                     date_labels = "%Y", expand = c(0.01, 0)) + 
        labs(x = NULL, y = "dynamic stability \n(higher is more unstable)", color = "rank") +
        theme_bw() + 
        theme(panel.grid.minor = element_line(size = 0.5))
    
    if (highlight_shifts)
    {
        portal_ds_plot <- portal_ds_plot + 
            geom_rect(data = data.frame(xmin = as.Date("1983-12-01"), xmax = as.Date("1984-07-01"), 
                                        ymin = -Inf, ymax = Inf), 
                      mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                      alpha = 0.3, inherit.aes = FALSE) + 
            #    geom_rect(data = data.frame(xmin = as.Date("1988-10-01"), xmax = as.Date("1996-01-01"), 
            #                                ymin = -Inf, ymax = Inf), 
            #              mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            #              alpha = 0.3, inherit.aes = FALSE) + 
            geom_rect(data = data.frame(xmin = as.Date("1998-09-01"), xmax = as.Date("1999-12-01"), 
                                        ymin = -Inf, ymax = Inf), 
                      mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                      alpha = 0.3, inherit.aes = FALSE) + 
            geom_rect(data = data.frame(xmin = as.Date("2009-06-01"), xmax = as.Date("2010-09-01"), 
                                        ymin = -Inf, ymax = Inf), 
                      mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                      alpha = 0.3, inherit.aes = FALSE)
    }
    
    # save output
    if (is.null(plot_file)) 
    {
        print(portal_ds_plot)
    } else {
        ggsave(plot_file, portal_ds_plot, 
               width = width, height = height)
    }
    return()
}