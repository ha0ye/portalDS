#### setup ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rEDM)
library(parallel)
# devtools::install_github("weecology/portalr")
# devtools::install_github("weecology/LDATS")

source(here::here("R", "dynamic_stability_functions.R"))

#### function to fit periodic spline through data, using yearday as input x ----

make_surrogate_annual_spline <- function(x, y, num_surr = 100)
{
    # filter out NA first
    n <- length(x)
    idx <- which(is.finite(x) & is.finite(y))
    x <- x[idx]
    y <- y[idx]
    
    xx <- c(x - 365, x, x + 365)
    yy <- c(y, y, y)
    seasonal_F <- smooth.spline(xx, yy)
    seasonal_cyc <- predict(seasonal_F, x)$y
    seasonal_resid <- y - seasonal_cyc
    
    vals <- matrix(unlist(
        lapply(seq(num_surr), function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        })
    ), ncol = num_surr)
    
    out <- matrix(NA, nrow = n, ncol = num_surr)
    out[idx, ] <- vals 
}

#### function for full dynamic stability analysis ----

compute_portal_dynamic_stability <- function(proc_data_file = "data/portal_block.RDS", 
                                             simplex_results_file = "output/portal_simplex_results.RDS", 
                                             ccm_results_file = "output/portal_ccm_results.RDS", 
                                             smap_matrices_file = "output/portal_smap_matrices.RDS", 
                                             eigenvalues_file = "output/portal_eigenvalues.RDS", 
                                             filter_q = NULL, 
                                             file_suffix = NULL, 
                                             num_cores = 2)
{

    if (is.null(file_suffix))
    {
        file_suffix <- ""
        if (!is.null(filter_q))
            file_suffix <- paste0("_", as.character(floor(filter_q * 100)))
    }
    proc_data_file <- sub(".RDS", paste0(file_suffix, ".RDS"), proc_data_file)
    simplex_results_file <- sub(".RDS", paste0(file_suffix, ".RDS"), simplex_results_file)
    ccm_results_file <- sub(".RDS", paste0(file_suffix, ".RDS"), ccm_results_file)
    smap_matrices_file <- sub(".RDS", paste0(file_suffix, ".RDS"), smap_matrices_file)
    eigenvalues_file <- sub(".RDS", paste0(file_suffix, ".RDS"), eigenvalues_file)
    
    # pre-processing
    if (!file.exists(proc_data_file))
    {
        # options here are:
        #   time = "all" (allow us to do correct interpolation and accouting for seasonality)
        #   level = "plot" (allow us to pull out abundances on the plots we want)
        #   effort = TRUE (so we can check effort)
        #   na_drop = TRUE (ignore periods where sampling did not occur)
        raw_rodent_data <- portalr::abundance(time = "all", 
                                              level = "plot", 
                                              effort = TRUE, 
                                              na_drop = TRUE)
        
        # summarize by each newmmonnumber, and for only the control plots we want
        block <- raw_rodent_data %>% 
            filter(plot %in% c(2, 4, 8, 11, 12, 14, 17, 22), 
                   censusdate < "2015-04-18") %>% 
            select(-treatment, -plot, -period) %>%
            group_by(newmoonnumber, censusdate) %>%
            summarize_all(sum) %>%
            ungroup()
        
        # check that effort is equal across samples
        stopifnot(length(unique(block$ntraps)) == 1)
        
        if (!is.null(filter_q))
        {
            species_list <- block %>%
                gather(species, abundance, BA:SO) %>%
                group_by(species) %>%
                summarize(quantile_q = quantile(abundance, 1 - filter_q)) %>%
                filter(quantile_q > 0) %>%
                pull(species)
            
            block <- block %>% 
                select(c("newmoonnumber", "censusdate", "ntraps", species_list))
        }
        
        # add in NAs for unsampled newmoonnumbers and interpolate
        block <- block %>%
            complete(newmoonnumber = full_seq(newmoonnumber, 1), fill = list(NA)) %>%
            mutate_at(vars(-newmoonnumber, -ntraps), forecast::na.interp) %>%
            mutate(censusdate = as.Date(as.numeric(censusdate), origin = "1970-01-01")) %>% 
            select(-newmoonnumber, -ntraps)
        
        saveRDS(block, proc_data_file)
    }
    block <- readRDS(proc_data_file)
    block_yearday <- yday(block$censusdate)
    block <- block %>% select(-censusdate) # drop the time column
    
    # identify best embedding dimension and generate surrogate data
    if (!file.exists(simplex_results_file))
    {
        simplex_results <- compute_simplex(block, E = 1:16)
        
        simplex_results$surrogate_data <- 
            map(seq(NROW(simplex_results)), function(i) {
                make_surrogate_annual_spline(block_yearday, 
                                             simplex_results[[i, "data"]]$abundance, 
                                             num_surr = 200)
            })
        saveRDS(simplex_results, simplex_results_file)
    }
    simplex_results <- readRDS(simplex_results_file)
    
    if (!file.exists(ccm_results_file))
    {
        max_ts_length <- max(map_int(simplex_results$data, NROW))
        lib_vec <- c(6, 12, 24, 40, 80, 140, 220, 320, max_ts_length)
        ccm_results <- compute_CCM(simplex_results, lib_sizes = lib_vec, 
                                   num_samples = 200, num_cores = num_cores)
        saveRDS(ccm_results, ccm_results_file)
    }
    ccm_results <- readRDS(ccm_results_file)
    
    if (!file.exists(smap_matrices_file))
    {
        ccm_links <- identify_ccm_links(ccm_results)
        smap_coeffs <- compute_smap_coeffs(block, ccm_links)
        smap_matrices <- compute_smap_matrices(smap_coeffs, ccm_links)
        saveRDS(smap_matrices, smap_matrices_file)
    }
    
    if (!file.exists(eigenvalues_file))
    {
        eigenvalues <- compute_eigenvalues(smap_matrices)
        saveRDS(eigenvalues, eigenvalues_file)
    }
}

#### function to produce eigenvalues plot ----

plot_eigenvalues <- function(proc_data_file = here::here("data", "portal_block.RDS"), 
                             eigenvalues_file = here::here("output", "portal_eigenvalues.RDS"), 
                             plot_file = NULL, width = 8, height = 4.5, highlight_shifts = FALSE)
{
    # read in data
    eigenvalues <- readRDS(eigenvalues_file)
    block <- readRDS(proc_data_file)
    portal_newmoons_table <- portalr::load_data()$newmoons_table %>%
        mutate(censusdate = as.Date(censusdate))
    
    # generate df for plotting
    stability_df <- data.frame(censusdate = block$censusdate, 
                               ds = map_dbl(eigenvalues, ~max(abs({.})))) %>%
        left_join(portal_newmoons_table, by = "censusdate")
    
    # construct plot
    portal_ds_plot <- stability_df %>%
        ggplot(aes(x = censusdate, y = ds)) + 
        geom_line() + 
        geom_hline(yintercept = 1.0, size = 1, linetype = 2) + 
        scale_x_date(breaks = seq(from = as.Date("1985-01-01"), 
                                  to = as.Date("2015-01-01"), 
                                  by = "5 years"), 
                     date_labels = "%Y", expand = c(0.01, 0)) + 
        labs(x = NULL, y = "dynamic stability \n(higher is more unstable)") +
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

#### function to produce smap coeffs plot ----

plot_smap_coeffs <- function(proc_data_file = here::here("data", "portal_block.RDS"), 
                             smap_matrices_file = here::here("output", "portal_smap_matrices.RDS"), 
                             plot_file = NULL, width = 6, height = NULL)
{
    # make data.frame of smap coefficients
    smap_matrices <- readRDS(smap_matrices_file)
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
    block <- readRDS(proc_data_file)
    stopifnot(max(smap_coeff_df$t) <= NROW(block))
    smap_coeff_df$censusdate <- as.Date(block$censusdate[smap_coeff_df$t])
    
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
    
    combined_plot <- cowplot::plot_grid(ts_plot, density_plot, nrow = 1, 
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

#### base run ----
#  using all species in the dataset and default options
proc_data_file <- here::here("data", "portal_block.RDS")
simplex_results_file <- here::here("output", "portal_simplex_results.RDS")
ccm_results_file <- here::here("output", "portal_ccm_results.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices.RDS")
smap_plot_file <- here::here("output", "portal_smap_values.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues.RDS")
eigenvalues_plot_file <- here::here("output", "portal_dynamic_stability.pdf")

compute_portal_dynamic_stability(proc_data_file, 
                                 simplex_results_file, 
                                 ccm_results_file, 
                                 smap_matrices_file, 
                                 eigenvalues_file)

plot_eigenvalues(proc_data_file, 
                 eigenvalues_file, 
                 eigenvalues_plot_file)

plot_smap_coeffs(proc_data_file, 
                 smap_matrices_file, 
                 smap_plot_file)

#### revised run (50%) ----
#  only species which appear at least 50% of the time)
#  DM, DO, DS, NA, OL, OT, PP

compute_portal_dynamic_stability(filter_q = 0.5, 
                                 num_cores = 1)

proc_data_file <- here::here("data", "portal_block_50.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices_50.RDS")
smap_plot_file <- here::here("output", "portal_smap_values_50.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues_50.RDS")
eigenvalues_plot_file <- here::here("output", "portal_dynamic_stability_50.pdf")

plot_eigenvalues(proc_data_file, 
                 eigenvalues_file, 
                 eigenvalues_plot_file)

plot_eigenvalues(proc_data_file, 
                 eigenvalues_file, 
                 plot_file = here::here("output", "portal_dynamic_stability_50_highlight.pdf"), 
                 highlight_shifts = TRUE)

plot_smap_coeffs(proc_data_file, 
                 smap_matrices_file, 
                 smap_plot_file)

#### revised run (33%) ----
#  only species which appear at least 33% of the time)
#  DM, DO, DS, NA, OL, OT, PB, PE, PF, PP, RM

compute_portal_dynamic_stability(filter_q = 1/3,
                                 num_cores = 1)

proc_data_file <- here::here("data", "portal_block_33.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices_33.RDS")
smap_plot_file <- here::here("output", "portal_smap_values_33.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues_33.RDS")
eigenvalues_plot_file <- here::here("output", "portal_dynamic_stability_33.pdf")

plot_eigenvalues(proc_data_file, 
                 eigenvalues_file, 
                 eigenvalues_plot_file)

plot_smap_coeffs(proc_data_file, 
                 smap_matrices_file, 
                 smap_plot_file)

#### run on LDA topics ----

proc_data_file <- here::here("data", "portal_block_lda.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices_lda.RDS")
smap_plot_file <- here::here("output", "portal_smap_values_lda.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues_lda.RDS")
plot_file <- here::here("output", "portal_dynamic_stability_lda.pdf")

# do our own block generation using LDA topics
if (!file.exists(proc_data_file))
{
    raw_rodent_data <- portalr::abundance(time = "all", 
                                          level = "plot", 
                                          effort = TRUE, 
                                          na_drop = TRUE)
    
    # summarize by each newmmonnumber, and for only the control plots we want
    block <- raw_rodent_data %>% 
        filter(plot %in% c(2, 4, 8, 11, 12, 14, 17, 22), 
               censusdate < "2015-04-18") %>% 
        select(-treatment, -plot, -period) %>%
        group_by(newmoonnumber, censusdate) %>%
        summarize_all(sum) %>%
        ungroup()
    
    # check that effort is equal across samples
    stopifnot(length(unique(block$ntraps)) == 1)
    
    # generate LDA
    LDA_models <- LDATS::parLDA(data = select(block, BA:SO), ntopics = 2:5, nseeds = 200, ncores = 2)
    best_LDA <- LDATS::LDA_select(LDA_models, correction = TRUE)
    topic_block <- data.frame(best_LDA@gamma)
    names(topic_block) <- paste0("Topic", seq(NCOL(topic_block)))
    topic_block <- cbind(select(block, -(BA:SO)), topic_block)
    
    # add in NAs for unsampled newmoonnumbers and interpolate
    block <- topic_block %>%
        complete(newmoonnumber = full_seq(newmoonnumber, 1), fill = list(NA)) %>%
        mutate_at(vars(-newmoonnumber, -ntraps), forecast::na.interp) %>%
        mutate(censusdate = as.Date(as.numeric(censusdate), origin = "1970-01-01")) %>% 
        select(-newmoonnumber, -ntraps)
    
    saveRDS(block, proc_data_file)
}

compute_portal_dynamic_stability(file_suffix = "_lda")

plot_eigenvalues(proc_data_file, 
                 eigenvalues_file, 
                 eigenvalues_plot_file)

plot_smap_coeffs(proc_data_file, 
                 smap_matrices_file, 
                 smap_plot_file)

