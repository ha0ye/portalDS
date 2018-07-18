smap_matrices <- readRDS("output/portal_smap_matrices_50.RDS")

eigenvalue_dist <- map_df(seq(smap_matrices), function(i) {
    m <- smap_matrices[[i]]
    if(any(is.na(m)))
        return(data.frame())
    eigenvalues <- sort(abs(eigen(m)$values), decreasing = TRUE)
    data.frame(eigenvalues = eigenvalues, t = i, rank = seq(eigenvalues))
})

eigenvalue_dist %>%
    filter(rank < 10) %>%
    ggplot(aes(x = t, y = eigenvalues, color = as.factor(rank))) + 
    geom_line() + 
    scale_color_viridis_d() + 
    theme_bw()
