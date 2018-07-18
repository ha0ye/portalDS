library(ggraph)
library(igraph)

source(here::here("R", "dynamic_stability_functions.R"))

make_network_from_ccm <- function(ccm_results, layout = "circle")
{
    ccm_results %>%
        identify_ccm_links() %>%
        arrange(target_column) %>%
        select(target_column, lib_column) %>%
        graph_from_data_frame() %>%
        ggraph(layout = layout) + 
        geom_edge_fan(edge_width = 0.5) + 
        geom_node_circle(aes(r = 0.05, fill = name)) + 
        theme_graph(foreground = "black", fg_text_colour = "white", 
                    background = "transparent") + 
        coord_fixed() + 
        theme(text = element_text(family = "Helvetica"),
              panel.border = element_rect(color = NA, fill = NA))
}