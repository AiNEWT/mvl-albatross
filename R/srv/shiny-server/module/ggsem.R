library(tidyverse)
library(tidygraph)
library(ggraph)
library(lavaan)

ggsem <- function(fit, layout = "kk") {

  # Extract standardized parameters
  params <- lavaan::standardizedSolution(fit)

  # Edge properties
  param_edges <- params %>% 
    filter(op %in% c("=~", "~", "~~"), lhs != rhs, pvalue < .10) %>%
    transmute(to = lhs,
              from = rhs,
              val = est.std,
              type = dplyr::case_when(
                op == "=~" ~ "loading",
                op == "~"  ~ "regression",
                op == "~~" ~ "correlation",
                TRUE ~ NA_character_))

  # Identify latent variables for nodes
  latent_nodes <- param_edges %>% 
    filter(type == "loading") %>% 
    distinct(to) %>% 
    transmute(metric = to, latent = TRUE)

  # Node properties
  param_nodes <- params %>% 
    filter(lhs == rhs) %>% 
    transmute(metric = lhs, e = est.std) %>% 
    left_join(latent_nodes) %>% 
    mutate(latent = if_else(is.na(latent), FALSE, latent))

  # Complete Graph Object
  param_graph <- tidygraph::tbl_graph(param_nodes, param_edges)

  # Plot
  ggraph(param_graph, layout = layout) +
    # Latent factor Nodes
    geom_node_point(aes(alpha = as.numeric(latent)),
                                     shape = 16, size = 5) +
    geom_node_point(aes(alpha = as.numeric(latent)),
                                     shape = 16, size = 4, color = "white") +
    # Observed Nodes
    geom_node_point(aes(alpha = as.numeric(!latent)),
                                     shape = 15, size = 5) +
    geom_node_point(aes(alpha = as.numeric(!latent)),
                                     shape = 15, size = 4, color = "white") +
    # Regression Paths (and text)
    geom_edge_link(aes(color = val, label = round(val, 2),
                       alpha = as.numeric(type == "regression")),
                   linetype = 1, angle_calc = "along", vjust = -.5,
                   arrow = arrow(20, unit(.3, "cm"), type = "closed")) +
    # Factor Loadings (no text)
    geom_edge_link(aes(color = val, alpha = as.numeric(type == "loading")),
                   linetype = 3, angle_calc = "along",
                   arrow = arrow(20, unit(.3, "cm"), ends = "first", type = "closed")) +
    # Correlation Paths (no text)
    geom_edge_link(aes(color = val, alpha = as.numeric(type == "correlation")),
                   linetype = 2, angle_calc = "along",
                   arrow = arrow(20, unit(.3, "cm"), type = "closed", ends = "both")) +
    # Node names
    geom_node_text(aes(label = metric),
                   nudge_y = .25, hjust = "inward") +
    # Node residual error
    geom_node_text(aes(label = sprintf("%.2f", e)),
                   nudge_y = -.1, size = 3) +
    # Scales and themes
    scale_alpha(guide = FALSE, range = c(0, 1)) +
    scale_edge_alpha(guide = FALSE, range = c(0, 1)) +
    scale_edge_colour_gradient2(guide = FALSE, low = "red", mid = "darkgray", high = "green") +
    scale_edge_linetype(guide = FALSE) +
    scale_size(guide = FALSE) +
    theme_graph()
}
