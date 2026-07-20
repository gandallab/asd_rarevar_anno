#===============================================================================
# GO cosine-similarity network plotting ("ball-and-stick" graphs)
#
# Shared code for the GR / MORPH / MEF2C GO-similarity graphs in the paper.
# Each call to plot_go_graph() plots one network (one subcluster or one regulon);
# the calling notebook builds the per-cluster graphs (go_knn_edges) and combines
# panels (patchwork), adds titles / layout, etc.
#===============================================================================

# Within-set kNN edges from a cosine-similarity matrix (deduped, undirected).
# Each gene links to its k_nn most-similar partners within the matrix.
go_knn_edges <- function(sim_mat, k_nn = 3) {
    genes <- rownames(sim_mat)
    map_dfr(genes, function(g) {
        others <- setdiff(genes, g)
        sims   <- sort(sim_mat[g, others], decreasing = TRUE)
        if (length(sims) == 0)
            return(tibble(from = character(), to = character(), similarity = numeric()))
        tibble(from       = g,
               to         = names(sims)[seq_len(min(k_nn, length(sims)))],
               similarity = sims[seq_len(min(k_nn, length(sims)))])
    }) %>%
        mutate(key = map2_chr(from, to, ~ paste(sort(c(.x, .y)), collapse = "_"))) %>%
        distinct(key, .keep_all = TRUE) %>%
        dplyr::select(-key)
}

# Plot one GO-similarity network (GR subcluster / MORPH subcluster / MEF2C regulon).
# A graphopt force-directed layout; nodes filled by GO topic, sized by log10(BF); edge
# width scales with GO cosine similarity; hub/TF genes get a bold outline + bold label.
#   net          : igraph object (built from go_knn_edges())
#   node_stats   : data frame with gene_name, GO_topic, and log_bf (= log10(BF)) columns
#   bold_genes   : character vector of genes to draw with a bold black outline + bold label
#   title        : panel title
#   size_range   : node size aesthetic range (e.g. c(1, 4) for GR/MORPH, c(2, 6) for MEF2C)
#   size_limits  : shared log10(BF) limits so node sizes are comparable across panels
#   edge_limits  : shared similarity limits so edge widths are comparable across panels
#   topic_legend : show the GO-topic fill legend (TRUE for the multi-topic MEF2C plot;
#                  FALSE for single-topic GR/MORPH, where it would be a redundant one entry)
#   seed         : layout seed (graphopt is stochastic)
#   xmult/ymult  : per-panel padding inside the panel (horizontal / vertical whitespace)
#   fill_values  : named vector of GO-topic -> colour (defaults to the project topic_colors)
plot_go_graph <- function(net, node_stats, bold_genes = character(0), title = NULL,
                          size_range = c(1, 4), size_limits = NULL, 
                          edge_range = c(0.025, 0.40), edge_limits = NULL, edge_color = "gray50",
                          topic_legend = FALSE, seed = 1214, xmult = 0.15, ymult = 0.15,
                          fill_values = topic_colors) {
    set.seed(seed)
    lay <- create_layout(net, layout = "graphopt", charge = 0.04)
    lay$x <- lay$x * 2
    lay$y <- lay$y * 2

    lay %>%
        left_join(node_stats, by = join_by(name == gene_name)) %>%
        mutate(hub_tf = name %in% bold_genes) %>%
        
        ggraph() +
        geom_edge_link0(aes(linewidth = similarity), color = edge_color) +
        geom_node_point(aes(size = log_bf, fill = GO_topic, color = hub_tf, stroke = hub_tf),
                        shape = 21) +
        geom_node_text(aes(label = name, fontface = if_else(hub_tf, "bold", "plain")),
                       vjust = -0.8, size = 1.6) +

        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray90"),
                           guide = "none") +   # hubs shown by black border + bold label
        scale_discrete_manual(aesthetics = "stroke",
                              values = c("TRUE" = 0.6, "FALSE" = 0.2),
                              guide = "none") +
        scale_fill_manual(values = fill_values, na.value = "gray70", name = "GO topic",
                          guide = if (topic_legend) "legend" else "none") +
        scale_size_continuous(range = size_range, limits = size_limits, name = "log10(BF)") +
        scale_edge_width_continuous(range = edge_range, limits = edge_limits) +
        # per-panel whitespace: bigger mult = more padding, network occupies less of the panel
        scale_x_continuous(expand = expansion(mult = xmult)) +
        scale_y_continuous(expand = expansion(mult = ymult)) +
        theme_void() +
        coord_cartesian(clip = "off") +
        labs(title = title) +
        theme(
            plot.title = element_text(size = 7, hjust = 0.5, face = "bold", margin = margin(b = 5)),
            legend.key.size = unit(0.25, "cm"),
            legend.title = element_text(size = 5.5),
            legend.text = element_text(size = 5)
        )
}
