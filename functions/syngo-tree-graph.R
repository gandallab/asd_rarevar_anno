#===============================================================================
# Build a pruned ggraph tree from a df_cmu (SynGO CMU) domain subset
# Relocated from 05_00_SynGO-ms-fig.qmd; used by 05_SynGO analysis and Fig3.
#===============================================================================

# Build a directed tbl_graph (parent → child) from a df_cmu domain subset.
# Pruning follows Bernie's criterion: show a node if Z_all > 2.8 OR it is a
# direct child of the root (n_gen <= 2). HotNet neighborhood nodes are always
# included regardless of Z_all.
make_syngo_tbl_graph <- function(df_ann, z_thresh = 2.8) {
    hotnet_ids <- df_ann$id[!is.na(df_ann$hotnet_group) & df_ann$hotnet_group != "None"]
    nodes <- df_ann %>%
        mutate(z_all_num = as.numeric(z_all)) %>%
        filter(n_gen <= 2 | (if (!is.null(z_thresh)) abs(z_all_num) > z_thresh else TRUE) | id %in% hotnet_ids) %>%
        mutate(
            adj_z      = as.numeric(adj_z_all),
            is_hotnet  = !is.na(hotnet_group) & hotnet_group != "None",
            # Full GO term name, stripping the "(GO:XXXXXXX)" or "(SYNGO:xxx)" suffix
            name_clean = str_remove(name, "\\s*\\([A-Z]+:[^)]+\\)$"),
            node_idx   = row_number()
        )

    edges <- nodes %>%
        filter(!is.na(parent_id), parent_id %in% nodes$id) %>%
        transmute(
            from = match(parent_id, nodes$id),
            to   = node_idx
        )

    # Use tidygraph's .N() accessor to attach node attributes to edges after
    # graph creation — pre-computing in transmute() doesn't expose them to ggraph
    tbl_graph(nodes = nodes, edges = edges, directed = TRUE) %>%
        activate(edges) %>%
        mutate(to_is_hotnet = .N()$is_hotnet[to]) %>%
        activate(nodes)
}
