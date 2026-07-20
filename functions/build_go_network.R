#===============================================================================
# build_go_network(): gene–gene GO cosine-similarity network construction
#
# Shared general network-construction helper used for both the module / regulon
# "hairball" plots (module-plots.qmd) and the GO-topic subcluster analyses
# (02_02). Construction ONLY, layout and plotting live in the calling scripts
#===============================================================================

# Cosine similarity between genes (rows of a gene × GO binary matrix).
# Genes with no GO annotation get an all-zero row -> similarity 0 to everything
cos_sim_genes <- function(mat) {
    norms <- sqrt(rowSums(mat^2))
    norms[norms == 0] <- NA
    mat_norm <- mat / norms
    sim <- mat_norm %*% t(mat_norm)
    sim[is.na(sim)] <- 0
    sim
}

# Build a gene–gene GO-similarity network from a character vector of genes.
#   genes      : character vector of gene symbols
#   go_mapping : tibble of GO annotation with columns `ensembl_gene_id`, `go_id`
#   gene_stats : tibble bridging gene symbols to Ensembl IDs and TADA FDR — must have
#                columns `gene_name`, `ensembl_gene_id`, `FDR`
#   keep_genes : genes to keep even if they share no GO terms (e.g. hubs / TFs)
#   sim_thresh : minimum cosine similarity to draw an edge (edges use similarity > sim_thresh)
#   k_nn       : number of nearest neighbors each gene connects to
#   risk_fdr   : FDR threshold below which a gene is kept even when isolated
# Returns a list:
#   sim_mat : full gene × gene cosine-similarity matrix, incl. no-GO genes (for clustering)
#   edges   : kNN edges (from, to, similarity) among genes sharing >= 1 GO term
#   nodes   : tibble of retained genes (connected genes + risk genes + keep_genes)
#   no_go   : genes with no GO annotation in `go_mapping`
build_go_network <- function(genes, go_mapping, gene_stats,
                             keep_genes = NULL, sim_thresh = 0, k_nn = 3, risk_fdr = 0.001) {
    genes <- unique(genes)

    # gene × GO binary matrix
    gene_go <- go_mapping %>%
        filter(ensembl_gene_id %in% {
            gene_stats %>% filter(gene_name %in% genes) %>% pull(ensembl_gene_id)
        }) %>%
        left_join(gene_stats %>% dplyr::select(ensembl_gene_id, gene_name),
                  by = "ensembl_gene_id") %>%
        distinct(gene_name, go_id)

    go_mat <- gene_go %>%
        mutate(present = 1L) %>%
        pivot_wider(id_cols = gene_name, names_from = go_id,
                    values_from = present, values_fill = 0L) %>%
        column_to_rownames("gene_name") %>%
        as.matrix()

    sim <- cos_sim_genes(go_mat)

    # genes with no GO annotation: add as zero-similarity nodes (self-similarity 1)
    no_go <- genes[!genes %in% rownames(sim)]
    if (length(no_go) > 0) {
        zero_rows <- matrix(0, length(no_go), nrow(sim), dimnames = list(no_go, rownames(sim)))
        zero_cols <- matrix(0, nrow(sim) + length(no_go), length(no_go),
                            dimnames = list(c(rownames(sim), no_go), no_go))
        sim <- cbind(rbind(sim, zero_rows), zero_cols)
        diag(sim) <- 1
    }

    # kNN edges: each gene links to its k_nn most-similar neighbors that
    # actually share >= 1 GO term (similarity strictly > sim_thresh)
    edges <- map_dfr(rownames(sim), function(g) {
        sims <- sort(sim[g, rownames(sim) != g], decreasing = TRUE)
        sims <- sims[sims > sim_thresh]
        if (length(sims) == 0) return(NULL)
        tibble(from = g, to = names(sims)[seq_len(min(k_nn, length(sims)))],
               similarity = sims[seq_len(min(k_nn, length(sims)))])
    }) %>%
        mutate(key = map2_chr(from, to, ~ paste(sort(c(.x, .y)), collapse = "_"))) %>%
        distinct(key, .keep_all = TRUE) %>%
        dplyr::select(-key)

    # keep connected genes + ASD risk genes + requested keep_genes (hubs / TFs)
    risk_genes <- gene_stats %>%
        filter(gene_name %in% genes, !is.na(FDR), FDR < risk_fdr) %>%
        pull(gene_name)
    keep <- union(unique(c(edges$from, edges$to)),
                  intersect(rownames(sim), c(risk_genes, keep_genes)))

    list(sim_mat = sim, edges = edges, nodes = tibble(name = keep), no_go = no_go)
}
