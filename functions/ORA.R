#' @title ORA_fisher
#' 
#' @description Runs overrepresentation analysis on two gene lists
#' 
#' @param gene_list1 first list of genes to test
#' @param gene_list2 second list of genes to test
#' @param background_list1 background for genes from list1
#' @param background_list2 background for genes from list2
#' @param alternative "greater" for over-representation; "less" for under-representation, "two-sided" for both
#' 
#' @export
ORA_fisher <- function(
    gene_list1,
    gene_list2,
    background_list1,
    background_list2,
    alternative = "greater" # default to one-tailed over-representation
) {
    
    # Ensure each list is unique
    gene_list1 <- unique(gene_list1)
    gene_list2 <- unique(gene_list2)
    background_list1 <- unique(background_list1)
    background_list2 <- unique(background_list2)
    
    # Overlap counts
    q <- length(intersect(gene_list1, gene_list2)) # overlap
    k <- length(intersect(gene_list2, background_list1)) # list2 in list1 background
    m <- length(intersect(gene_list1, background_list2)) # list1 in list2 background
    t <- length(intersect(background_list1, background_list2)) # universe
    
    ## Build 2x2 contingency table
    ##           inList1     !List1
    ## inList2      q         k - q
    ## !inList2   m - q    t - m - k + q
    
    mat <- matrix(
        c(q,
          k - q,
          m - q,
          t - m - k + q),
        nrow = 2,
        byrow = TRUE
    )
    
    # Run Fisher exact test
    fisher_out <- fisher.test(mat, conf.int = TRUE, alternative = alternative)
    
    # Return output as tibble
    tibble(
        odds_ratio = unname(fisher_out$estimate),
        ci_lower = fisher_out$conf.int[1],
        ci_upper = fisher_out$conf.int[2],
        fold_enrichment = (q / m) / (k / t), # observed fraction over expected fraction
        p_value = fisher_out$p.value,
        n_overlap = q,
        n_list1 = m,
        n_list2 = k,
        n_background = t
    )
    
    
}
