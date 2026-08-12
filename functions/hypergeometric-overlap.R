#' calculate_odds_ratio
#' 
#' Calculates the odds ratio and its 95% confidence interval (CI) for a 2x2 table.
#' Uses Woolf's method for SE calculation.
#'  
#' The odds ratio is the ratio of the odds of an event in one group to the odds in another
#' 
#' @param q Numeric. Overlap size - 1
#' @param m Numeric. Size of set 1
#' @param n Numeric. Complement of set 1 (N - m, where N is total population size)
#' @param k Numeric. Size of set 2
#' 
#' @return A list containing:
#'         \itemize{
#'              \item \strong{odds_ratio}: Numeric. The hypergeometric odds ratio.
#'              \item \strong{ci_lower}: Numeric. The lower bound of the 95% CI.
#'              \item \strong{ci_upper}: Numeric. The upper bound of the 95% CI.
#'         }
#' @export

calculate_odds_ratio <- function(q, m, n, k) {
    
    # Total population size
    N <- m + n
    
    # Set 2x2 contingency table components
    a <- q
    b <- k - q
    c <- m - q
    d <- N - (a + b + c)
    
    # Add pseudocount to avoid zero division
    a_adj <- a + 0.5
    b_adj <- b + 0.5
    c_adj <- c + 0.5
    d_adj <- d + 0.5
    
    # Calculate the odds_ratio & log transform
    odds_ratio <- (a_adj * d_adj) / (b_adj * c_adj)
    log_odds_ratio <- log(odds_ratio)
    
    # Calculate standard error of the log odds ratio (LOR) (Woolf's formula)
    se_log_odds_ratio <- sqrt(1/a_adj + 1/b_adj + 1/c_adj + 1/d_adj)
    
    # Calculate 95% CI for LOR (using Z = 1.96)
    z_score <- qnorm(0.975) 
    
    log_or_lower <- log_odds_ratio - z_score * se_log_odds_ratio
    log_or_upper <- log_odds_ratio + z_score * se_log_odds_ratio
    
    # Exponentiate
    ci_lower <- exp(log_or_lower)
    ci_upper <- exp(log_or_upper)
    
    # Return results
    return(list(
        odds_ratio = odds_ratio,
        ci_lower = ci_lower,
        ci_upper = ci_upper
    ))
    
}


#' calculate_hypergeometric_overlap
#' 
#' Calculates the hypergeometric overlap between two gene lists given a total gene pool ("universe")
#' 
#' @param gene_list1 A character vector of gene names (symbols) or IDs
#' @param gene_list2 A (different) character vector of gene names (symbols) or IDs. 
#'                   Gene identifier must match between gene_list1 and gene_list2 (i.e., both symbols or IDs)
#' @param gene_universe The total intersecting gene pool that the gene lists were derived from
#' 
#' @return A list containing the following:
#'         \itemize{
#'              \item \strong{intersect_genes}: Character. A vector representing the genes shared between the two lists
#'              \item \strong{p_value}: Numeric. The hypergeometric p-value of the overlap
#'              \item \strong{odds_ratio}: Numeric. The hypergeometric odds ratio (derived from calculate_odds_ratio)
#'         }
#' @export

calculate_hypergeometric_overlap <- function(gene_list1, gene_list2, gene_universe) {
    
    # Ensure that both gene lists are contained in the universe
    gene_list1 <- intersect(gene_list1, gene_universe) %>% unique
    gene_list2 <- intersect(gene_list2, gene_universe) %>% unique
    
    # Identify the overlap between the provided gene lists
    intersect_genes <- intersect(gene_list1, gene_list2)
    
    # Set variables for the hypergeometric test
    q <- length(intersect_genes)
    m <- length(gene_list2)
    n <- length(gene_universe) - m
    k <- length(gene_list1)
    
    # Run hypergeometric test
    log_pval <- phyper(q = q - 1, 
                       m = m, 
                       n = n, 
                       k = k,
                       lower.tail = FALSE,
                       log.p = TRUE
    )
    p_value <- exp(log_pval)
    
    # Calculate hypergeometric odds ratio
    odds_ratio <- calculate_odds_ratio(q, m, n, k)
    
    # Return a list of common genes and the hypergeometric p-value
    return(
        list("intersect_genes" = intersect_genes, 
             "p_value" = p_value, 
             "odds_ratio" = odds_ratio)
    )
    
}
