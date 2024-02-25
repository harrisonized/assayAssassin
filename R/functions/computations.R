import::here(magrittr, '%>%')
import::here(rlang, 'syms')
import::here(dplyr, 'group_by', 'summarize')
import::here(tidyr, 'pivot_wider')

## Functions
## ct_thresholds_from_results
## compute_dct_table


#' Compute dCT Table
#' 
ct_thresholds_from_results <- function(results) {
    ct_thresholds <- results[, c('plate_id', 'ct_threshold')] %>%
        group_by(plate_id) %>%
        summarize(ct_threshold=first(ct_threshold))
    ct_thresholds <- list2env(setNames(
        nm=ct_thresholds[['plate_id']],  # keys
        object=as.list(ct_thresholds[['ct_threshold']])  # values
    ))
return(ct_thresholds)
}


#' Compute dCT Table
#' 
dct_table_from_results <- function(
    results,
    index_cols=c('plate_id', 'sample_id', 'tissue', 'gene'),
    sample_genes=c('Dnase1l1'),
    control_genes=c('Actin', 'Hprt')
) {

    ct_long <- results %>%
        group_by(!!!syms(index_cols)) %>%
        summarize(mean_ct=mean(ct, na.rm=TRUE),
                  stdev_ct=sd(ct, na.rm=TRUE),
                  .groups = 'drop')
    ct_wide <- pivot_wider(
        ct_long,
        names_from=c(gene),
        values_from=c(mean_ct, stdev_ct)
    )

    for (sample_gene in sample_genes) {
        for (control_gene in control_genes) {
            sample_gene <- tolower(sample_gene)
            control_gene <- tolower(control_gene)
            dct_col <- paste('dct', sample_gene, control_gene, sep='_')
            fold_change_col <- paste('fold_change', sample_gene, control_gene, sep='_')
            sample_gene_col <- paste('mean_ct', sample_gene, sep='_')
            control_gene_col <- paste('mean_ct', control_gene, sep='_')

            ct_wide[[dct_col]] = ct_wide[[sample_gene_col]] - ct_wide[[control_gene_col]]
            ct_wide[[fold_change_col]] = 2^(-ct_wide[[dct_col]])
        }
    }

    return(ct_wide)  # aka dct_table
}
