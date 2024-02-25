import::here(magrittr, '%>%')
import::here(rlang, 'syms')
import::here(dplyr, 'group_by', 'summarize')
import::here(tidyr, 'pivot_wider')


#' Compute dCT Table
#' 
compute_dct_table <- function(
    results,
    index_cols=c('sample_id', 'tissue', 'gene')
) {

    ct_long <- results %>%
        group_by(!!!syms(index_cols)) %>%
        summarize(mean_ct=mean(ct, na.rm=TRUE),
                  stdev_ct=sd(ct, na.rm=TRUE),
                  .groups = 'drop')
    ct_long <- ct_long[(ct_long[['tissue']]!=''),]

    ct_wide <- pivot_wider(
        ct_long,
        names_from=c(gene),
        values_from=c(mean_ct, stdev_ct)
    )

    # TODO: derive this from the colnames
    sample_genes <- c('dnase1l1')
    control_genes <- c('actin', 'hprt')

    for (sample_gene in sample_genes) {
        for (control_gene in control_genes) {
            dct_col <- paste('dct', sample_gene, control_gene, sep='_')
            fold_change_col <- paste('fold_change', sample_gene, control_gene, sep='_')
            sample_gene_col <- paste('mean_ct', sample_gene, sep='_')
            control_gene_col <- paste('mean_ct', control_gene, sep='_')

            ct_wide[[dct_col]] = ct_wide[[sample_gene_col]] - ct_wide[[control_gene_col]]
            ct_wide[[fold_change_col]] = 2^(-ct_wide[[dct_col]])
        }
    }

    return(ct_wide)
}
