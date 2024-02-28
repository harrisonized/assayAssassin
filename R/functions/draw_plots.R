import::here(rlang, 'sym')
import::here(dplyr, 'group_by', 'summarize')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'df_to_plate', 'reset_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'savefig', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_heatmap', 'plot_scatter', 'plot_lines', 'plot_dots_and_bars',
    .character_only=TRUE)

## Functions
## draw_cq_conf_scatter
## draw_cq_conf_dots
## draw_heatmaps
## draw_amp_curves
## draw_fold_changes


#' Draw CQ Conf Scatter
#' 
draw_cq_conf_scatter <- function(
    results,
    dirpath,
    value='ct',
    troubleshooting=FALSE,
    showfig=FALSE
) {
    plate_ids <- sort(unique(results[['plate_id']]))
    for (plate_id in plate_ids) {

        fig <- plot_scatter(
            results[(results['plate_id']==plate_id), ],
            x=value,
            y='cq_conf',
            color='tissue',
            xlabel=value,
            ylabel='cq confidence',
            title=paste('Quality Control for', value)
        )
        if (showfig) { print(fig) }
        savefig(file.path(dirpath, 'qc', plate_id, paste0('scatter-cq_conf-', value, '-', plate_id, '.png')),
                dpi=400,
                troubleshooting=troubleshooting)
    }
}


#' Draw Cq Conf Dots
#' 
draw_cq_conf_dots <- function(
    results,
    genes=c('Dnase1l1', 'Actin', 'Hprt'),
    dirpath,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    for (gene in genes) {

        tmp <- results[(results['gene']==tolower(gene)), ]

        # compute barsize
        barsize <- tmp %>%
            group_by(tissue) %>%
            summarize(min_bar=min(!!sym('cq_conf'), na.rm=TRUE),
                      mean_bar=mean(!!sym('cq_conf'), na.rm=TRUE),
                      max_bar=max(!!sym('cq_conf'), na.rm=TRUE))
        tmp <- merge(
            tmp, barsize,
            by='tissue', all.x=TRUE, all.y=FALSE 
        )

        fig <- plot_dots_and_bars(
            tmp[order(tmp[['mean_bar']], decreasing = TRUE),],  # sort
            x='tissue',
            y='cq_conf',
            color='sample_id',
            ymin='min_bar',
            ymax='max_bar',
            log_y=FALSE
        )
        if (showfig) { print(fig) }
        savefig(file.path(wd, opt[['figures-dir']], 'qc', paste0('cq_conf-', gene, '.png')),
                width=1600, dpi=400,
                troubleshooting=troubleshooting)
    }
}


#' Draw Heatmaps
#' 
draw_heatmaps <- function(
    results,
    dirpath,
    value='ct',
    min_cq_conf=0.75,
    troubleshooting=FALSE,
    showfig=FALSE
) {
    plate_ids <- sort(unique(results[['plate_id']]))
    for (plate_id in plate_ids) {
        tmp <- results
        tmp[(tmp[['cq_conf']] <= min_cq_conf), 'ct'] <- NA
        plate <- df_to_plate(
            tmp[(tmp['plate_id']==plate_id), ],
            value=value,
            num_wells=96
        )

        fig <- plot_heatmap(plate,
            show_xlabel=FALSE,
            show_ylabel=FALSE,
            title=paste(value, 'for plate', plate_id),
            annotations=TRUE, digits=2)
        if (showfig) { print(fig) }
        savefig(file.path(dirpath, 'qc', plate_id, paste0('heatmap-', value, '-', plate_id, '.png')),
                dpi=400,
                troubleshooting=troubleshooting)
    }
}


#' Draw Amp Curves
#' 
draw_amp_curves <- function(
    amp_data,
    dirpath,
    ct_thresholds=NULL,
    metadata_cols=c("tissue", "gene"),
    troubleshooting=FALSE,
    showfig=FALSE
) {
    plate_ids <- sort(unique(amp_data[['plate_id']]))
    combinations <- expand.grid(
        plate_id=plate_ids,
        colname=metadata_cols,
        stringsAsFactors=FALSE
    )

    for (row in rownames(combinations)) {
        plate_id <- combinations[row, c('plate_id')]
        group <- combinations[row, c('colname')]
        ct_threshold <- ct_thresholds[[plate_id]]

        subset <- amp_data[(amp_data['plate_id']==plate_id), ]

        fig <- plot_lines(
            subset,
            x='cycle',
            y='delta_rn',
            group='row_id',
            color=group,
            thresholds=ct_threshold,
            xlabel='Cycle Number',
            ylabel='Delta Rn',
            title=paste('Amplification Data for plate', plate_id)
        )
        if (showfig) { print(fig) }
        savefig(file.path(dirpath, 'qc', plate_id, 'amp_data', paste0('delta_rn-', group, '-', plate_id, '.png')),
                height=1000, width=1600, dpi=400,
                troubleshooting=troubleshooting)
    }
}


#' Draw Fold Changes
#' 
draw_fold_changes <- function(
    dct_table,
    sample_genes=c('Dnase1l1'),
    control_genes=c('Actin', 'Hprt'),
    dirpath,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    combinations <- expand.grid(
        sample_gene=sample_genes,
        control_gene=control_genes,
        stringsAsFactors=FALSE
    )

    for (row in rownames(combinations)) {
        sample_gene_name <- combinations[row, c('sample_gene')]
        control_gene_name <- combinations[row, c('control_gene')]
        fold_change_col <- paste(
            'fold_change', tolower(sample_gene_name), tolower(control_gene_name), sep='_')

        # compute barsize
        barsize <- dct_table %>%
            group_by(tissue) %>%
            summarize(min_fold_change=min(!!sym(fold_change_col), na.rm=TRUE),
                      mean_fold_change=mean(!!sym(fold_change_col), na.rm=TRUE),
                      max_fold_change=max(!!sym(fold_change_col), na.rm=TRUE))
        tmp <- merge(
            dct_table, barsize,
            by='tissue', all.x=TRUE, all.y=FALSE 
        )

        fig <- plot_dots_and_bars(
            tmp[order(tmp[['mean_fold_change']], decreasing = TRUE),],  # sort
            x='tissue',
            y=fold_change_col,
            color='sample_id',
            title=paste("Relative Expression of", sample_gene_name, "vs.", control_gene_name)
        )
        if (showfig) { print(fig) }
        savefig(file.path(wd, opt[['figures-dir']], paste0(fold_change_col, '.png')),
                width=1600, dpi=400,
                troubleshooting=troubleshooting)
    }
    
}
