import::here(rlang, 'sym')
import::here(dplyr, 'group_by', 'summarize')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'df_to_plate', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'savefig', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_heatmap', 'plot_amp_curves', 'plot_fold_change', .character_only=TRUE)

## Functions
## draw_ct_heatmaps
## draw_amp_curves
## draw_fold_changes


#' Draw CT Heatmaps
#' 
draw_ct_heatmaps <- function(
    results,
    plate_ids,
    dirpath,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    for (plate_id in plate_ids) {

        plate <- df_to_plate(
            results[(results['plate_id']==plate_id), ],
            value='ct'
        )

        fig <- plot_heatmap(plate,
            show_xlabel=FALSE,
            show_ylabel=FALSE,
            title=paste('CT for Plate ', plate_id), annotations=TRUE, digits=2)
        if (showfig) { print(fig) }
        savefig(file.path(dirpath, plate_id, paste0('heatmap-ct-', plate_id, '.png')),
                dpi=400,
                troubleshooting=troubleshooting)
    }
}


#' Draw Amp Curves
#' 
draw_amp_curves <- function(
    amp_data,
    plate_ids,
    dirpath,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    amp_data <- amp_data[(amp_data[['gene']]!='') & (amp_data[['cq_conf']] > 0.5), ]

    for (plate_id in plate_ids) {
        for (group in c('gene', 'tissue')) {

            ct_threshold <- results[(results['plate_id']==plate_id), 'ct_threshold'][[1]]
            subset <- amp_data[(amp_data['plate_id']==plate_id), ]
            
            fig <- plot_amp_curves(subset, ct_threshold, color=group, plate_id= plate_id)
            if (showfig) { print(fig) }
            savefig(file.path(dirpath, plate_id, paste0('delta_rn-', group, '-', plate_id, '.png')),
                    height=1000, width=1600, dpi=400,
                    troubleshooting=troubleshooting)
        }
    }
}


draw_fold_changes <- function(
    ct_wide,
    control_genes=c('Actin', 'Hprt'),
    dirpath,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    for (control_gene_name in control_genes) {
        fold_change_col <- paste0('fold_change_dnase1l1_', tolower(control_gene_name))

        ct_wide <- ct_wide[order(ct_wide[[fold_change_col]], decreasing = TRUE),]
        ct_wide <- reset_index(ct_wide, drop=TRUE)

        # manual filter
        ct_wide <- ct_wide[(ct_wide[fold_change_col] <= 1), ]

        # barsize
        barsize <- ct_wide %>%
            group_by(tissue) %>%
            summarize(min_fold_change=min(!!sym(fold_change_col), na.rm=TRUE),
                      max_fold_change=max(!!sym(fold_change_col), na.rm=TRUE))
        tmp <- merge(ct_wide, barsize,
            by='tissue', all.x=TRUE, all.y=FALSE 
        )

        tmp <- tmp[!is.na(tmp['sample_id']), ]

        fig <- plot_fold_change(
            tmp,
            x='tissue',
            y=fold_change_col,
            color='sample_id',
            title=paste("Relative Expression of Dnase1l1 vs.", control_gene_name)
        )
        if (showfig) { print(fig) }
        savefig(file.path(wd, opt[['figures-dir']], paste0(fold_change_col, '.png')),
                width=1600, dpi=400,
                troubleshooting=troubleshooting)
    }
}
