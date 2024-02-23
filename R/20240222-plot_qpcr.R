## Basic qPCR analysis

wd = dirname(this.path::here())  # wd = '~/github/R/leige-waffle'
suppressPackageStartupMessages(library('dplyr'))
# suppressPackageStartupMessages(library('ggplot2'))
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(readxl, 'read_excel')
import::from(tidyr, 'pivot_wider')
import::from(ggplot2,
    'ggplot', 'aes', 'theme', 'labs', 'geom_rect', 'geom_jitter', 'geom_line', 'geom_hline',
    'guide_axis', 'scale_x_discrete', 'scale_y_continuous', 'scale_fill_manual', 'element_text',
    'stage', 'ylim')

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'df_to_plate', 'set_index', 'reset_index', 'rev_df', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', 'read_excel_or_csv', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'move_list_item_to_start', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_bar', 'plot_heatmap', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input"), default='data/input',
                metavar='data/input',
                type="character",help="path/to/input/dir"),
   
    make_option(c("-o", "--output"), default="figures/20240222-third-pass",
                metavar="figures", type="character",
                help="set the output directory for the data"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("20240222-plot_qpcr-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read data

log_print(paste(Sys.time(), 'Reading data...'))

input_dirs <- items_in_a_not_b(
    list.dirs(file.path(wd, opt[['input']]), full.names=FALSE), ''
)  # folders are samples under this scheme

result_dfs <- new.env()
amp_data_dfs <- new.env()
for (input_dir in input_dirs) {

    sample_name <- input_dir  # may not be true in future experiments
    qpcr_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='xls')[[1]]

    # plate setup file
    plate_setup_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='csv')
    plate_setup <- read_excel_or_csv(plate_setup_file)
    plate_setup['sample_name'] <- sample_name
    plate_setup <- plate_setup[, 
        move_list_item_to_start(colnames(plate_setup), 'sample_name')
    ]

    # result
    result <- read_excel(qpcr_file, sheet='Results', skip=46, n_max=96)
    colnames(result) <- unname(sapply(colnames(result), title_to_snake_case))
    result[['ct']] <- unname(sapply(result[['ct']], function(x) gsub("Undetermined", NA, x)))
    result[['ct']] <- as.numeric(result[['ct']])
    result <- merge(plate_setup, result[c('well', 'ct', 'ct_threshold', 'amp_score', 'cq_conf', 'tm1')],
        by='well', all.x=TRUE, all.y=FALSE
    )
    result_dfs[[sample_name]] <- result

    # amp_data
    amp_data <- read_excel(qpcr_file, sheet='Amplification Data', skip=46)
    colnames(amp_data) <- unname(sapply(colnames(amp_data), title_to_snake_case))
    amp_data <- merge(plate_setup, amp_data[c('well', 'cycle', 'rn', 'delta_rn')],
        by='well', all.x=TRUE, all.y=FALSE
    )
    amp_data_dfs[[sample_name]] <- amp_data

}
result_df <- do.call(rbind, as.list(result_dfs))
amp_data_df <- do.call(rbind, as.list(amp_data_dfs))


# ----------------------------------------------------------------------
# QC

for (input_dir in input_dirs) {

    result_subset <- result_df[(result_df['sample_name']==input_dir), ]

    # plot heatmap
    plate <- df_to_plate(result_subset, value='ct')
    fig <- plot_heatmap(plate,
        show_xlabel=FALSE,
        show_ylabel=FALSE,
        title=paste('CT for', input_dir), annotations=TRUE, digits=2)

    savefig(file.path(wd, opt[['output']], input_dir, paste0('heatmap-', input_dir, '.png')),
            dpi=400,
            troubleshooting=troubleshooting)
}


# ----------------------------------------------------------------------
# Plot Amplification Curves

amp_data_df <- within(amp_data_df,
    sample_id <- paste(sample_name, tissue, gene, well_position, sep=', ')
)
amp_data_df <- amp_data_df[(amp_data_df[['gene']]!=''), ]

for (input_dir in input_dirs) {
    sample_name <- input_dir

    ct_threshold <- result_df[(result_df['sample_name']==sample_name), 'ct_threshold'][[1]]

    ggplot(amp_data_df[
               (amp_data_df['sample_name']==sample_name) &
               (amp_data_df['delta_rn'] > 0), ],
       aes(x=.data[['cycle']], y=.data[['delta_rn']],
           group=sample_id,
           colour=gene)) +
    geom_line(alpha=0.7, size=0.5) +
    geom_hline(yintercept=ct_threshold, colour='red') +
    labs(x='Cycle Number', y='Delta Rn', title=paste("Amplification Curve for", sample_name)) +
    scale_y_continuous(trans='log10', labels = function(x) round(x, 5), limits = c(0.00001, 15))

    savefig(file.path(wd, opt[['output']], input_dir, paste0('delta_rn-', input_dir, '.png')),
            height=1000, width=1600, dpi=400,
            troubleshooting=troubleshooting)

}


# ----------------------------------------------------------------------
# Compute dCT Table

log_print(paste(Sys.time(), 'Computing ddct...'))

index_cols <- c('sample_name', 'tissue', 'gene')
ct_long <- result_df %>%
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

# ----------------------------------------------------------------------
# Plot Fold Change

log_print(paste(Sys.time(), 'Plotting fold change...'))

control_genes <- c('Actin', 'Hprt')

for (control_gene_name in control_genes) {

    fold_change_col <- paste0('fold_change_dnase1l1_', tolower(control_gene_name))

    ct_wide <- ct_wide[order(ct_wide[[fold_change_col]], decreasing = TRUE),]
    ct_wide <- reset_index(ct_wide, drop=TRUE)

    # manual filters
    ct_wide <- ct_wide[(ct_wide['stdev_ct_dnase1l1'] <= 2), ]
    ct_wide <- ct_wide[(ct_wide[fold_change_col] <= 1), ]
    ct_wide <- ct_wide[(ct_wide['tissue'] != 'pc'), ]

    # barsize
    barsize <- ct_wide %>%
        group_by(tissue) %>%
        summarize(min_fold_change=min(!!sym(fold_change_col), na.rm=TRUE),
                  max_fold_change=max(!!sym(fold_change_col), na.rm=TRUE))
    tmp <- merge(ct_wide, barsize,
        by='tissue', all.x=TRUE, all.y=FALSE 
    )

    tmp <- tmp[!is.na(tmp['sample_name']), ]

    ggplot(
        tmp,
        aes(x=reorder(.data[['tissue']], .data[[fold_change_col]], decreasing=TRUE),
            y=.data[[fold_change_col]]), na.rm=TRUE) +
        # see: https://stackoverflow.com/questions/32642856/how-do-i-use-geom-rect-with-discrete-axis-values
        geom_rect(
            aes(xmin = stage(
                    reorder(.data[['tissue']], .data[[fold_change_col]], decreasing=TRUE),
                    after_scale = xmin-0.45),
                xmax = stage(
                    reorder(.data[['tissue']], .data[[fold_change_col]], decreasing=TRUE),
                    after_scale = xmax+0.45),
                ymin = .data[['min_fold_change']],
                ymax = .data[['max_fold_change']]),
                fill="#7f7f7f",
                stat='identity', inherit.aes=TRUE) +
        geom_jitter(aes(colour=.data[['sample_name']])) +
        labs(x='Tissue', y='Fold Change', title=paste("Relative Expression of Dnase1l1 vs.", control_gene_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_y_continuous(trans='log10', labels = function(x) round(x, 5))

    savefig(file.path(wd, opt[['output']], paste0(fold_change_col, '.png')),
            width=1600, dpi=400,
            troubleshooting=troubleshooting)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
