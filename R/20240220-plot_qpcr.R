## Basic qPCR analysis

wd = dirname(this.path::here())  # wd = '~/github/R/qpcr-analysis'
suppressPackageStartupMessages(library('dplyr'))
# suppressPackageStartupMessages(library('ggplot2'))
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(readxl, 'read_excel')
import::from(tidyr, 'pivot_wider')
import::from(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',  'geom_bar', 'geom_jitter',
    'guide_axis', 'scale_x_discrete', 'scale_y_log10', 'element_text')

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'df_to_plate', 'set_index', 'reset_index', 'rev_df', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', 'read_excel_or_csv', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
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
   
    make_option(c("-o", "--output"), default="figures/20240220-second-pass",
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
log <- log_open(paste0("20240215-plot_qpcr-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read data

log_print(paste(Sys.time(), 'Reading data...'))

input_dirs <- items_in_a_not_b(
    list.dirs(file.path(wd, opt[['input']]), full.names=FALSE), ''
)  # sample_names

dfs <- new.env()
for (input_dir in input_dirs) {

    # qpcr data
    qpcr_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='xls')
    qpcr_data <- read_excel(qpcr_file, sheet='Results', skip=46, n_max=96)
    colnames(qpcr_data) <- unname(sapply(colnames(qpcr_data), title_to_snake_case))
    qpcr_data[['ct']] <- unname(sapply(qpcr_data[['ct']], function(x) gsub("Undetermined", NA, x)))
    qpcr_data[['ct']] <- as.numeric(qpcr_data[['ct']])

    # plate setup file
    plate_setup_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='csv')
    plate_setup <- read_excel_or_csv(plate_setup_file)
    plate_setup['sample_name'] <- input_dir

    df <- merge(
        plate_setup[, c('sample_name', 'well', 'well_position', 'tissue', 'gene', 'task')],
        qpcr_data[c('well', 'ct')],
        by='well',
        all.x=TRUE,
        all.y=FALSE
    )

    dfs[[input_dir]] <- df
}
dfs <- as.list(dfs)
df <- do.call(rbind, dfs)

# ----------------------------------------------------------------------
# QC

for (input_dir in input_dirs) {

    subset <- df[(df['sample_name']==input_dir), ]

    # plot heatmap
    plate <- df_to_plate(subset, value='ct')
    fig <- plot_heatmap(plate, annotations=TRUE, digits=2)
    savefig(file.path(wd, opt[['output']], input_dir, paste0('heatmap-', input_dir, '.png')),
            dpi=400,
            troubleshooting=troubleshooting)
}


# ----------------------------------------------------------------------
# Compute dCT Table

log_print(paste(Sys.time(), 'Computing ddct...'))

index_cols <- c('sample_name', 'tissue', 'gene')
ct_long <- df %>%
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

ct_wide <- ct_wide[order(ct_wide[['fold_change_dnase1l1_actin']], decreasing = TRUE),]
ct_wide <- reset_index(ct_wide)
ct_wide <- ct_wide[!is.na(ct_wide['sample_name']), ]


# manually filter
ct_wide <- ct_wide[(ct_wide['fold_change_dnase1l1_actin'] <= 1), ]
ct_wide <- ct_wide[(ct_wide['tissue'] != 'pc'), ]  # low quality


# ----------------------------------------------------------------------
# Plot

log_print(paste(Sys.time(), 'Plotting...'))

ggplot(
    ct_wide,
    aes(x=reorder(.data[['tissue']], .data[['fold_change_dnase1l1_actin']], decreasing=TRUE),
        y=.data[['fold_change_dnase1l1_actin']],
        fill=.data[['tissue']])) +
    geom_bar(stat="identity", position = "dodge") +
    geom_jitter() +
    labs(x='Tissue', y='Fold Change', title="Relative Expression of Dnase1l1 vs. Actin") +
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none') +
    scale_y_log10()

savefig(file.path(wd, opt[['output']], paste0('relative_expression.png')),
        dpi=400,
        troubleshooting=troubleshooting)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
