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
    'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', 'read_excel_or_csv', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_bar', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input"), default='data/input',
                metavar='data/input',
                type="character",help="path/to/input/dir"),
   
    make_option(c("-o", "--output"), default="figures/20240215-first-pass",
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

input_dirs <- items_in_a_not_b(list.dirs(file.path(wd, opt[['input']]), full.names=FALSE), '')

dfs <- new.env()
for (input_dir in input_dirs) {

    # qpcr data
    qpcr_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='xls')
    qpcr_data <- read_excel(qpcr_file, sheet='Results', skip=46, n_max=96)
    colnames(qpcr_data) <- unname(sapply(colnames(qpcr_data), title_to_snake_case))
    qpcr_data[['ct']] <- as.numeric(qpcr_data[['ct']])

    # plate setup file
    plate_setup_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='csv')
    plate_setup <- read_excel_or_csv(plate_setup_file)
    plate_setup['mouse_id'] <- input_dir

    df <- merge(
        plate_setup[, c('mouse_id', 'well', 'well_position', 'tissue', 'gene', 'task')],
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
# Compute CT Table

log_print(paste(Sys.time(), 'Computing ddct...'))

ct_table <- df %>%
    group_by(mouse_id, tissue, gene) %>%
    summarize(mean_ct=mean(ct, na.rm=TRUE))
ct_table <- ct_table[(ct_table[['tissue']]!=''),]

dct_table <- pivot_wider(ct_table, names_from=c(gene), values_from=mean_ct)
dct_table['housekeeping'] <- rowMeans(dct_table[c('actin', 'hprt')], na.rm = TRUE)

dct_table <- merge(
    dct_table,
    dct_table %>%
        group_by(mouse_id) %>%
        summarize(
            mean_dnase1l1=mean(dnase1l1, na.rm=TRUE),
            mean_housekeeping=mean(housekeeping, na.rm=TRUE)
        ),
    by='mouse_id',
    all.x=TRUE,
    all.y=FALSE
)

dct_table[['ddct']] <- (
    (dct_table[['mean_housekeeping']] - dct_table[['housekeeping']]) -
    (dct_table[['mean_dnase1l1']] - dct_table[['dnase1l1']])
    
)
dct_table[['fold_change']] <- 2^(-dct_table[['ddct']])
dct_table <- dct_table[(dct_table['fold_change'] <= 100), ]

dct_table <- dct_table[order(dct_table[['fold_change']], decreasing = TRUE),]
dct_table <- reset_index(dct_table)
dct_table <- dct_table[!is.na(dct_table['mouse_id']), ]


# manually filter
dct_table <- dct_table[(dct_table['tissue'] != 'pc'), ]


# ----------------------------------------------------------------------
# Plot

log_print(paste(Sys.time(), 'Plotting...'))

fig <- plot_bar(dct_table,
    x='tissue',
    y='fold_change'
)

ggplot(
    dct_table,
    aes(x=reorder(.data[['tissue']], .data[['fold_change']], decreasing=TRUE),
        y=.data[['fold_change']],
        fill=.data[['tissue']])) +
    geom_bar(stat="identity", position = "dodge") +
    geom_jitter() +
    labs(x='Tissue', y='Fold Change', title="Relative Expression") +
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none') +
    scale_y_log10()

savefig(file.path(wd, opt[['output']], paste0('relative_expression.png')),
        dpi=400,
        troubleshooting=troubleshooting)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
