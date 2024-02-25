## Basic qPCR analysis

wd = dirname(this.path::here())  # wd = '~/github/R/leige-waffle'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('ggplot2'))
library('zeallot')  # %<-%
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(tidyr, 'pivot_wider')

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'df_to_plate', 'set_index', 'reset_index', 'rev_df', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_bar', 'plot_heatmap', 'plot_fold_change', .character_only=TRUE)

import::from(file.path(wd, 'R', 'functions', 'reader.R'),
    'read_qpcr', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'draw_plots.R'),
    'draw_ct_heatmaps', 'draw_amp_curves', 'draw_fold_changes', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_dct_table', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/input',
                metavar='data/input',
                type="character",help="path/to/input/dir"),
   
    make_option(c("-o", "--output-dir"), default="data/output",
                metavar="data/output", type="character",
                help="set the output directory for the data"),

    make_option(c("-f", "--figures-dir"), default="figures/output",
                metavar="figures/output", type="character",
                help="set the output directory for the figures"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("plot_qpcr-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read data

log_print(paste(Sys.time(), 'Reading data...'))

c(results, amp_data, plate_ids) %<-% read_qpcr(
    file.path(wd, opt[['input-dir']])
)

# qc
draw_ct_heatmaps(
    results,
    plate_ids,
    dirpath=file.path(wd, opt[['figures-dir']]),
    troubleshooting=troubleshooting,
    showfig=troubleshooting
)

# ----------------------------------------------------------------------
# Plot Amplification Curves

amp_data <- within(amp_data,
    row_id <- paste(sample_id, tissue, gene, well_position, sep=', ')
)
amp_data <- merge(
    amp_data,
    results[, c('well', 'sample_id', 'cq_conf')],
    by=c('well', 'sample_id'), all.x=TRUE, all.y=FALSE 
)

draw_amp_curves(
    amp_data,
    plate_ids,
    dirpath=file.path(wd, opt[['figures-dir']]),
    troubleshooting=troubleshooting,
    showfig=troubleshooting
)

# ----------------------------------------------------------------------
# Plot Fold Change

log_print(paste(Sys.time(), 'Plotting fold change...'))

ct_wide <- compute_dct_table(results)
ct_wide <- ct_wide[(ct_wide['sample_id'] != '12176'), ]

# manual filters
ct_wide <- ct_wide[(ct_wide['stdev_ct_dnase1l1'] <= 2), ]
ct_wide <- ct_wide[(ct_wide['tissue'] != 'pc'), ]

draw_fold_changes(
    ct_wide,
    dirpath=file.path(wd, opt[['figures-dir']]),
    troubleshooting=troubleshooting,
    showfig=troubleshooting
)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
