## Basic qPCR analysis

wd = dirname(this.path::here())  # wd = '~/github/R/leige-waffle'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('ggplot2'))
library('zeallot')  # %<-%
library('optparse')
library('logr')
import::from(magrittr, '%>%')

import::from(file.path(wd, 'R', 'functions', 'reader.R'),
    'read_qpcr', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'draw_plots.R'),
    'draw_ct_heatmaps', 'draw_amp_curves', 'draw_fold_changes', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'computations.R'),
    'ct_thresholds_from_results', 'dct_table_from_results', .character_only=TRUE)


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

    make_option(c("-c", "--min-cq-conf"), default=0.75,
                metavar="0.75", type="numeric",
                help="set the cq_conf cutoff level"),

    make_option(c("-d", "--drop-plates"), default="12176",
                metavar="12176", type="character",
                help="semicolon(;)-separated list of unwanted plate_ids"),

    make_option(c("-s", "--sample-genes"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="semicolon(;)-separated list of sample genes"),

    make_option(c("-g", "--control-genes"), default="Actin;Hprt",
                metavar="Actin;Hprt", type="character",
                help="semicolon(;)-separated list of control genes"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

drop_plates <- strsplit(opt[['drop-plates']], ';')[[1]]
sample_genes <- strsplit(opt[['sample-genes']], ';')[[1]]
control_genes <- strsplit(opt[['control-genes']], ';')[[1]]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("plot_qpcr-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read data

log_print(paste(Sys.time(), 'Reading data...'))

c(results, amp_data, plate_ids, metadata_cols) %<-% read_qpcr(
    file.path(wd, opt[['input-dir']])
)
ct_thresholds <- ct_thresholds_from_results(results)


# ----------------------------------------------------------------------
# Draw Ct Heatmaps

log_print(paste(Sys.time(), 'Drawing CT heatmaps...'))

draw_ct_heatmaps(
    results[(results[['cq_conf']] > opt[['min-cq-conf']]), ],
    dirpath=file.path(wd, opt[['figures-dir']]),
    troubleshooting=troubleshooting,
    showfig=troubleshooting
)

# ----------------------------------------------------------------------
# Draw Amplification Curves

log_print(paste(Sys.time(), 'Drawing amp curves...'))

amp_data <- within(amp_data,
    row_id <- paste(sample_id, tissue, gene, well_position, sep=', ')
)
amp_data <- merge(
    amp_data,
    results[, c('sample_id', 'well', 'cq_conf')],
    by=c('sample_id', 'well'), all.x=TRUE, all.y=FALSE 
)

draw_amp_curves(
    amp_data[(amp_data[['cq_conf']] > opt[['min-cq-conf']]), ],
    dirpath=file.path(wd, opt[['figures-dir']]),
    ct_thresholds=ct_thresholds,
    metadata_cols=c(metadata_cols, 'cq_conf'),
    troubleshooting=troubleshooting,
    showfig=troubleshooting
)

# ----------------------------------------------------------------------
# Plot Fold Changes

log_print(paste(Sys.time(), 'Plotting fold change...'))

dct_table <- dct_table_from_results(
    results[(results[['ignore']]==0),],
    index_cols=c('plate_id', 'sample_id', 'tissue', 'gene'),
    sample_genes=sample_genes,
    control_genes=control_genes
)

# save dct_table
if (!troubleshooting) {
    output_dir <- file.path(wd, opt[['output-dir']])
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }
    filepath = file.path(output_dir, 'dct_table.csv')
    write.table(dct_table, file = filepath, row.names = FALSE, sep = ',' )
}

# filter low quality
for (sample_gene in sample_genes) {
    dct_table <- dct_table[
        (dct_table[paste0('stdev_ct_', tolower(sample_gene))] <= 2),
    ]  # sLN
}

# exclude plates
dct_table <- dct_table[
    !(dct_table[['plate_id']] %in% drop_plates) &
    (dct_table[['plate_id']] != drop_plates),
]
dct_table <- dct_table[!is.na(dct_table['sample_id']), ]

draw_fold_changes(
    dct_table,
    sample_genes=sample_genes,
    control_genes=control_genes,
    dirpath=file.path(wd, opt[['figures-dir']]),
    troubleshooting=troubleshooting,
    showfig=troubleshooting
)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
