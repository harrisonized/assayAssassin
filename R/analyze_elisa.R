## Basic ELISA analysis

wd = dirname(this.path::here())  # wd = '~/github/R/assayAssassin'
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('ggplot2'))
library('zeallot')  # %<-%
library('optparse')
suppressPackageStartupMessages(library('logr'))
import::from(magrittr, '%>%')
import::from(plyr, 'rbind.fill')

import::from(file.path(wd, 'R', 'functions', 'reader.R'),
    'read_elisa', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', 'read_tsv_from_text', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'smelt', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'save_fig', 'plot_violin', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/elisa/input',
                metavar='data/elisa/input',
                type="character",help="path/to/input/dir"),
   
    make_option(c("-o", "--output-dir"), default="data/elisa/output",
                metavar="data/elisa/output", type="character",
                help="set the output directory for the data"),

    make_option(c("-f", "--figures-dir"), default="figures/elisa/output",
                metavar="figures/elisa/output", type="character",
                help="set the output directory for the figures"),

    make_option(c("-d", "--drop-plates"), default="",
                metavar="", type="character",
                help="semicolon(;)-separated list of unwanted plate_ids"),

    make_option(c("-s", "--save-html"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="save html files in addition to PNG"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

drop_plates <- strsplit(opt[['drop-plates']], ';')[[1]]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("analyze_elisa-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read data

log_print(paste(Sys.time(), 'Reading data...'))
input_path <- file.path(wd, opt[['input-dir']])
df <- read_elisa(input_path)


# ----------------------------------------------------------------------
# Plot

log_print(paste(Sys.time(), 'Plotting...'))
fig <- plot_violin(
    df[(df[['develop_time']]=='45 min'), ],
    x='treatment_length', y='abs', group_by='genotype',
    ylabel='Absorbance',
    hover_data=c(
        'plate_id','sample_id', 'genotype', 'treatment', 'treatment_length', 'abs'
    )
)

# save
if (!troubleshooting) {
    save_fig(
        fig=fig,
        dirpath=file.path(wd, opt[['figures-dir']]),
        filename='violin-abs_by_age',
        save_html=opt[['save-html']]
    )
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
