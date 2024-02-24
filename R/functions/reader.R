import::here(readxl, 'read_excel')
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', 'read_excel_or_csv', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'move_list_item_to_start', .character_only=TRUE)

## Functions
## read_qpcr


#' Read QPCR
#' 
#' @description
#' See the README for the input directory structure. (Note: TODO)
#'
read_qpcr <- function(input) {

    input_dirs <- items_in_a_not_b(list.dirs(input, full.names=FALSE), '')

    result_dfs <- new.env()
    amp_data_dfs <- new.env()

    for (input_dir in input_dirs) {

        sample_id <- input_dir  # may not be true in future experiments
        qpcr_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='xls')[[1]]


        # plate setup file
        plate_setup_file <- list_files(file.path(wd, opt[['input']], input_dir), ext='csv')
        plate_setup <- read_excel_or_csv(plate_setup_file)  # change to read_csv
        plate_setup['sample_id'] <- sample_id
        plate_setup <- plate_setup[, 
            move_list_item_to_start(colnames(plate_setup), 'sample_id')
        ]

        # result
        result <- read_excel(qpcr_file, sheet='Results', skip=46, n_max=96)
        colnames(result) <- unname(sapply(colnames(result), title_to_snake_case))
        result[['ct']] <- unname(sapply(result[['ct']], function(x) gsub("Undetermined", NA, x)))
        result[['ct']] <- as.numeric(result[['ct']])
        result <- merge(plate_setup, result[c('well', 'ct', 'ct_threshold', 'amp_score', 'cq_conf', 'tm1')],
            by='well', all.x=TRUE, all.y=FALSE
        )
        result_dfs[[sample_id]] <- result

        # amp_data
        amp_data <- read_excel(qpcr_file, sheet='Amplification Data', skip=46)
        colnames(amp_data) <- unname(sapply(colnames(amp_data), title_to_snake_case))
        amp_data <- merge(plate_setup, amp_data[c('well', 'cycle', 'rn', 'delta_rn')],
            by='well', all.x=TRUE, all.y=FALSE
        )
        amp_data_dfs[[sample_id]] <- amp_data

    }
    result_df <- do.call(rbind, as.list(result_dfs))
    amp_data_df <- do.call(rbind, as.list(amp_data_dfs))

    return( list(result_df, amp_data_df) )
}
