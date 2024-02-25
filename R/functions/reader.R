import::here(readxl, 'read_excel')
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', 'read_excel_or_csv', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'title_to_snake_case', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', 'move_list_item_to_start', .character_only=TRUE)

## Functions
## read_qpcr


#' Read QPCR
#' 
#' @description
#' TODO: Add input directory structure to README.
#'
read_qpcr <- function(input_path) {

    all_results <- new.env()
    all_amp_data <- new.env()

    plate_ids <- items_in_a_not_b(list.dirs(input_path, full.names=FALSE), '')
    for (plate_id in plate_ids) {

        # TODO: try/catch for broken files
        qpcr_file <- list_files(file.path(input_path, plate_id), ext='xls')[[1]]


        # plate setup file
        plate_setup <- read.csv(
            file.path(input_path, plate_id, 'plate-setup.csv'),
            header=TRUE, check.names=FALSE
        )
        plate_setup[['sample_id']] <- as.character(plate_setup[['sample_id']])
        plate_setup['plate_id'] <- plate_id
        plate_setup <- plate_setup[, 
            move_list_item_to_start(colnames(plate_setup), 'plate_id')
        ]
        

        # results
        results <- read_excel(qpcr_file, sheet='Results', skip=46, n_max=96)
        colnames(results) <- unname(sapply(colnames(results), title_to_snake_case))
        results[['ct']] <- unname(sapply(results[['ct']], function(x) gsub("Undetermined", NA, x)))
        results[['ct']] <- as.numeric(results[['ct']])
        results <- merge(
            plate_setup,
            results[c('well', 'ct', 'ct_threshold', 'amp_score', 'cq_conf', 'tm1')],
            by='well', all.x=TRUE, all.y=FALSE
        )
        all_results[[plate_id]] <- results


        # amp_data
        amp_data <- read_excel(qpcr_file, sheet='Amplification Data', skip=46)
        colnames(amp_data) <- unname(sapply(colnames(amp_data), title_to_snake_case))
        amp_data <- merge(
            plate_setup,
            amp_data[c('well', 'cycle', 'rn', 'delta_rn')],
            by='well', all.x=TRUE, all.y=FALSE
        )
        all_amp_data[[plate_id]] <- amp_data
        
    }

    results <- do.call(rbind, as.list(all_results))
    amp_data <- do.call(rbind, as.list(all_amp_data))

    return( list(results, amp_data, plate_ids) )
}
