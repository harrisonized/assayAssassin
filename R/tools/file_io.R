import::here(readxl, 'read_excel')
import::here(ggplot2, 'ggsave', 'last_plot')
import::here(grid, 'grid.newpage', 'grid.draw')

## Functions
## list_files
## read_excel_or_csv
## read_tsv_from_text
## savefig


#' List all files with a specific extension
#' 
#' @description
#' This is a thin wrapper around [list.files()].
#' 
#' @references
#' \href{https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching}{StackOverflow post}
#' 
#' @export
list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: 
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


#' Switch case to read excel or csv based on the extension
#'
#' @description Mainly used to simplify scripts
#' 
#' @export
read_excel_or_csv <- function(filepath) {
    ext=tools::file_ext(filepath)
    if (ext == 'xlsx') {
        df <- read_excel(filepath)
    } else if (ext == 'csv') {
        df <- read.csv(filepath, header=TRUE, check.names=FALSE)
    } else {
        log_print(paste(Sys.time(), 'Please enter a xlsx or csv file.'))
        stop()
    }
    return(df)
}


#' Parse tab-delimited text line-by-line
#'
#' @description
#' This is useful for parsing data exported from plate readers, as it allows you to retrieve a table internal to a spreadsheet.
#' Note that there is a Mac-specific bug that prevents readLines from reading files encoded in UTF-16.
#' 
#' @usage
#' # for 96-well plates:
#' df <- read_csv_from_text(
#'   file_path,
#'   skiprows=3, nrows=8,
#'   skipcols=2, ncols=12,
#'   index=LETTERS[1:8],
#'   columns=seq(1, 12),
#'   numeric=TRUE
#' )
#' 
#' @export
read_tsv_from_text <- function(
    file_path,
    encoding='UTF-8', sep='\t',
    skiprows=0, nrows=NULL,
    skipcols=0, ncols=NULL,
    index=NULL,
    columns=NULL,
    numeric=FALSE
) {

    con = file(file_path, encoding=encoding)
    withCallingHandlers({
        rawData <- readLines(con)
    }, warning = function(w) {
        # print(conditionMessage(w))
        if (startsWith(conditionMessage(w), "incomplete final line")) {
            invokeRestart( "muffleWarning" )
        }
    })
    close(con)
    
    # autodetermine ranges if not specified
    if(is.null(nrows)) {
        nrows <- length(rawData)
    }
    if(is.null(ncols)) {
        rowArr = unlist(strsplit(rawData[1+skiprows], split='\t'))
        ncols = length(rowArr)-skipcols
    }
    
    # instantiate empty dataframe and append row-by-row
    df <- data.frame(matrix(ncol=ncols, nrow=0))
    for (row in rawData[(1+skiprows):(nrows+skiprows)]) {
        rowArr <- unlist(strsplit(row, split='\t'))
        df[nrow(df) + 1,] = rowArr[(1+skipcols):(ncols+skipcols)]
    }
    
    # rename columns
    colnames(df) <- columns
    rownames(df) <- index

    if(numeric) {
        df[] <- lapply(df, function(x) as.numeric(as.character(x)))
    }

    return(df)
}


#' Save Figure
#' 
#' @description Switch case to reduce the number of lines in the main script
#' 
#' @export
savefig <- function(
    filepath,
    fig=NULL,
    height=800, width=1200, dpi=300, units="px", scaling=0.5,
    makedir=TRUE,
    troubleshooting=FALSE,
    lib='ggplot',  # choose: ggplot, grid
    default_ext = '.png'
) {
    if (!troubleshooting) {

        # make directory
        dirpath <- dirname(filepath)
        if (makedir && !dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }

        # add file extension if not included
        if (tools::file_ext(filepath)=='') {
            filepath <- paste0(filepath, default_ext)
        }

        if (lib=='ggplot') {
            if (!inherits(fig, "ggplot")) {
                fig <- last_plot()
            }
            ggsave(
                filepath,
                plot=fig,
                height=height, width=width, dpi=dpi, units=units, scaling=scaling
            )
        } else if (lib=='grid') {
            png(filepath,
                height=height, width=width, res=dpi, units=units
            )
            grid.newpage()
            grid.draw(fig$gtable)
            dev.off()
        } else {
            warning(paste0("lib='", lib, "' not found"))
        }
    }
}
