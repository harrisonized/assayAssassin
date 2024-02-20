import::here(readxl, 'read_excel')
import::here(ggplot2, 'ggsave', 'last_plot')
import::here(grid, 'grid.newpage', 'grid.draw')

## Functions
## list_files
## read_excel_or_csv
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


#' Save Figure
#' 
#' @description Switch case to reduce the number of lines in the main script
#' 
#' @export
savefig <- function(
    filepath,
    fig=NULL,
    height=800, width=1200, dpi=300, units="px", scaling=0.5,
    makedir=FALSE,
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
