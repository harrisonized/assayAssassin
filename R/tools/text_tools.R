## Functions
## title_to_snake_case

#' Standardize Space-separated Titles
#' 
#' @description Converts "Column Title" to column_title
#' 
#' @examples
#' title_to_snake_case('Column Title')
#' 
#' @export
title_to_snake_case <- function(text) {
    return(tolower(
        paste(
            unlist(strsplit(text, '[ ]')), collapse='_')
        )
    )
}
