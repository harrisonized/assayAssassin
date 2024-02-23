## Functions
## items_in_a_not_b


#' Return items unique to vector a
#' 
#' @description
#' Return all the items found in a and not b. Useful for filtering dataframe columns.
#' 
#' @param a list or vector
#' @param b list or vector
#' @return Returns the filtered list or vector (same type as a)
#' 
#' @references
#' \href{https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list}{Stack Overflow}
#' 
#' @examples
#' items_in_a_not_b(c('a', 'b', 'c', '1', '2'), c('1', '2'))
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' Move List Items to Start
#' 
move_list_item_to_start <- function(input_list, item) {
    output_list = c(item, items_in_a_not_b(input_list, item))
    return(output_list)
}
