import::here(rlang, 'sym')
import::here(dplyr, 'group_by', 'summarize')
import::here(tidyr, 'pivot_wider')

## Functions
## set_index
## reset_index
## rev_df
## smelt
## df_to_plate


#' Set Index
#' 
#' @description
#' Uses the values in a column to set an index. Needs to be unique.
#' Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.set_index.html}{set_index}.
#' 
#' @export
set_index <- function(df, colname, drop=TRUE) {

    rownames(df) <- df[[colname]]
    if (drop) {
        df <- df[, !names(df)==colname]
    }
    return(df)
}


#' Reset Index
#' 
#' @description
#' Moves the values in the index to a column. Resets the index to the default integer index.
#' Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.reset_index.html}{reset_index}.
#' 
#' @param df a dataframe
#' @param index_name select a new name for the index column
#' @param drop if TRUE, does not copy the index values to the new column
#' @return Returns a dataframe with index renamed and new integer index 
#' 
#' @examples
#' reset_index(mtcars, index_name='model')
#' 
#' @references
#' \href{https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column}{StackOverflow post}
#' 
#' @export
reset_index <- function(df, index_name='index', drop=FALSE) {
    if (!drop) {
        df <- cbind(index = rownames(df), df)
        colnames(df)[colnames(df) == "index"] = index_name
    }
    rownames(df) <- 1:nrow(df)
    return (df)
}


#' Reverse the order of a dataframe
#' 
#' @param df a dataframe
#' @param how 'row' or 'col'
#' @return Returns the reversed dataframe.
#' 
#' @examples
#' rev_df(mtcars, how='row')
#' rev_df(mtcars, how='col')
#' 
#' @export
rev_df <- function(df, how='row') {
    if (how == 'row') {
        return(df[dim(df)[1]:1,])
    } else if (how == 'col') {
        return(rev(df))
    } else {
        stop("Choose how='row' or how='col'")
    }
}


#' Special Melt
#' 
#' @description
#' Convert a dataframe from a table to long format
#' This is an alternative to melt that doesn't throw errors
#' 
#' @examples
#' smelt(mtcars[, c('cyl', 'mpg')])
#' 
#' @references
#' See: \href{https://stackoverflow.com/questions/28355377/how-to-add-index-of-a-list-item-after-melt-in-r}{Stack Overflow link}
#' 
#' @export
smelt <- function(
   df,
   rowname='row',
   colname='col',
   valname='val'
) {
   melted <- transform(stack(setNames(df, colnames(df))), id=rownames(df))
   colnames(melted) <- c(valname, colname, rowname)
   return(rev(melted))
}


#' Dataframe to Plate
#' 
#' @description
#' Reshape a column dataframe to a plate format
#' 
#' @export
df_to_plate <- function(df, well_id='well_position', value='ct') {

    tmp <- df[, c(well_id, value)] %>%
        group_by(!!sym(well_id)) %>%
        summarize(!!(paste0('mean_', value)) := mean(!!sym(value), na.rm=TRUE))
    tmp[['row']] = lapply(tmp[well_id], function(x) gsub("[^a-zA-Z]", "", x))[[1]]
    tmp[['col']] = as.numeric(lapply(tmp[well_id], function(x) gsub("[^0-9]", "", x))[[1]])
    tmp <- tmp[order(tmp$row, tmp$col),]

    plate <- as.data.frame(pivot_wider(tmp, id_cols='row', names_from='col', values_from=paste0('mean_', value)))
    plate <- set_index(plate, colname='row')

    return(plate)
}
