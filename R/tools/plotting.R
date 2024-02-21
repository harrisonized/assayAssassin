import::here(rlang, 'sym')
import::here(ggplot2,
    'ggplot', 'aes', 'aes_string', 'theme', 'labs', 
    'geom_bar', 'geom_tile', 'geom_text', 'coord_fixed',
    'guide_axis', 'scale_x_discrete', 'scale_y_discrete', 'scale_fill_gradient',
    'element_text', 'element_blank')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'smelt', .character_only=TRUE)


## Functions
## plot_bar

#' Plot Bar
#' 
#' @description Plot a simple bar plot.
#' 
#' @export
plot_bar <- function(
    df,
    x='cell_type',
    y='value',
    group.by=NULL,  # gene of interest
    fill=NULL,  # steelblue, overwrites group.by
    xlabel=NULL,
    ylabel="Number of Genes",
    title="Number of Cells",
    xaxis_angle=45,
    legend_position='bottom',
    sort=TRUE
) {

    # color
    if (is.null(group.by)) { color <- NULL } else { color <- sym(group.by) }
    if (is.null(fill)) {
        bar <- geom_bar(stat="identity")
    } else {
        bar <- geom_bar(stat="identity", fill=fill)
    }

    # plot base
    if (sort) {
        base_plot <- ggplot(data=df,
            aes(x=reorder(.data[[x]], .data[[y]], decreasing=TRUE),
                y=.data[[y]], fill=!!color ))
    } else {
        base_plot <- ggplot(data=df,
            aes(x=.data[[x]], y=.data[[y]], fill=!!color ))
    }

    fig <- base_plot +
        bar +
        labs(x=xlabel, y=ylabel, title=title) +
        scale_x_discrete( guide=guide_axis(angle = xaxis_angle) ) +
        theme(legend.position=legend_position)

    return(fig)
}


#' Plot a heatmap
#'
#' @usage
#' plot_heatmap(
#'   df,
#'   title="Raw Data",
#'   show_xlabel=TRUE,
#'   show_ylabel=TRUE,
#'   annotations=TRUE,
#'   scientific_notation=FALSE,
#'   digits=0
#' )
#' 
#' @section Vignette:
#' See `vignettes/plot_heatmap.Rmd`
#' 
#' @export
plot_heatmap <- function(
   df,
   x='col',
   y='row',
   fill='val',
   title=NULL,
   show_xlabel=TRUE,
   show_ylabel=TRUE,
   annotations=FALSE,
   scientific_notation=FALSE,
   digits=1
) {
    
    tab <- smelt(rev_df(df))  # reshape

    # axis labels
    if (show_xlabel) {
        xlabel = element_text()
    } else {
        xlabel = element_blank()
    }
    if (show_ylabel) {
        ylabel = element_text()
    } else {
        ylabel = element_blank()
    }

    # annotations
    if (annotations) {
        label = 'label'
        if (scientific_notation) {
            tab['label'] = lapply(
                tab['val'], 
                function(x) formatC(x, format='e', digits=2)
            )
        } else {
            tab['label'] = lapply(
                tab['val'], 
                function(x) round(x, digits)
            )
        }
    } else {
        label = NULL
    }

    # plot
    fig <- ggplot(tab, aes_string(x=x, y=y, fill=fill)) +
        geom_tile(color="white", lwd=0.3, linetype=1) +
        coord_fixed(expand=TRUE) +
        scale_y_discrete(limits=rev) +
        labs(title = title) +
        theme(plot.title = element_text(size = 10),
                       axis.title.x = xlabel,
                       axis.title.y = ylabel) +
        scale_fill_gradient(low="#FFF8F8", high="#A50026") +
        if (annotations) {
            geom_text(
                aes_string(label=label),
                color = 'black',
                size = 2
            )
        }

    return (fig)
}

