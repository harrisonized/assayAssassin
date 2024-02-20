import::here(rlang, 'sym')
import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'labs',  'geom_bar',
    'guide_axis', 'scale_x_discrete')

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
