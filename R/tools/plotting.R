import::here(rlang, 'sym')
import::here(ggplot2,
    'ggplot', 'aes', 'aes_string', 'theme', 'labs',
    'geom_bar', 'geom_line', 'geom_hline', 'geom_tile', 'geom_rect', 'geom_jitter', 
    'geom_text', 'coord_fixed', 'guide_axis',
    'scale_x_discrete', 'scale_y_continuous', 'scale_y_discrete', 'scale_fill_gradient',
    'element_text', 'element_blank')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rev_df', 'smelt', .character_only=TRUE)

## Functions
## plot_bar
## plot_heatmap
## plot_amp_curves
## plot_fold_change


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
    xlabel=NULL,
    ylabel=NULL,
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
        labs(x=xlabel, y=ylabel, title = title) +
        theme(plot.title = element_text(size = 10),
                       axis.title.x = xlabel,
                       axis.title.y = ylabel) +
        scale_fill_gradient(low="#FFF8F8", high="#A50026") +
        if (annotations) {
            geom_text(
                aes_string(label=label),
                color = 'black',
                size = 2,
                na.rm=TRUE
            )
        }

    return(fig)
}


#' Plot Amp Curves
#' 
plot_amp_curves <- function(
    amp_data,
    ct_threshold=NA,
    color='row_id',
    plate_id='1'
) {
    
    fig <- ggplot(amp_data[(amp_data['delta_rn'] > 0), ],
        aes(x=.data[['cycle']], y=.data[['delta_rn']],
            group=row_id,
            colour=!!sym(color))) +
        geom_line(alpha=0.7, size=0.5) +
        geom_hline(yintercept=ct_threshold, colour='red') +
        labs(x='Cycle Number', y='Delta Rn', title=paste("Amplification Curve for Plate", plate_id)) +
        scale_y_continuous(trans='log10', labels = function(x) round(x, 5), limits = c(0.00001, 15))

    return(fig)
}


#' Plot Fold Change
#' 
plot_fold_change <- function(
    df,
    x='tissue',
    y='fold_change_dnase1l1_actin',
    color='sample_id',
    xlabel='Tissue',
    ylabel='Fold Change',
    title=NULL
) {

    fig <- ggplot(
        df,
        aes(x=reorder(.data[[x]], .data[[y]], decreasing=TRUE),
            y=.data[[y]]), na.rm=TRUE) +
        # see: https://stackoverflow.com/questions/32642856/how-do-i-use-geom-rect-with-discrete-axis-values
        geom_rect(
            aes(xmin = stage(
                    reorder(.data[[x]], .data[[y]], decreasing=TRUE),
                    after_scale = xmin-0.45),
                xmax = stage(
                    reorder(.data[[x]], .data[[y]], decreasing=TRUE),
                    after_scale = xmax+0.45),
                ymin = .data[['min_fold_change']],
                ymax = .data[['max_fold_change']]),
                fill="#7f7f7f",
                stat='identity', inherit.aes=TRUE) +
        geom_jitter(aes(colour=.data[[color]])) +
        labs(x=xlabel, y=ylabel, title=title) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_y_continuous(trans='log10', labels = function(x) round(x, 5))

    return(fig)
}
