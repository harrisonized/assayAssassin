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
## plot_heatmap
## plot_scatter
## plot_lines
## plot_dots_and_bars


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
            tab['label']
        } else {
            tab['label'] = lapply(
                tab['val'], 
                function(x) round(x, digits)
            )
        }
        # tab[['label']] <- as.character(tab[['label']])
    } else {
        label = NULL
    }

    fig <- ggplot(tab, aes(x=.data[[x]], y=.data[[y]], fill=.data[[fill]])) +
        geom_tile(color="white", lwd=0.3, linetype=1, na.rm=FALSE) +
        scale_y_discrete(limits=rev) +
        coord_fixed(expand=TRUE) +
        labs(title = title) +
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


#' Plot Scatter
#' 
plot_scatter <- function(
    df,
    x,
    y,
    color=NULL,  # NULL for no groups
    xlabel=NULL,
    ylabel=NULL,
    title=NULL,
    alpha=0.7,
    point_size=0.5,
    log_x=FALSE,
    log_y=FALSE,
    legend_large_circle=TRUE
) {

    # group.by
    if (is.null(color)) { colour <- NULL } else { colour <- sym(color) }

    # scale axes
    if (log_x) {
        scale_x <- scale_x_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    } else { scale_x <- list() }
    if (log_y) {
        scale_y <- scale_y_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    } else { scale_y <- list() }

    if (legend_large_circle) {
        guide <- guides( colour=guide_legend(override.aes=list(size=4L,alpha=1)) )
    } else {
        guide <- list()
    }

    # plot
    fig <- ggplot(df,
              aes(x=.data[[x]], y=.data[[y]],
                  colour=!!colour) ) +
       geom_point(alpha=alpha, size=point_size, na.rm=TRUE) +
       labs(x=xlabel, y=ylabel, title=title) + 
       guide +
       scale_x +
       scale_y

    return(fig)
}


#' Plot Lines
#' 
plot_lines <- function(
    df,
    x,
    y,
    group='gene',
    color='row_id',
    thresholds=NA,
    xlabel=NULL,
    ylabel=NULL,
    title=NULL,
    log_y=TRUE
) {
    if (log_y) {
        scale_y <- scale_y_continuous(
            trans='log10',
            labels = function(x) round(x, 5),
            limits = c(0.00001, 15)
        )
    } else { scale_y <- list() }

    if (!is.na(thresholds)) {
        hline <- geom_hline(yintercept=thresholds, colour='red')
    } else {
        hline <- list()
    }
    
    fig <- ggplot(df[(df[[y]] > 0), ],
        aes(x=.data[[x]], y=.data[[y]],
            group=.data[[group]],
            colour=.data[[color]])) +
        geom_line(alpha=0.7, size=0.5, na.rm=TRUE) +
        hline +
        labs(x=xlabel, y=ylabel, title=title) +
        scale_y

    return(fig)
}


#' Plot Dots and Bars
#'
#' @description
#' Like a violin plot, but with a bar that spans the dots instead
#' 
plot_dots_and_bars <- function(
    df,
    x='tissue',
    y='fold_change_dnase1l1_actin',
    color='sample_id',
    ymin='min_fold_change',
    ymax='max_fold_change',
    xlabel='Tissue',
    ylabel='Fold Change',
    title=NULL,
    log_y=TRUE
) {

    if (log_y) {
        scale_y <- scale_y_continuous(trans='log10', labels = function(x) round(x, 5))
    } else {
        scale_y <- list()
    }


    fig <- ggplot(
        df,
        aes(x=reorder(.data[[x]], .data[[y]], decreasing=TRUE),
            y=.data[[y]]), na.rm=TRUE) +
        # see: https://stackoverflow.com/questions/32642856/how-do-i-use-geom-rect-with-discrete-axis-values
        geom_rect(
            aes(xmin = stage(
                    .data[[x]],
                    after_scale = xmin-0.45),
                xmax = stage(
                    .data[[x]],
                    after_scale = xmax+0.45),
                ymin = .data[[ymin]],
                ymax = .data[[ymax]]),
                fill="#7f7f7f",
                stat='identity', inherit.aes=TRUE) +
        geom_jitter(aes(colour=.data[[color]])) +
        labs(x=xlabel, y=ylabel, title=title) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_y

    return(fig)
}
