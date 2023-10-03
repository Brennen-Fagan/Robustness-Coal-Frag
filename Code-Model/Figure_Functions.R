

# Functions: ##################################################################
ggplotTheme <- function(
    ..., base_size = 22
) {
    ggplot2::theme_bw(base_size) + ggplot2::theme(...)
}

# https://github.com/tidyverse/ggplot2/wiki/
# Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(
    ..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")
) {

    plots <- list(...)
    position <- match.arg(position)
    g <- ggplot2::ggplotGrob(
        plots[[1]] + ggplot2::theme(legend.position = position)
    )$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + ggplot2::theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(
        position,
        "bottom" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 1,
            heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)),
        "right" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 2,
            widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth))
    )

    # grid::grid.newpage()
    # grid::grid.draw(combined)

    # return gtable invisibly
    invisible(combined)
}

plotting_summary_statistics_4part <- function(
    data, # With columns xmax, M, N, exponent
    palette = c("#000000", # "black",
                "#FF0000", # "red"
                "#0000FF", # "blue"
                "#AA33AA", # "purple"
                "#AAAA00" # "yellow"
    ),
    trim_thresholds = NULL,  # Requires two values
    trim_multipliers = NULL, # Requires two values
    #times_to_animate = NULL, # Requires two values
    reduce_time_scale_by = 10000,
    reduce_xmax_scale_by = 10000 ,
    # reduce_N_scale_by = 10000
    pointsize = 0.5, pointalpha = 0.2
) {
    # Perform trimming for better axes. #########################################
    if (!is.null(trim_thresholds)) {
        if (is.null(trim_multiplier)) {
            xmax_exceed <- quantile(data$xmax, p = trim_thresholds,
                                    na.rm = TRUE) * c(0.9, 1.1)
            expt_exceed <- quantile(data$`KS-MLE:exponent`, p = trim_thresholds,
                                    na.rm = TRUE) * c(0.9, 1.1)
        } else {
            xmax_exceed <- quantile(data$xmax, p = trim_thresholds,
                                    na.rm = TRUE) * trim_multiplier
            expt_exceed <- quantile(data$`KS-MLE:exponent`, p = trim_thresholds,
                                    na.rm = TRUE) * trim_multiplier
        }

        data <- data %>% dplyr::mutate(
            xmax = ifelse(xmax > xmax_exceed[2], Inf, xmax),
            `KS-MLE:exponent` = ifelse(`KS-MLE:exponent` > expt_exceed[2], Inf, `KS-MLE:exponent`),
            xmax = ifelse(xmax < xmax_exceed[1], -Inf, xmax),
            `KS-MLE:exponent` = ifelse(`KS-MLE:exponent` < expt_exceed[1], -Inf, `KS-MLE:exponent`)
        )
    }

    # Identify area to color. ###################################################
    # if (!is.null(times_to_animate)) {
    #     data$highlighted <- ifelse(
    #         times_to_animate[1] <= data$time &
    #             data$time <= times_to_animate[2],
    #         data$Simulation,
    #         FALSE
    #     )
    # } else {
        data$highlighted <- data$Simulation
    # }
    data$highlighted <- factor(data$highlighted,
                               levels = c(FALSE, unique(data$Simulation)))

    # Prep heatmap bounding curves. #############################################
    data$xmax_norm <- data$xmax / data$M
    data$N_norm <- data$N / data$M

    curves <- data %>% dplyr::group_by(
        xmax_norm
    ) %>% dplyr::summarise(
        N_norm_hi = max(N_norm) + 0.005,
        N_norm_lo = min(N_norm) - 0.005
    ) %>% dplyr::ungroup(
    ) %>% tidyr::gather(
        key = "which", value = "Boundary", N_norm_hi, N_norm_lo
    ) %>% dplyr::filter(
        xmax_norm <= 1, xmax_norm >= 0,
        Boundary <= 1, Boundary >= 0
    )

    # Prep variable to animate over. ############################################
    # if (!is.null(times_to_animate)) {
    #     dynamicdata <- data %>% dplyr::filter(
    #         times_to_animate[1] <= time, time <= times_to_animate[2]
    #     )
    # } else {
        # dynamicdata <- data
    # }
    # dynamicdata$time_to_animate <- dynamicdata$time

    # Prep labels. ##############################################################
    reduce_time_by_word <- strsplit(gsub(formatC(reduce_time_scale_by, digits = 1),
                                         pattern = "e+0", replacement = "0^",
                                         fixed = TRUE), '^', fixed = TRUE)[[1]]
    reduce_xmax_by_word <- strsplit(gsub(formatC(reduce_xmax_scale_by, digits = 1),
                                         pattern = "e+0", replacement = "0^",
                                         fixed = TRUE), '^', fixed = TRUE)[[1]]
    # reduce_N_by_word <- strsplit(gsub(formatC(reduce_N_scale_by, digits = 1),
    #                                   pattern = "e+0", replacement = "0^",
    #                                   fixed = TRUE), '^', fixed = TRUE)[[1]]
    # https://github.com/thomasp85/gganimate/issues/332
    xmax_label <- as.expression(bquote(
        'k'[max] * ', ' %*% .(reduce_xmax_by_word[1])^.(reduce_xmax_by_word[2])
    ))
    time_label <- as.expression(bquote(
        't, ' %*% .(reduce_time_by_word[1])^.(reduce_time_by_word[2])
    ))

    # # OVERRIDE:
    # xmax_label <- "Largest Size, 1000s"
    # time_label <- "Time, 10000s"


    # Subfigure c+d) Exponent. ####################################################
    figure_expoksmle <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = time/reduce_time_scale_by,
            y = `KS-MLE:exponent`,#exponent,
            color = highlighted,
            shape = file
        )
    ) + ggplot2::geom_point(
        alpha = pointalpha,
        size = pointsize
    ) + ggplot2::labs(
        y = expression(KS-MLE~alpha), #"Power-law Exponent",#expression(MLE~alpha),
        x = ""
    # ) + ggplot2::annotate(
    #     "text",
    #     x = 12500/reduce_time_scale_by, y = 2.975, label = "c)", size = 6
    # ) + ggplot2::geom_point(
    #     data = dynamicdata,
    #     mapping = ggplot2::aes(
    #         x = time_to_animate/reduce_time_scale_by,
    #         y = exponent,
    #         color = highlighted,
    #         shape = file
    #     ),
    #     size = 4 # Bigger than the background points
    ) + ggplot2::scale_shape_manual(
        values = 15:18 # c(0:2,5)
    # ) + gganimate::transition_manual(
        # time_to_animate
    )

    figure_expomle <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = time/reduce_time_scale_by,
            y = `MLE:exponent`,#exponent,
            color = highlighted,
            shape = file
        )
    ) + ggplot2::geom_point(
        alpha = pointalpha,
        size = pointsize
    ) + ggplot2::labs(
        y = expression(MLE~alpha),#"Power-law Exponent",#expression(MLE~alpha),
        x = time_label
    # ) + ggplot2::annotate(
    #     "text",
    #     x = 12500/reduce_time_scale_by, y = 2.975, label = "d)", size = 6
        # ) + ggplot2::geom_point(
        #     data = dynamicdata,
        #     mapping = ggplot2::aes(
        #         x = time_to_animate/reduce_time_scale_by,
        #         y = exponent,
        #         color = highlighted,
        #         shape = file
        #     ),
        #     size = 4 # Bigger than the background points
    ) + ggplot2::scale_shape_manual(
        values = 15:18 # c(0:2,5)
        # ) + gganimate::transition_manual(
        # time_to_animate
    )

    # Subfigure b) xmax. ########################################################
    figure_xmax <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = time/reduce_time_scale_by,
            y = xmax/reduce_xmax_scale_by,
            color = highlighted,
            shape = file
        )
    ) + ggplot2::geom_point(
        alpha = pointalpha,
        size = pointsize
    ) + ggplot2::labs(
        # y = bquote(
        #   'k'[max] * ', ' %*% .(reduce_xmax_by_word[1])^.(reduce_xmax_by_word[2])
        # ),
        y = xmax_label,
        # x = bquote(
        #   't, ' %*% .(reduce_time_by_word[1])^.(reduce_time_by_word[2])
        # )
        x = time_label
    # ) + ggplot2::annotate(
    #     "text",
    #     x = 12500/reduce_time_scale_by,
    #     y = 1500/reduce_xmax_scale_by,#15500/reduce_xmax_scale_by,
    #     label = "b)", size = 6
    # ) + ggplot2::geom_point(
    #     data = dynamicdata,
    #     mapping = ggplot2::aes(
    #         x = time_to_animate/reduce_time_scale_by,
    #         y = xmax/reduce_xmax_scale_by,
    #         color = highlighted,
    #         shape = file
    #     ),
    #     size = 4 # Bigger than the background points
    ) + ggplot2::scale_shape_manual(
        values = 15:18 # c(0:2,5)
    # ) + gganimate::transition_manual(
        # time_to_animate
    )

    # Subfigure a) state space heatmap. #########################################
    figure_space <- ggplot2::ggplot(
        data, ggplot2::aes(x = xmax_norm, y = N_norm)
    # ) + ggplot2::geom_point(#ggplot2::geom_line(
    #     data = curves,
    #     mapping = ggplot2::aes(
    #         x = xmax_norm,
    #         y = Boundary,
    #         group = which
    #     ),
    #     color = "red",
    #     size = 0.5, # originally 1.3
    #     inherit.aes = FALSE
    ) + ggplot2::geom_bin2d( # Background Image Done.
        bins = 100
    # ) + ggplot2::geom_point( # Foreground animated dot:
    #     data = dynamicdata,
    #     mapping = ggplot2::aes(
    #         x = xmax_norm,
    #         y = N_norm,
    #         color = highlighted,
    #         shape = file
    #     ),
    #     size = 4,
    #     inherit.aes = FALSE
    ) + ggplot2::scale_shape_manual(
        values = 15:18 # c(0:2,5)
    ) + ggplot2::scale_fill_gradient2(
        low = "white", high = "black",
        trans = "log10", name = "", #"a) Count",
        breaks = scales::trans_breaks("log10", function(x) floor(10^x))
    ) + ggplot2::scale_x_continuous(
        name = expression('k'[max] / M),
        #name = "Normalised Largest Size"
        #breaks = c(0, 0.5, 1)
    ) + ggplot2::scale_y_continuous(
        name = expression(N/M),
        #name = "Normalised Number of Groups"
        #breaks = c(0, 0.5, 1)
    ) + ggplotTheme(
        # legend.position = c(0.87, 0.23),
        legend.position = c(0.88, 0.80),
        legend.background = ggplot2::element_rect(fill = "transparent", #"white",
                                                  color = "transparent"),
        #legend.spacing.y = ggplot2::unit(0, "mm"),
        # legend.margin = ggplot2::margin(-10, 0, -12, 0),
        # legend.text = ggplot2::element_text(size = 10),
        # legend.title = ggplot2::element_text(size = 12),
        # legend.key.height = ggplot2::unit(4, 'mm'),
        # panel.spacing.x = unit(4.5, "mm")#,
        # strip.background = ggplot2::element_blank(),
        # strip.text.x = ggplot2::element_blank()
    ) + ggplot2::guides(
        color = FALSE,
        shape = FALSE#,
        # fill = FALSE
        # ) + ggplot2::annotate(
        #   "text",
        #   x = .1, y = .1, label = "a)", size = 6
    # ) + gganimate::transition_manual(
    #     frames = time_to_animate
        )

    # Subfigure a) add curved arrows. ##########################################
    # figure_space <- figure_space + ggplot2::geom_path(
    #   data = data.frame(
    #     x = c(12, 6,
    #           seq(from = 0, to = 6, by = 1/3)) * 1000 / reduce_xmax_by,
    #     y = c(62.5, 67.5,
    #           3 * (((seq(from = 0, to = 6, by = 1/3)) - 6) / 6) ^ 4 + 57) * 1000 / reduce_N_by,
    #     group = c(1, 1,
    #               rep(2, length(seq(from = 0, to = 6, by = 1/3))))
    #   ),
    #   arrow = ggplot2::arrow(type = "closed", angle = 25,
    #                          length = unit(0.15, "inches")),
    #   mapping = ggplot2::aes(x = x, y = y, group = group),
    #   size = 1.3, color = if (!is.null(palette)) {
    #     palette[1]
    #   } else {
    #     "black"
    #   }
    # )
    #
    # Adding details. ###########################################################
    # if (!is.null(times_to_animate) ) {
    #     figure_expo <- figure_expo + ggplot2::geom_vline(
    #         xintercept = times_to_animate/reduce_time_by
    #     )
    #     figure_xmax <- figure_xmax + ggplot2::geom_vline(
    #         xintercept = times_to_animate/reduce_time_by
    #     )
    # }

    if (!is.null(palette)) {
        figure_expomle <- figure_expomle + ggplot2::scale_colour_manual(
            values = palette
        )
        figure_expoksmle <- figure_expoksmle + ggplot2::scale_colour_manual(
            values = palette
        )
        figure_xmax <- figure_xmax + ggplot2::scale_colour_manual(
            values = palette
        )

        figure_space <- figure_space + ggplot2::scale_color_manual(
            name = NULL, values = palette
        )
    }

    if (!figures_legends) {
        figure_expoksmle <- figure_expoksmle + ggplotTheme(legend.position = "none")
        figure_expomle <- figure_expomle + ggplotTheme(legend.position = "none")
        figure_xmax <- figure_xmax + ggplotTheme(legend.position = "none")
        #
        # figure <- gridExtra::arrangeGrob(
        #     figure_space,
        #     gridExtra::arrangeGrob(
        #       figure_expo,
        #       figure_xmax,
        #       ncol = 1, nrow = 2),
        #     ncol = 2, nrow = 1
        #   )
    } else {
        figure_expoksmle <- figure_expoksmle + ggplotTheme(
            legend.position = "bottom",
            legend.direction = "vertical"
        )
        figure_expomle <- figure_expomle + ggplotTheme(
            legend.position = "bottom",
            legend.direction = "vertical"
        )
        figure_xmax <- figure_xmax + ggplotTheme(
            legend.position = "bottom",
            legend.direction = "vertical"
        )
        #
        #   figure <- gridExtra::arrangeGrob(
        #       figure_space,
        #       grid_arrange_shared_legend(
        #         figure_expo,
        #         figure_xmax,
        #         ncol = 1, nrow = 2
        #       ),
        #       ncol = 2, nrow = 1
        #     )
    }
    return(list(
        mle = figure_expomle,
        ksmle = figure_expoksmle,
        xmax = figure_xmax,
        space = figure_space
    ))
}

