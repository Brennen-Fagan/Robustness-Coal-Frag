# Load and Average (pdf + ccdf)
# Borrowing from CoFr_PaperImages_20_07_09_rnotebook.Rmd
imagedpi <- 350#"retina"
imageunits <- "mm"
imageheight <- 297
imagewidth <- 420
imageheight_small <- 120
imagewidth_small <- 148
imagepath <- file.path(directory, "Images2023-08-10")

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(VGAM)
library(R.matlab)

directory <- file.path('.')
data_folders <- dir(directory, pattern = "Data", full.names = TRUE)
analyses_folders <- dir(file.path(directory, "Analyses"),
                        pattern = "Data", full.names = TRUE)

data_folders_targets <- dir(data_folders[7])

for (data_folders_target in data_folders_targets) {
  # data_folders_target <- "SteadyState_ImpAll_v020"

  # Functions: #################################################################
  plotting_time_average <- function(data, palette = c("black", "red")) {
    data <- data %>% dplyr::mutate(
      # Empirical = grepl(x = Solution, pattern = "Empirical"),
      Fitted = grepl(x = Solution, pattern = "Fitted", fixed = TRUE),
      Complementary = grepl(x = Solution,
                            pattern = "Complementary", fixed = TRUE),
      Complementary = factor(ifelse(Complementary,
                                    "At Least Group Size",
                                    "Equal to Group Size"))
    )

    plot_text <- tibble::tibble(
      Group_Size = c(5E4, 5E4, 1E2, 1E2),#20,
      Value = c(6E3, 7E3, 3E4, 1E4),
      #   min(data %>% dplyr::filter(
      #   Empirical == TRUE, Value > 0
      # ) %>% dplyr::pull(Value)),
      Complementary = factor(ifelse(c(TRUE, FALSE, TRUE, FALSE),
                                    "At Least Group Size",
                                    "Equal to Group Size")),
      text = c(
        '"a) " * sum(n[i], i >= k, M)',
        '"b) n"[k]',
        # "#(Group \u2265 Size)", "#(Groups of Size)"
        "alpha: 2.5 * ', ' * Slope: 1.5",
        "alpha: 2.5 * ', ' * Slope: 2.5"
      )
    )

    figure <- ggplot2::ggplot(
      data %>% dplyr::filter(!Fitted),
      ggplot2::aes(
        x = Group_Size,
        y = Value,
        color = Fitted
      ),
    ) + ggplot2::geom_point(
      size = 4
    ) + ggplot2::geom_line(
      data = data %>% dplyr::filter(Fitted),
      size = 3
    ) + ggplot2::geom_abline(
      data = data.frame(
        slope = -2.5 + c(1, 0),
        intercept = c(5, 5),
        Complementary = factor(ifelse(c(TRUE, FALSE),
                                      "At Least Group Size",
                                      "Equal to Group Size"))
      ),
      mapping = ggplot2::aes(
        slope = slope,
        intercept = intercept
      ),
      size = 2, linetype = "dashed", color = "blue"
    ) + ggplot2::geom_text(
      data = plot_text,
      mapping = ggplot2::aes(
        x = Group_Size,
        y = Value,
        label = text
      ),
      size = 6,
      inherit.aes = FALSE,
      parse = TRUE
    ) + ggplot2::facet_wrap(
      Complementary ~ ., nrow = 2, ncol = 1,
      scales = "free_y"#, switch = "y"
    ) + ggplot2::labs(
      #x = "Group Size",
      x = "Size, k",
      y = "Time Averaged Number"
      #y = "Time Averaged Number of Groups"
    ) + ggplot2::scale_x_log10(
    ) + ggplot2::scale_y_log10(
    )

    if (!is.null(palette)) {
      figure <- figure + ggplot2::scale_colour_manual(
        values = palette
      )
    }

    if (!figures_legends) {
      figure <- figure + ggplotTheme(
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank()
      )
    } else {
      figure <- figure + ggplotTheme(
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank()
      )
    }

    return(figure)
  }

  # Values for pmf evaluated at N.
  solution_powerlaw <- function(N, exponent, min = 1, max = NULL,
                                complementary = FALSE) {
    if (is.null(max)) {
      if (!complementary) {
        retVal <- (N ^ -(exponent)) / VGAM::zeta(x = exponent, shift = min)
      } else {
        retVal <- (VGAM::zeta(x = exponent, shift = N)
                   / VGAM::zeta(x = exponent, shift = min))
      }
    }  else {
      vals <- (min:max)^(-exponent)
      retVal <- vals[N - min + 1]/sum(vals)
      retVal[is.na(retVal)] <- 0
      if (complementary) {
        retVal <- rev(cumsum(rev(retVal)))
      }
    }
    return(retVal)
  }

  # Values for solution unnormalised for 1:N.
  solution_analytic <- function(N, v, barrier = NULL,
                                complementary = FALSE) {
    if (is.null(barrier)) {
      retVal <- ((1/((2 - v) * sqrt(pi))
                  * exp(lgamma(1:N - 1 / 2) - lgamma(1:N))
                  * (4*(1 - v)/(2 - v) ^ 2) ^ (1:N - 1)
                  * (1:N) ^ (-2))) * N
      if (complementary) {
        retVal <- rev(cumsum(rev(retVal)))
      }
    } else {
      retVal <- solution_powerlaw(1:N, exponent = 2.5, min = 1, max = barrier,
                                  complementary = complementary) * N
    }
    return(retVal)
  }

  # Load: ######################################################################
  # Load only sparse files? If false, loads only non-sparse files.
  sparse_mode <- 1
  file_mat_pattern <- if (sparse_mode) '.*_S[.]mat$' else '.*[^(_S)][.]mat$'

  ### Summary Statistics: ######################################################
  files_summarystat <- dir(file.path(
    analyses_folders, data_folders_target
  ), full.names = TRUE)
  res_phy <- NULL
  for (file in files_summarystat) {
    load(file)
    if(is.null(res_phy)) {
      res_phy <- result_storage_subset
    } else {
      res_phy <- rbind(res_phy, result_storage_subset)
    }
  }

  ### Matrices: ################################################################
  # Load matrices to perform column summations.
  files_1 <- dir(file.path(data_folders, data_folders_target),
                 full.names = TRUE, pattern = file_mat_pattern)

  figure_1_Pop <- unique(unlist(lapply(strsplit(unlist(lapply(
    files_1,
    basename)), "_"), function(x) x[[5]])))
  figure_1_Pop <- as.numeric(figure_1_Pop)
  stopifnot(!is.na(figure_1_Pop))

  figure_1_bar <- 10000

  figure_1_matrix <- matrix(data = 0,
                            nrow = 1,
                            ncol = figure_1_Pop)
  figure_1_rows <- 0

  for (fi in files_1) {
    d <-
      if (substr(fi, nchar(fi) - 3,
                 nchar(fi)) == ".mat") {
        R.matlab::readMat(fi)
      } else {
        NULL
      }
    if (is.null(d)) next

    for (mat in d) {
      figure_1_matrix[1, ] <- figure_1_matrix[1, ] + colSums(mat)
      figure_1_rows <- figure_1_rows+ nrow(mat)
    }
  }

  ### Summary Statistics: ######################################################
  figure_1_median_exponent <- median(res_phy %>% dplyr::filter(
    method_type == "Basic Truncated MLE"
  ) %>% pull(exponent))

  figure_1_median_target <- which.min(abs(res_phy %>% dplyr::filter(
    method_type == "Basic Truncated MLE"
  ) %>% pull(exponent) - figure_1_median_exponent))

  figure_1_median_expAndMax <- (res_phy %>% dplyr::filter(
    method_type == "Basic Truncated MLE"
  ))[figure_1_median_target, c(3, 5)]

  ### To Ggplot: ###############################################################
  # Convert the loaded matrix into a data frame for ggplotting.
  figure_1_data <- data.frame(
    Group_Size = 1:(figure_1_Pop),
    Empirical = (figure_1_matrix / figure_1_rows)[1, ],
    Empirical_Complementary = rev(cumsum(rev(
      figure_1_matrix / figure_1_rows
    ))),
    # Analytical = solution_analytic(figure_1_Pop_bar, v = figure_1_PFr,
    #                                barrier = 10^4),
    # Complementary = solution_analytic(figure_1_Pop_bar, v = figure_1_PFr,
    #                                   barrier = 10^4, complementary = TRUE),
    Fitted = solution_powerlaw(
      1:figure_1_Pop,
      exponent = figure_1_median_expAndMax[[2]],
      min = 1, max = figure_1_median_expAndMax[[1]],
      complementary = FALSE) * figure_1_Pop,
    Fitted_Complementary = solution_powerlaw(
      1:figure_1_Pop,
      exponent = figure_1_median_expAndMax[[2]],
      min = 1, max = figure_1_median_expAndMax[[1]],
      complementary = TRUE) * figure_1_Pop,
    stringsAsFactors = FALSE
  ) %>% tidyr::gather("Solution", "Value",
                      Empirical, Empirical_Complementary,
                      # Analytical, Complementary,
                      Fitted, Fitted_Complementary
  )

  ### Run and Save #############################################################
  figure_1 <- plotting_time_average(figure_1_data)

  # plot_number <- 1
  # if (plot_number %in% plots_to_print) {
  ggplot2::ggsave(
    plot = figure_1,
    filename = paste0("Figure", data_folders_target, "_MeanMedian.png"),
    path = imagepath,
    dpi = imagedpi,
    width = imagewidth_small,
    height = imageheight_small,
    units = imageunits,
    limitsize = FALSE
  )
  # }

}
