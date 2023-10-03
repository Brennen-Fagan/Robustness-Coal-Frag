# Warning: If the libraries are not available, this script attempts to install!
# Settings: ###################################################################

# Location: Code-Model
# https://stackoverflow.com/a/36276269
directory <- getSrcDirectory(function(dummy) {dummy})

source(file.path(directory, "Figure_Functions.R"))

#directory <- '.'
# Location: Root.
directory <- file.path(directory, "..")

figures_legends <- FALSE

# Trim Proportion:
# Trim y scale data (xmax, exponent) if
#    y_dat < trim_multiplier[1] * quantile(y_dat, p = trim_thresholds[1]) ||
#    y_dat > trim_multiplier[2] * quantile(y_dat, p = trim_thresholds[2])
# Doing this prevents extreme data from being plotted and distorting the axes.
# Instead of deleting these data, we place them on the axes (set to +/-Inf).
trim_multiplier <- c(0.1, 1.1)
trim_thresholds <- c(0.1, 0.99)
# trim_multiplier <-  NULL
# trim_thresholds <- NULL

video_suffix <- Sys.Date()

Method_loaded <- (Methods <- c(
  "Basic Truncated MLE",
  "CSN + ZXX MLE and KS: Lower Bound"
))[1:2]

Method_conversion <- c(
  "Basic Truncated MLE" = 'MLE',
  "CSN + ZXX MLE and KS: Lower Bound" = 'KS-MLE'
)

plotOutput <- TRUE
labelOutput <- TRUE

chunking <- 5000

# Fix X11 dependency?
# https://stackoverflow.com/a/19916865
options(bitmapType = 'cairo')

# Adjust libraries: ###########################################################
print("Loading libraries")
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath)
}
.libPaths(c(librarypath, .libPaths()))

packages <- c(
  "ggplot2", "grid", "gridExtra", "scales", # Plotting Packages
  #"gganimate", "av",                        # Animation Packages
  "R.matlab",                               # Data loading
  "dplyr", "tidyr" #,                       # Data manipulation
)

for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, lib = librarypath,
                     repos = 'https://cloud.r-project.org',
                     dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

# Folders: ####################################################################
folder_imag <- file.path(directory, paste0("Images", video_suffix))
if (!dir.exists(folder_imag)) {
  dir.create(folder_imag, showWarnings = FALSE)
}

print(paste(
  folder_imag, if (dir.exists(folder_imag)) {"exists"} else {"DNE."}
))

# Loading: ####################################################################
files_analyses <- c(
  grep(
    x = dir(
      file.path(directory, "Analyses", "Data23-07-19"),
      full.names = TRUE, recursive = TRUE
    ),
    pattern = "_10000_", # Population Size.
    value = TRUE, fixed = TRUE
  )
)

res_phy <- NULL
for (f in files_analyses) {
  print(f)
  fname <- basename(f)
  load(f) # loads "result_storage_subset"

  # Filter to data we plan to plot.
  result_storage_subset <-
    result_storage_subset %>% dplyr::mutate(
      method_type = Method_conversion[as.character(method_type)]
      # ) %>% dplyr::filter(
      #     method_type == Method_loaded
    ) %>% dplyr::mutate(
      file = fname,
      # M = 10000,
      time = time * 10, # Sample Time to Event Time.
      Simulation = gsub("^.+Run[_]", "", file),
      #kernel =  substring(file, 29, 30)
      #     n1 = n_1
      # ) %>% dplyr::select(
      #     -n_1
    ) %>% dplyr::rename(
      N = clusters
    ) %>% dplyr::select(
      -kernel
    ) %>% tidyr::separate_wider_delim(
      file, delim = "_", names = c(
        NA, NA, "kernel", "PFrag", "M", NA, NA, NA, NA, NA, NA
      ), cols_remove = FALSE
    )

  for (chunk in 1:(nrow(result_storage_subset)/chunking)) {
    # Splitn1d2 has problems with xmax (one says inf, one says only group size)
    # and with the error, for similar reasons.
    rss_exponent <- result_storage_subset[
      ((chunk - 1) * chunking + 1) : (chunk * chunking),
    ] %>% dplyr::select(
      -xmin, -xmax, -error
    ) %>% tidyr::spread(
      key = method_type,
      value = exponent,
      drop = TRUE#,
      # sep = ':'
    )

    rss_xmin <- result_storage_subset[
      ((chunk - 1) * chunking + 1) : (chunk * chunking),
    ] %>% dplyr::select(
      -exponent, -xmax, -error
    ) %>% tidyr::spread(
      key = method_type,
      value = xmin,
      drop = TRUE#,
      # sep = ':'
    )

    rss_xmax <- result_storage_subset[
      ((chunk - 1) * chunking + 1) : (chunk * chunking),
    ] %>% dplyr::select(
      -exponent, -xmin, -error
    ) %>% tidyr::spread(
      key = method_type,
      value = xmax,
      drop = TRUE#,
      # sep = ':'
    )

    result_storage_subset_target <- rss_exponent %>% dplyr::full_join(
      rss_xmin, by = c('time', #'xmax',
                       'M', #'n1',
                       'N', 'Simulation', "PFrag",
                       #'error',
                       'kernel', 'file'),
      suffix = c(':exponent', ':xmin')
    ) %>% dplyr::full_join(
      rss_xmax, by = c('time', #'xmax',
                       'M', #'n1',
                       'N', 'Simulation', "PFrag",
                       #'error',
                       'kernel', 'file'),
      # suffix = c('', ':xmax')
    ) %>% dplyr::rename(
      `KS-MLE:xmax` = `KS-MLE`,
      `MLE:xmax` = `MLE`,
    ) %>% dplyr::arrange(time)

    if (chunk == 1) {
      result_storage_subset_out <- result_storage_subset_target

      # Create storage ahead of time
      result_storage_subset_out <- result_storage_subset_out[
        rep(1:nrow(result_storage_subset_out),
            nrow(result_storage_subset)/chunking),
      ]
    } else {
      result_storage_subset_out[
        ((chunk - 1) * (chunking / 2) + 1) : (chunk * (chunking / 2)),
      ] <- result_storage_subset_target
    }
  }

  result_storage_subset <- result_storage_subset_out

  # Save.
  res_phy <- rbind(
    res_phy,
    result_storage_subset
  )

  print(paste(fname, "stored."))
}
if (exists("result_storage_subset")) {
  rm("result_storage_subset")
  rm("rss_exponent")
  rm("rss_xmin")
  rm("result_storage_subset_out")
} else if (!exists("result_storage_subset")) {
  stop("result_storage_subset DNE: was data ever loaded?")
}

res_phy <- res_phy %>% dplyr::rowwise() %>% dplyr::mutate(
  M = as.numeric(M),
  xmax = min(`KS-MLE:xmax`, `MLE:xmax`)
)

# Run and Save: ###############################################################
# Create figures separately.
kerns <- unique(res_phy$kernel)
PFrags <- unique(res_phy$PFrag)
if (plotOutput) {
  for (kern in kerns) {
    print(paste("Graph: ", kern))
    for (pfrag in PFrags) {

      tempplot <- plotting_summary_statistics_4part(
        data = res_phy %>% dplyr::filter(
          kernel == kern, PFrag == pfrag
        ), # With columns xmax, M, N, exponent
        palette = c("#000000", # "black",
                    "#FF0000", # "red"
                    "#0000FF", # "blue"
                    "#AA33AA", # "purple"
                    "#AAAA00" # "yellow"
        ),
        trim_thresholds = trim_thresholds,  # Requires two values
        trim_multipliers = trim_multiplier, # Requires two values
        #times_to_animate = NULL, # Animate everything
        reduce_time_scale_by = 10000,
        reduce_xmax_scale_by = 100 #,
        # reduce_N_scale_by = 10000
      )

      for (plotindex in 1:length(tempplot)) {
        if(labelOutput) {
          # tags <- paste0(c("K[", "F["),
          #                strsplit(kern, split = "")[[1]],
          #                "]", collapse = "~")
          tags <- dplyr::case_match(kern,
                                    "BBa2b4"       ~ "Beta-Bin.(2, 4)",
                                    "BBa3b3"       ~ "Beta-Bin.(3, 3)",
                                    "BBa4b2"       ~ "Beta-Bin.(4, 2)",
                                    "Fragn1d2"     ~ "Shatter 1/2",
                                    "Fragn1d10"    ~ "Shatter 1/10",
                                    "Fragn9d10"    ~ "Shatter 9/10",
                                    "Splitn1d2"    ~ "Halve",
                                    "ImpAll"       ~ "Stick-break Both",
                                    "ImpComb"      ~ "Stick-break Coal.",
                                    "Partitioning" ~ "Partition UAR",
                                    "PL110"        ~ "CRP(1.10)",
                                    "PL120"        ~ "CRP(1.20)",
                                    "PL130"        ~ "CRP(1.30)",
                                    "PL140"        ~ "CRP(1.40)",
                                    "PL160"        ~ "CRP(1.60)",
                                    "PL170"        ~ "CRP(1.70)",
                                    "PL180"        ~ "CRP(1.80)",
                                    "PL190"        ~ "CRP(1.90)",
          )
          tempplot[[plotindex]] <- tempplot[[plotindex]] + ggplot2::labs(
            tag = tags#str2lang(tags)
          )
          tempplot[[plotindex]]$theme$plot.tag.position <- c(0.17, 0.06)
        }

        ggplot2::ggsave(
          filename = paste0("Figure_SteadyState_",
                            kern, '_', pfrag, '_',
                            names(tempplot)[plotindex],
                            '.png'),
          plot = tempplot[[plotindex]],
          path = folder_imag,
          width = 5.5, height = 4,
          dpi = 'retina'
        )
      }
    }

  }
}

# res_phy %>% dplyr::group_by(
#   kernel, PFrag, M
# ) %>% dplyr::summarise(
#   N_mean = mean(N),
#   N_median = median(N),
#   xmax_mean = mean(xmax),
#   xmax_median = median(xmax),
#   ksmle_025 = quantile(`KS-MLE:exponent`[`KS-MLE:xmin` > 1 & `KS-MLE:exponent` != 0], p = 0.025),
#   ksmle_500 = quantile(`KS-MLE:exponent`[`KS-MLE:xmin` > 1 & `KS-MLE:exponent` != 0], p = 0.5),
#   ksmle_975 = quantile(`KS-MLE:exponent`[`KS-MLE:xmin` > 1 & `KS-MLE:exponent` != 0], p = 0.975),
#   ksmle_0s = sum(`KS-MLE:xmin` > 1 & `KS-MLE:exponent` == 0),
#   ksmle_Count = sum(`KS-MLE:xmin` > 1),
#   ksmle_0Percentage = ksmle_0s / ksmle_Count * 100,
#   mle_025 = quantile(`MLE:exponent`, p = 0.025),
#   mle_500 = quantile(`MLE:exponent`, p = 0.5),
#   mle_975 = quantile(`MLE:exponent`, p = 0.975),
#   mle_Count = dplyr::n()
# ) %>% dplyr::mutate(
#   dplyr::across(N_mean:mle_Count,
#                 function(x) formatC(x, width = 5))
#   ) %>% write.table(file = "temp.tsv", sep = "\t", row.names = FALSE)
