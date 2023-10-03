# Warning: If the libraries are not available, this script attempts to install!
# Settings: ###################################################################

# Location:
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

plotOutput <- T

chunking <- 500

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
# Files Total (Shattering) Fragmentation, v = 0.01, M = 1E4
files_analyses <- c(
  grep(
    x = dir(file.path(directory, "Analyses", "Data20-08-25"), full.names = TRUE),
    pattern = "_10000_",
    value = TRUE, fixed = TRUE
  )
)

files_data <- c(
  dir(file.path(
    directory, "Data20-08-25", "MixedKernel_MM_Total"
  ), full.names = TRUE),
  dir(file.path(
    directory, "Data20-08-25", "MixedKernel_MM_PL120"
  ), full.names = TRUE),
  dir(file.path(
    directory, "Data20-08-25", "MixedKernel_MM_PL150"
  ), full.names = TRUE),
  dir(file.path(
    directory, "Data20-08-25", "MixedKernel_MM_PL180"
  ), full.names = TRUE),
  dir(file.path(
    directory, "Data20-08-25", "MixedKernel_MM_PL190"
  ), full.names = TRUE),
  dir(file.path(
    directory, "Data20-08-25", "MixedKernel_MM_PL195"
  ), full.names = TRUE)
)

res_phy <- NULL
for (f in files_analyses) {
  print(f)
  fname <- basename(f)
  load(f) # loads "result_storage_subset"

  # Not all kernels appear to have been properly saved.
  fname_run <- strsplit(fname, split = "_")

  if (fname_run[[1]][7] == "Run") {
    #fname_run = fname_run[[1]][8]
    fname_join <- paste(fname_run[[1]][3:8], sep = "", collapse = "_")
    fkern <- paste(fname_run[[1]][2:5], sep = "", collapse = "_")
  } else if (fname_run[[1]][6] == "Run") {
    #fname_run = fname_run[[1]][7]
    fname_join <- paste(fname_run[[1]][3:7], sep = "", collapse = "_")
    fkern <- paste(fname_run[[1]][2:4], sep = "", collapse = "_")
  } else {
    stop("Error in file loaded. Does not conform to expected convention.")
  }


  # Filter to data we plan to plot.
  result_storage_subset <-
    result_storage_subset %>% dplyr::mutate(
      method_type = Method_conversion[as.character(method_type)]
      # ) %>% dplyr::filter(
      #     method_type == Method_loaded
    ) %>% dplyr::mutate(
      file = fname,
      M = 10000,
      time = time * 10, # Sample Time to Event Time.
      Simulation = gsub("^.+Run[_]", "", file),
      xmax = ifelse(is.infinite(xmax),
                    ifelse(time == lag(time), lag(xmax), lead(xmax)),
                    xmax),
      kernel = ifelse(is.na(kernel), fkern, kernel)
    )  %>% dplyr::select(
      -error
    ) %>% dplyr::arrange(
      time
    )

  n_1_Flag <- 'n_1' %in% colnames(result_storage_subset)

  if (n_1_Flag) {
    result_storage_subset <- result_storage_subset %>% dplyr::mutate(
      n1 = n_1
    ) %>% dplyr::select(
      -n_1
    )

    bywords <- c('time', 'xmax', 'M', 'n1', 'N', 'Simulation',
                 #'error',
                 'kernel', 'file')
  } else {

    bywords <- c('time', 'xmax', 'M', 'Simulation',
                 #'error',
                 'kernel', 'file')
  }

  print("Chunking")
  for (chunk in 1:(nrow(result_storage_subset)/chunking)) {
    rss_exponent <- result_storage_subset[
      ((chunk - 1) * chunking + 1) : (chunk * chunking),
    ] %>% dplyr::select(
      -xmin
    ) %>% tidyr::spread(
      key = method_type,
      value = exponent,
      drop = TRUE#,
      # sep = ':'
    )

    rss_xmin <- result_storage_subset[
      ((chunk - 1) * chunking + 1) : (chunk * chunking),
    ] %>% dplyr::select(
      -exponent
    ) %>% tidyr::spread(
      key = method_type,
      value = xmin,
      drop = TRUE#,
      # sep = ':'
    )

    if (chunk == 1) {
      result_storage_subset_out <- rss_exponent %>% dplyr::full_join(
        rss_xmin, by = bywords,
        suffix = c(':exponent', ':xmin')
      ) %>% dplyr::arrange(time)

      # Create storage ahead of time
      result_storage_subset_out <- result_storage_subset_out[
        rep(1:nrow(result_storage_subset_out),
            nrow(result_storage_subset)/chunking),
      ]
    } else {
      result_storage_subset_out[
        ((chunk - 1) * (chunking / 2) + 1) : (chunk * (chunking / 2)),
      ] <- rss_exponent %>% dplyr::full_join(
        rss_xmin, by = bywords,
        suffix = c('exponent', 'xmin')
      ) %>% dplyr::arrange(time)
    }
  }

  result_storage_subset <- result_storage_subset_out

  rows <- unique(result_storage_subset$time)/10

  print(paste("Length of rows:", length(rows)))
  if (!n_1_Flag) {
    # Detect appropriate file from files_data to load.
    print(fname_run)

    file_data <- grep(pattern = paste0("_", fname_join, "_"),
                      #paste0("_10000_Run_", fname_run),
                      x = files_data,
                      value = TRUE, fixed = TRUE)

    file_data_matrix <- R.matlab::readMat(
      file.path(file_data)
    )$stitch # Name of the Matrix in MATLAB format.

    # Add relevant values to result_storage_subset.
    result_storage_subset$N <- rowSums(file_data_matrix[rows, ])
    result_storage_subset$n1 <- file_data_matrix[rows, 1]
  }

  # Save.
  res_phy <- rbind(
    res_phy,
    result_storage_subset
  )

  print(paste(fname, "stored."))
}
if (exists("result_storage_subset")) { #&& exists("file_data_matrix")) {
  rm("result_storage_subset")
  #rm("file_data_matrix")
  rm("rss_exponent")
  rm("rss_xmin")
  rm("result_storage_subset_out")
} else if (!exists("result_storage_subset")) {
  stop("result_storage_subset DNE: was data ever loaded?")
  # } else if (!exists("file_data_matrix")) {
  #     stop("file_data_matrix DNE: was data ever loaded?")
}

# Run and Save: ###############################################################
# Create figures separately.
kerns <- unique(res_phy$kernel)
if (plotOutput) {
  for (kern in kerns) {
    print(paste("Graph: ", kern))

    tempplot <- plotting_summary_statistics_4part(
      data = res_phy %>% dplyr::filter(
        kernel == kern
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
      ggplot2::ggsave(
        filename = paste0("Figure_GelShatter_",
                          kern, '_',
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

# Additional Workings: ########################################################

kernel_conversion = c(
  "Coalescence", "Coal. and Frag.", "Fragmentation",
  "Accretion", "Accr. and Attr.", "Attrition",
  "1.20", "1.50", "1.80"
)
names(kernel_conversion) <- kerns[c(26, 27, 31, 11, 17, 23, 2, 3, 4)]

res_phy_temp <- res_phy %>% dplyr::filter(
  grepl("ImpCo", kernel, fixed = TRUE)
  | grepl("BB1", kernel, fixed = TRUE)
  | grepl("Accr3Unq1", kernel, fixed = TRUE)
  | grepl("Attr3Unq1", kernel, fixed = TRUE)
  | grepl("AtAc3Unq1", kernel, fixed = TRUE)
  | grepl("PL", kernel, fixed = TRUE)
) %>% dplyr::mutate(
  Panel = dplyr::case_when(
    grepl("PL", kernel, fixed = TRUE) ~ "c) Power-law fragmentation",
    grepl("Unq", kernel, fixed = TRUE) ~ "a) Accretion and attrition",
    TRUE ~ "b) Stick-breaking variants"
  ),
  Variant = factor(kernel_conversion[kernel], ordered = TRUE,
                   levels = kernel_conversion)
)

# variant_violins_cyclic, 600 wide, 800 tall
tempplot <- ggplot2::ggplot(
  res_phy_temp %>% dplyr::filter(
    # `KS-MLE:exponent` > 1, `KS-MLE:exponent` < 6
  ),
  ggplot2::aes(x = Variant, y = `KS-MLE:exponent`, group = Variant)
) + ggplot2::geom_violin(
  adjust = 0.125,
  draw_quantiles = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
) + ggplot2::facet_wrap(
  . ~ Panel, nrow = 3, scales = "free_x"
) + ggplotTheme() + ggplot2::coord_cartesian(ylim = c(0.5, 3.5))

# method_violins_gelshatter, 900 x 600
ggplot(
  rbind(
    res_phy %>% filter(
      kernel == "MixedKernel_MM_Total",
      `KS-MLE:exponent`<7, `KS-MLE:exponent`>0.5, `KS-MLE:xmin` == 1
    ) %>% rename(
      `KS-MLE:exponent, \nxmin = 1` = `KS-MLE:exponent`
    ) %>% gather(
      key = "Method", value = "Exponent",
      `MLE:exponent`,
      `KS-MLE:exponent, \nxmin = 1`,
    ),
    res_phy %>% filter(
      kernel == "MixedKernel_MM_Total",
      `KS-MLE:exponent`<7, `KS-MLE:exponent`>0.5, `KS-MLE:xmin` > 1
    ) %>% rename(
      `KS-MLE:exponent, \nxmin > 1` = `KS-MLE:exponent`
    ) %>% gather(
      key = "Method", value = "Exponent",
      `MLE:exponent`,
      `KS-MLE:exponent, \nxmin > 1`
    )
  ),
  ggplot2::aes(
    x = factor(Method, ordered = TRUE,
               levels = c('MLE:exponent',
                          'KS-MLE:exponent, \nxmin = 1',
                          'KS-MLE:exponent, \nxmin > 1')),
    y = Exponent
  )
) + geom_violin(
  adjust = 0.125,
  draw_quantiles = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
) + ggplotTheme(
  legend.position = "none"
) + ggplot2::scale_x_discrete(name = "Method")

res_phy %>% filter(
  kernel == "MixedKernel_MM_Total",
  `KS-MLE:exponent`<7, `KS-MLE:exponent`>0.5, `KS-MLE:xmin` == 1
) %>% dplyr::pull(`KS-MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy %>% filter(
  kernel == "MixedKernel_MM_Total",
  `KS-MLE:exponent`<7, `KS-MLE:exponent`>0.5, `KS-MLE:xmin` > 1
) %>% dplyr::pull(`KS-MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy %>% filter(
  kernel == "MixedKernel_MM_Total"
) %>% dplyr::pull(`MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)

res_phy %>% group_by(
  file
) %>% mutate(
  signumXmax = sign(xmax - lag(xmax))
) %>% ungroup(
) %>% group_by(
  kernel
) %>% summarise(
  meanK = mean(signumXmax, na.rm = TRUE)
) -> cyclicity

tempplotshort <- plotting_summary_statistics_4part(
  data = res_phy %>% dplyr::filter(
    kernel == kerns[27], time < 1E4
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
  reduce_xmax_scale_by = 100,
  # reduce_N_scale_by = 10000
  pointsize = 0.8, pointalpha = 0.8
)
ggplot2::ggsave(
  filename = paste0("Figure_GelShatterShort_",
                    kerns[27], '_',
                    names(tempplotshort)[1],
                    '.png'),
  plot = tempplotshort[[1]],
  path = folder_imag,
  width = 5.5, height = 4,
  dpi = 'retina'
)
ggplot2::ggsave(
  filename = paste0("Figure_GelShatterShort_",
                    kerns[27], '_',
                    names(tempplotshort)[2],
                    '.png'),
  plot = tempplotshort[[2]],
  path = folder_imag,
  width = 5.5, height = 4,
  dpi = 'retina'
)

for (kindex in c(28:33)) {
  tempplotshort <- plotting_summary_statistics_4part(
    data = res_phy %>% dplyr::filter(
      kernel == kerns[kindex], time < 1E4
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
    reduce_xmax_scale_by = 100,
    # reduce_N_scale_by = 10000
    pointsize = 0.8, pointalpha = 0.8
  )
  ggplot2::ggsave(
    filename = paste0("Figure_GelShatterShort_",
                      kerns[kindex], '_',
                      names(tempplotshort)[1],
                      '.png'),
    plot = tempplotshort[[1]],
    path = folder_imag,
    width = 5.5, height = 4,
    dpi = 'retina'
  )
  ggplot2::ggsave(
    filename = paste0("Figure_GelShatterShort_",
                      kerns[kindex], '_',
                      names(tempplotshort)[2],
                      '.png'),
    plot = tempplotshort[[2]],
    path = folder_imag,
    width = 5.5, height = 4,
    dpi = 'retina'
  )
}

res_phy %>% filter(
  kernel == "MixedKernel_MM_Total", `MLE:exponent` > 0, `MLE:exponent` < 100
) %>% group_by(kernel) %>% pull(`MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy_temp %>% filter(
  Variant == "Coalescence", `MLE:exponent` > 0, `MLE:exponent` < 100
) %>% group_by(Variant) %>% pull(`MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy_temp %>% filter(
  Variant == "Coal. and Frag.", `MLE:exponent` > 0, `MLE:exponent` < 100
) %>% group_by(Variant) %>% pull(`MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy_temp %>% filter(
  Variant == "1.80", `MLE:exponent` > 0, `MLE:exponent` < 100
) %>% group_by(Variant) %>% pull(`MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)

res_phy %>% filter(
  kernel == "MixedKernel_MM_Total",
  `KS-MLE:exponent` > 0, `KS-MLE:exponent` < 100, `KS-MLE:xmin` > 1
) %>% group_by(kernel) %>% pull(`KS-MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy_temp %>% filter(
  Variant == "Coalescence",
  `KS-MLE:exponent` > 0, `KS-MLE:exponent` < 100, `KS-MLE:xmin` > 1
) %>% group_by(Variant) %>% pull(`KS-MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy_temp %>% filter(
  Variant == "Coal. and Frag.",
  `KS-MLE:exponent` > 0, `KS-MLE:exponent` < 100, `KS-MLE:xmin` > 1
) %>% group_by(Variant) %>% pull(`KS-MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
res_phy_temp %>% filter(
  Variant == "1.80",
  `KS-MLE:exponent` > 0, `KS-MLE:exponent` < 100, `KS-MLE:xmin` > 1
) %>% group_by(Variant) %>% pull(`KS-MLE:exponent`) %>% quantile(
  p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
)
