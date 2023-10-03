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
      c(file.path(directory, "Analyses", "Data19-07-08"),
        file.path(directory, "Analyses", "Data20-06-03"),
        file.path(directory, "Analyses", "Data23-04-18")
      ),
      full.names = TRUE, recursive = TRUE
    ),
    pattern = "_MM", # Population Size.
    value = TRUE, fixed = TRUE
  )
)

# Mixed loading scheme due to mixed filename types.
res_phy <- NULL
for (f in files_analyses) {
  print(f)
  fname <- basename(f)
  load(f) # loads "result_storage_subset"

  fname_split <- strsplit(fname, "_")[[1]]
  if(fname_split[3] == "MM") {
    Pop <- as.numeric(fname_split[6])
  } else {
    Pop <- as.numeric(fname_split[4])
  }

  # Filter to data we plan to plot.
  result_storage_subset <-
    result_storage_subset %>% dplyr::mutate(
      method_type = Method_conversion[as.character(method_type)]
      # ) %>% dplyr::filter(
      #     method_type == Method_loaded
    ) %>% dplyr::mutate(
      file = fname,
      M = Pop,
      time = time * 10, # Sample Time to Event Time.
      Simulation = gsub("^.+Run[_]", "", file),
      kernel = substring(file, 29, 30)
      #     n1 = n_1
      # ) %>% dplyr::select(
      #     -n_1
    ) %>% dplyr::rename(
      N = clusters
    )

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
        rss_xmin, by = c('time', 'xmax', 'M', #'n1',
                         'N', 'Simulation',
                         'error', 'kernel', 'file'),
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
        rss_xmin, by = c('time', 'xmax', 'M', #'n1',
                         'N', 'Simulation',
                         'error', 'kernel', 'file'),
        suffix = c('exponent', 'xmin')
      ) %>% dplyr::arrange(time)
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

# Retrieve Fragmentation Data from File Name.
res_phy$Fragmentation <- as.numeric(paste0("0.", regmatches(
  x = res_phy$file,
  m = regexpr( text = res_phy$file, pattern = "(?<=v0)[0-9]+", perl = T)
)))

# r is the variable that predicts the cyclicity.
res_phy$r <- with(res_phy, Fragmentation * M / (1 - Fragmentation))

# Retrieve the cyclicity:
load(file.path(directory, "summarystats2-2020-07-10.RData"))
res_phy <- res_phy %>% dplyr::left_join(
  summarystats2 %>% dplyr::select(Pop, PFr, Avg_Signum_Xmax_over_simulation),
  by = c("M" = "Pop", "Fragmentation" = "PFr")
)

### Note: Why not calculate the cyclicity? ############
# In principle, this would be fine to do! In practice,#
# our simulations enforced fragmentation for speed    #
# when there was nothing else that the simulation     #
# could do. To control for that (especially when the  #
# fragmentation probability was low), we needed to    #
# estimate how many steps we had skipped (which was   #
# not recorded in the original simulation). Hence,    #
# we load the results of that calculation instead.    #
#######################################################

# Run and Save: ###############################################################
# Create figures separately.
groups <- unique(interaction(res_phy$kernel, res_phy$M, res_phy$Fragmentation))
if (plotOutput) {
  for (group in groups) {
    print(paste("Graph: ", group))

    target <- res_phy %>% dplyr::filter(
      interaction(kernel, M, Fragmentation) == group
    )

    target_details <- target %>% dplyr::select(
      kernel, M, Fragmentation, r, Avg_Signum_Xmax_over_simulation
      ) %>% dplyr::distinct()

    tempplot <- plotting_summary_statistics_4part(
      data = target, # With columns xmax, M, N, exponent
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
        parameter_r <- res_phy %>% dplyr::filter(
          interaction(kernel, M, Fragmentation) == group
        ) %>% dplyr::pull(
          r
        ) %>% unique %>% formatC(digits = 4)

        parameter_cycl <- res_phy %>% dplyr::filter(
          interaction(kernel, M, Fragmentation) == group
        ) %>% dplyr::pull(
          Avg_Signum_Xmax_over_simulation
        ) %>% unique %>% formatC(digits = 4)

        tags <- paste0("(", parameter_r, ", ", parameter_cycl, ")")
        tempplot[[plotindex]] <- tempplot[[plotindex]] + ggplot2::labs(
          tag = tags#str2lang(tags)
        )
        tempplot[[plotindex]]$theme$plot.tag.position <- c(0.20, 0.06)
        tempplot[[plotindex]]$theme$plot.tag$size <- rel(0.75)
      }

      ggplot2::ggsave(
        filename = with(
          target_details,
          paste0("Figure_SteadyState_rK_",
                 kernel, '_', M, "_", Fragmentation, "_",
                 names(tempplot)[plotindex],
                 '.png')
          ),
        plot = tempplot[[plotindex]],
        path = folder_imag,
        width = 5.5, height = 4,
        dpi = 'retina'
      )
    }

  }
}
