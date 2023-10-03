# Warning: If the libraries are not available, this script attempts to install!
# Settings: ###################################################################

# Location:
# https://stackoverflow.com/a/36276269
directory <- getSrcDirectory(function(dummy) {dummy})

source(file.path(directory, "Figure_Functions.R"))

#directory <- '.'
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

chunking <- 5000

plotOutput <- 1

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
# Files Total (Shattering) Fragmentation, v = 0.2, M = 1E4
files_analyses <- grep(
    x = dir(file.path(directory, "Analyses", "Data20-10-10"),
            full.names = TRUE),
    pattern = "_10000_",
    value = TRUE, fixed = TRUE
)#"Analyses18-07-06-18-09-09"))
files_data <- dir(file.path(
    directory, "Data20-10-10", "MixedKernel_MM_Total_v02"
), full.names = TRUE)

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
            M = 10000,
            time = time * 10, # Sample Time to Event Time.
            Simulation = gsub("^.+Run[_]", "", file),
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
                rss_xmin, by = c('time', 'xmax', 'M', 'Simulation',
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
                    rss_xmin, by = c('time', 'xmax', 'M', 'Simulation',
                                     'error', 'kernel', 'file'),
                    suffix = c('exponent', 'xmin')
                ) %>% dplyr::arrange(time)
        }
    }

    result_storage_subset <- result_storage_subset_out

    # Detect appropriate file from files_data to load.
    fname_run <- strsplit(fname, split = "_")
    print(fname_run)
    if (fname_run[[1]][7] == "Run") {
        fname_run = fname_run[[1]][8]
    } else {
        stop("Error in file loaded. Does not conform to expected convention.")
    }

    file_data <- grep(pattern = paste0("10000_Run_", fname_run),
                      x = files_data,
                      value = TRUE, fixed = TRUE)

    file_data_matrix <- R.matlab::readMat(
        file.path(file_data)
    )$stitch # Name of the Matrix in MATLAB format.

    # Add relevant values to result_storage_subset.
    result_storage_subset$N <- rowSums(file_data_matrix)
    result_storage_subset$n1 <- file_data_matrix[, 1]

    # Save.
    res_phy <- rbind(
        res_phy,
        result_storage_subset
    )

    print(paste(fname, "stored."))
}
if (exists("result_storage_subset") && exists("file_data_matrix")) {
    rm("result_storage_subset")
    rm("file_data_matrix")
    rm("rss_exponent")
    rm("rss_xmin")
    rm("result_storage_subset_out")
} else if (exists("result_storage_subset")) {
    stop("result_storage_subset DNE: was data ever loaded?")
} else if (exists("file_data_matrix")) {
    stop("file_data_matrix DNE: was data ever loaded?")
}

# Run and Save: ###############################################################
# Create figures separately.
if (plotOutput) {
    tempplot <- plotting_summary_statistics_4part(
        data = res_phy, # With columns xmax, M, N, exponent
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
            filename = paste0("Figure_SteadyState_",
                              names(tempplot)[plotindex],
                              '.pdf'),
            plot = tempplot[[plotindex]],
            path = folder_imag,
            width = 5.5, height = 4,
            dpi = 'retina'
        )
    }
}



print("-------------------------------------------------------------------")
for (sim in unique(res_phy$Simulation)) {
    temp <- (res_phy %>% dplyr::filter(Simulation == sim))$`MLE:exponent`
    print(paste0(sim, ', MLE: Exponent'))
    print(tseries::adf.test(temp))
    print(tseries::kpss.test(temp))

    temp <- (res_phy %>% dplyr::filter(Simulation == sim))$`KS-MLE:exponent`
    print(paste0(sim, ', KS-MLE: Exponent'))
    print(tseries::adf.test(temp))
    print(tseries::kpss.test(temp))

    temp <- diff((res_phy %>% dplyr::filter(Simulation == sim))$`MLE:exponent`)
    print(paste0(sim, ', diff(MLE: Exponent)'))
    print(tseries::adf.test(temp))
    print(tseries::kpss.test(temp))

    temp <- diff((res_phy %>% dplyr::filter(Simulation == sim))$`KS-MLE:exponent`)
    print(paste0(sim, ', diff(KS-MLE: Exponent)'))
    print(tseries::adf.test(temp))
    print(tseries::kpss.test(temp))
    print("-------------------------------------------------------------------")
}

# A natural first step is examine the distributional properties of our fitted summary statistics.
# To begin, we make a note about the presence of two regions in the KS-MLE plot in contrast to the MLE plot.
ggplot2::ggplot(
    res_phy %>% filter(`KS-MLE:exponent`<5, `KS-MLE:exponent`>1),
    ggplot2::aes(x = `KS-MLE:xmin`,
                 y = `KS-MLE:exponent`,
                 color = `KS-MLE:xmin`)
) + ggplot2::geom_bin2d()
# This is due to the sensitivity of the KS-MLE to the calculated $x_{\min}$ for which the power-law distribution begins.
# In practice, the higher values arise when $x_{\min} = 1$ is preferred by the algorithm, while lower values arise when $x_{\min} \geq 2$.
ggplot2::ggplot(
    res_phy %>% filter(`KS-MLE:exponent`<5, `KS-MLE:exponent`>1),
    ggplot2::aes(x = `MLE:exponent`,
                 y = `KS-MLE:exponent`,
                 color = `KS-MLE:xmin`,
                 group = `KS-MLE:xmin`)
) + ggplot2::geom_bin2d()
# As one would expect, all cases where KS-MLE has $x_{\min} = 1$ has the MLE equal to the KS-MLE.
with(res_phy %>% filter(`KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` <= 1),
     {any(`KS-MLE:exponent` != `MLE:exponent`)})
# On the other hand, there is not an obvious relationship between the MLE exponent and the KS-MLE exponent when $x_{\min}$ is otherwise, although there is still an expected correlation ($0.454$).
with(res_phy %>% filter(`KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` > 1),
     {cor(`KS-MLE:exponent`,`MLE:exponent`)})
ggplot2::ggplot(
    res_phy %>% filter(`KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` >= 2),
    ggplot2::aes(x = `MLE:exponent`,
                 y = `KS-MLE:exponent`)
) + ggplot2::geom_bin2d()
# The KS-MLE and MLE are then equal $61.6\%$ of the time.
sum(res_phy$`KS-MLE:xmin` == 1)/nrow(res_phy)
# The distributions of the KS-MLE and MLE exponents can be seen in the violin plots of Figure \ref{fig}
ggplot(
    rbind(
        res_phy %>% filter(
            `KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` == 1
        ) %>% rename(
            `KS-MLE:exponent, \nxmin = 1` = `KS-MLE:exponent`
        ) %>% gather(
            key = "Method", value = "Exponent",
            `MLE:exponent`,
            `KS-MLE:exponent, \nxmin = 1`,
        ),
        res_phy %>% filter(
            `KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` > 1
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
with(
    res_phy,
    c(quantile(`MLE:exponent`,
               p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)),
      "Mean" = mean(`MLE:exponent`))
)
with(
    res_phy %>% filter(
        `KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` == 1
    ),
    c(quantile(`KS-MLE:exponent`,
               p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)),
      "Mean" = mean(`KS-MLE:exponent`))
)
with(
    res_phy %>% filter(
        `KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` > 1
    ),
    c(quantile(`KS-MLE:exponent`,
               p = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)),
      "Mean" = mean(`KS-MLE:exponent`))
)
# %Divergence from normality appears to occur in the tails only
with(
    res_phy %>% filter(
        `KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` > 1
    ),
    {qqnorm(`KS-MLE:exponent`); qqline(`KS-MLE:exponent`)}
)
with(
    res_phy %>% filter(
        `KS-MLE:exponent`<5, `KS-MLE:exponent`>1, `KS-MLE:xmin` == 1
    ),
    {qqnorm(`KS-MLE:exponent`); qqline(`KS-MLE:exponent`)}
)
with(
    res_phy %>% filter(
        `KS-MLE:exponent`<5, `KS-MLE:exponent`>1
    ),
    {qqnorm(`MLE:exponent`); qqline(`MLE:exponent`)}
)


sum(res_phy$`KS-MLE:exponent` <= 2.5) / nrow(res_phy)
sum((res_phy %>% filter(`KS-MLE:exponent` <= 2.5))$`KS-MLE:exponent` < 1) /
    nrow((res_phy %>% filter(`KS-MLE:exponent` <= 2.5)))
