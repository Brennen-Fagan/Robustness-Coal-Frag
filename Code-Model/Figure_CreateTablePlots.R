# Warning: If the libraries are not available, this script attempts to install!
# Settings: ###################################################################

# Location: Code-Model
# https://stackoverflow.com/a/36276269
directory <- getSrcDirectory(function(dummy) {dummy})

source(file.path(directory, "Figure_Functions.R"))

#directory <- '.'
# Location: Root.
directory <- file.path(directory, "..")

Method_loaded <- (Methods <- c(
  "Basic Truncated MLE",
  "CSN + ZXX MLE and KS: Lower Bound"
))[1:2]

Method_conversion <- c(
  "Basic Truncated MLE" = 'MLE',
  "CSN + ZXX MLE and KS: Lower Bound" = 'KS-MLE'
)

# Memory vs Time trade-off.
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

# Loading: ####################################################################

### Analyses, so that we have the exponent regimes: ###########################
files_analyses <- c(
  grep(
    x = dir(
      file.path(directory, "Analyses", "Data20-08-25"),
      full.names = TRUE, recursive = TRUE
    ),
    pattern = "_10000_", # Population Size.
    value = TRUE, fixed = TRUE
  ),
  grep(
    x = dir(
      file.path(directory, "Analyses", "Data20-10-10"),
      full.names = TRUE, recursive = TRUE
    ),
    pattern = "_10000_", # Population Size.
    value = TRUE, fixed = TRUE
  ),
  grep(
    x = dir(
      file.path(directory, "Analyses", "Data23-07-19"),
      full.names = TRUE, recursive = TRUE
    ),
    pattern = "_10000_", # Population Size.
    value = TRUE, fixed = TRUE
  )
)

# Remover unused PFrag = 30%'s and MMBar, CM, MC, and CC kernels.
files_analyses <- files_analyses[!grepl(pattern = "v030", files_analyses)]
files_analyses <- files_analyses[!grepl(pattern = "_MMBar", files_analyses)]
files_analyses <- files_analyses[!grepl(pattern = "_CM_", files_analyses)]
files_analyses <- files_analyses[!grepl(pattern = "_MC_", files_analyses)]
files_analyses <- files_analyses[!grepl(pattern = "_CC_", files_analyses)]
files_analyses <- files_analyses[!grepl(pattern = "_Split", files_analyses)]

res_phy <- NULL
for (f in files_analyses) {
  #print(f)
  fname <- basename(f)
  load(f) # loads "result_storage_subset"

  retrievePFrag <-
    regexec(pattern = "(?<=_v)[0-9]+(?=_)",text = fname, perl = TRUE)[[1]]


  # Filter to data we plan to plot.
  result_storage_subset <-
    result_storage_subset %>% dplyr::mutate(
      method_type = Method_conversion[as.character(method_type)]
      # ) %>% dplyr::filter(
      #     method_type == Method_loaded
    ) %>% dplyr::mutate(
      file = fname,
      M = 10000,
      PFrag =
        if(retrievePFrag == -1) {
          0.01 # Following old standard from before PFrag exploration
        } else {
          as.numeric(paste0("0.", substr(
            fname, retrievePFrag,
            retrievePFrag + attr(retrievePFrag, "match.length") - 1
          ))) * 10
        },
      time = time * 10, # Sample Time to Event Time.
      Simulation = gsub("^.+Run[_]", "", file),
      Kernel = strsplit(fname, split = "_", fixed = TRUE)[[1]][3],
      Variant = strsplit(fname, split = "_", fixed = TRUE)[[1]][4],
      error = ifelse(
        is.infinite(xmax) & error == "NULL : ", FALSE,
        # A specific error that we account for versus the others we don't.
        error
      ),
      xmax = ifelse(xmax > M, M, xmax)
      # Rare inconsistency where xmax = Inf is reported.
    )

  if("clusters" %in% colnames(result_storage_subset)) {
    result_storage_subset <- result_storage_subset %>% dplyr::rename(
      N = clusters
    )
  }


  if(!"N" %in% colnames(result_storage_subset) ||
     !"n_1" %in% colnames(result_storage_subset) ) {
    # Now rare old standard from before PFrag exploration
    # Need to load the original files to extract.

    file_candidates <- dir(
      file.path(
        directory, if (grepl("Data20-08-25", f, fixed = TRUE)) {
          "Data20-08-25"
        } else if (grepl("Data20-10-10", f, fixed = TRUE)) {
          "Data20-10-10"
        } else if (grepl("Data23-07-19", f, fixed = TRUE)) {
          "Data23-07-19"
        }
      ),
      recursive = TRUE,
      pattern = paste(
        strsplit(fname, split = "_", fixed = TRUE)[[1]][1:8],
        collapse = "_"
      ),
      full.names = TRUE
    )

    if (length(file_candidates) == 0) {
      warning(paste(
        "File", f, "is missing one of n_1 or N and no substitutes found."
      ))
      next
    } else if (length(file_candidates) > 1) {
      warning("More than 1 candidate found, proceeding with 1st.")
      file_candidates <- file_candidates[1]
    }

    # Load the candidates, Arrange for Implicit Times, Identify N, and Join.

    populationMatrix <- R.matlab::readMat(file_candidates)[[1]] # Old Standard.
    populationNs <- rowSums(populationMatrix)
    result_storage_subset <- dplyr::left_join(
      result_storage_subset, data.frame(
        time = sort(unique(result_storage_subset$time)),
        N = populationNs,
        n_1 = populationMatrix[, 1]
      ),
      by = c("time",
             if("N" %in% colnames(result_storage_subset)) "N",
             if("n_1" %in% colnames(result_storage_subset)) "n_1"
      )
    )

    # Related, this older standard also had an error code 1 for MLE.
    result_storage_subset_errors <- with(result_storage_subset,
                                         error == 1 & N == 1 & xmin == 0)
    if (sum(result_storage_subset_errors) > 0) {
      result_storage_subset$error[result_storage_subset_errors] <- 0
      result_storage_subset$xmin[result_storage_subset_errors] <-
        result_storage_subset$xmax[result_storage_subset_errors]
    }
  }

  result_storage_subset <- result_storage_subset %>% dplyr::arrange(time)

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

    to_store <- rss_exponent %>% dplyr::full_join(
      rss_xmin, by = c('time', 'xmax', 'M', 'n_1', 'PFrag',
                       'N', 'Simulation',
                       'error', 'Kernel', "Variant", 'kernel', 'file'),
      suffix = c(':exponent', ':xmin')
    ) %>% dplyr::arrange(time)

    if (nrow(to_store) != chunking / 2) {
      warning(paste(
        "File", f, "has too many rows:",
        nrow(to_store), "received,",
        chunking / 2, "expected."
      ))
    }

    if (chunk == 1) {
      result_storage_subset_out <- to_store

      # Create storage ahead of time
      result_storage_subset_out <- result_storage_subset_out[
        rep(1:nrow(result_storage_subset_out),
            nrow(result_storage_subset)/chunking),
      ]
    } else {
      result_storage_subset_out[
        ((chunk - 1) * (chunking / 2) + 1) : (chunk * (chunking / 2)),
      ] <- to_store
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

### SummaryStats2, so that we have the cyclicities: ###########################

files_cyclicities <- dir(directory, pattern = "summarystats2-2023",
                         full.names = TRUE)

# The original 2020-07-10 version has a different format.
cyclicity <- load(file.path(directory, "summarystats2-2020-07-10.RData"))
cyclicity <- get(cyclicity) %>% dplyr::rename(
  M = Pop, PFrag = PFr
) %>% dplyr::mutate(
  kernel = "MixedKernel_MM_Total",
  source = "summarystats2-2020-07-10.RData"
)

for(f in files_cyclicities) {
  load(f)
  cyclicity <- dplyr::bind_rows(cyclicity,
                                summarystats2 %>% dplyr::mutate(
                                  source = basename(f)
                                ))
}

# Unify formats. Note this is done by inspection, so results may vary!
cyclicity <- cyclicity %>% dplyr::mutate(
  kernel = ifelse(
    nchar(kernel) < 13,
    paste("MixedKernel", kernel, "v020", sep = "_"),
    kernel
  )
)

powerlaws <- res_phy %>% dplyr::mutate(
  kernel = ifelse(
    is.na(kernel),
    paste("MixedKernel", Kernel, Variant, sep = "_"),
    kernel
  )
) %>% dplyr::group_by(
  kernel, PFrag, M
) %>% dplyr::summarise(
  N_mean = mean(N),
  N_median = median(N),
  xmax_mean = mean(xmax),
  xmax_median = median(xmax),
  ksmle_025 = quantile(`KS-MLE:exponent`[`KS-MLE:xmin` > 1 &
                                           `KS-MLE:exponent` != 0], p = 0.025),
  ksmle_500 = quantile(`KS-MLE:exponent`[`KS-MLE:xmin` > 1 &
                                           `KS-MLE:exponent` != 0], p = 0.5),
  ksmle_975 = quantile(`KS-MLE:exponent`[`KS-MLE:xmin` > 1 &
                                           `KS-MLE:exponent` != 0], p = 0.975),
  ksmle_0s = sum(`KS-MLE:xmin` > 1 & `KS-MLE:exponent` == 0),
  ksmle_Count = sum(`KS-MLE:xmin` > 1),
  ksmle_0Percentage = ksmle_0s / ksmle_Count * 100,
  mle_025 = quantile(`MLE:exponent`, p = 0.025),
  mle_500 = quantile(`MLE:exponent`, p = 0.5),
  mle_975 = quantile(`MLE:exponent`, p = 0.975),
  mle_Count = dplyr::n(),
  .groups = "drop"
)

# Filtering out unused simulations.
cyclicity <- cyclicity %>% dplyr::filter(
  M == 10000, PFrag %in% c(0.01, 0.20),
  !grepl(kernel, pattern = "[_](MC|CC|CM)[_]"),
  !grepl(kernel, pattern = "Bar")
  )

powerlaws <- powerlaws %>% dplyr::filter(
  M == 10000, PFrag %in% c(0.01, 0.20),
  !grepl(kernel, pattern = "[_](MC|CC|CM)[_]"),
  !grepl(kernel, pattern = "Bar")
  )

temp <- dplyr::full_join(
  cyclicity, powerlaws
) %>% dplyr::mutate(
  # Create a consistent shorthand.
  Shorthand = ifelse(grepl("_", kernel, fixed = TRUE),
                     unlist(lapply(strsplit(kernel, split = "_", fixed = TRUE),
                                   function(x) x[3])),
                     kernel),
  Shorthand = ifelse(Shorthand == "v020",
                     unlist(lapply(strsplit(kernel, split = "_", fixed = TRUE),
                                   function(x) x[2])),
                     Shorthand),
  Shorthand = ifelse(Shorthand == "Total",
                     "Base",
                     gsub("(Total|Frag)", "", x = Shorthand)
                     ),
  Shorthand = dplyr::case_when(
    Shorthand == "BB4AND2" | Shorthand == "BBa4b2" ~ "BB(4,2)",
    Shorthand == "BB3AND3" | Shorthand == "BBa3b3" ~ "BB(3,3)",
    Shorthand == "BB2AND4" | Shorthand == "BBa2b4" ~ "BB(2,4)",
    Shorthand == "ImpAll"| Shorthand == "ImpCoFr" ~ "ICF",
    Shorthand == "ImpCo" | Shorthand == "ImpComb" ~ "IC",
    Shorthand == "ImpFr" | Shorthand == "BB1AND1" ~ "IF",
    Shorthand == "n1d10" | Shorthand == "1OF10" ~ "1/10",
    Shorthand == "n1d2" | Shorthand == "1OF2" ~ "1/2",
    Shorthand == "n9d10" | Shorthand == "9OF10" ~ "9/10",
    Shorthand == "Splitn1d2" ~ "Half",
    Shorthand == "Partitioning" ~ "Part",
    TRUE ~ Shorthand
  ),
  Shorthand = gsub("(Accr|Ac)", "AC", Shorthand),
  Shorthand = gsub("(Attr|At)", "ER", Shorthand),
  Shorthand = gsub("Unq0", "", Shorthand),
  Shorthand = gsub("Unq1", "U", Shorthand),
  Family = dplyr::case_when(
    grepl("(AC|ER)", Shorthand) ~ "AC, ER",
    grepl("(/|Half)", Shorthand) ~ "Partial Frag.",
    grepl("(Part|I|BB)", Shorthand) ~ "Imperfect",
    grepl("PL", Shorthand) ~ "Power Law",
    grepl("Base", Shorthand) ~ "Base"
  )
)

# To prevent doubling, we'll prefer earlier files to later.
# (Neither should be more accurate, some were just used as double checks.)
temp <- temp %>% dplyr::group_by(
  Shorthand, PFrag
) %>% dplyr::arrange(
  source
) %>% dplyr::summarise(
  dplyr::across(dplyr::everything(), .fns = function(x, y) x[1]),
  .groups = "drop"
) %>% dplyr::mutate(
  PFrag = dplyr::case_when(
    PFrag == 0.01 ~ "Base: Gel-Shatter",
    PFrag == 0.20 ~ "Base: Steady-State",
  ),
  PFrag = factor(PFrag,
                 levels = c("Base: Steady-State", "Base: Gel-Shatter"),
                 ordered = TRUE)
)

removed_shorthands <- c(
  "ER1U", "ER1", "AC1U", "AC1", "ERAC1", "ERAC1U", "ERAC3U", "ERAC9U",
  "PL195"
  )

summaryplot_MLE <- ggplot2::ggplot(
  temp %>% dplyr::filter(
    Family != "Base", Family != "Imperfect",
    !Shorthand %in% removed_shorthands
  ), ggplot2::aes(
    x = Avg_Signum_Xmax_over_simulation,
    ymin = mle_025, y = mle_500, ymax = mle_975, label = Shorthand
  )
) + ggplot2::geom_hline(
  yintercept = 2.5, linetype = "dashed"
) + ggplot2::geom_point(
  data = temp %>% dplyr::filter(
    Family != "Base", Family != "Imperfect",
    !Shorthand %in% removed_shorthands
  ) %>% dplyr::select(-PFrag),
  alpha = 0.2
) + ggplot2::geom_point(
  data = temp %>% dplyr::filter(
    Family == "Base"
  ) %>% dplyr::select(-PFrag, -Family),
  alpha = 0.2, color = "red"
) + ggplot2::geom_linerange(
) + ggplot2::geom_linerange(
  data = temp %>% dplyr::filter(
    Family == "Base"
  ) %>% dplyr::select(-Family),
  color = "red"
) + ggrepel::geom_text_repel(
  size = 3.3, min.segment.length = 0.2, nudge_x = 0.05, nudge_y = -0.5
) + ggrepel::geom_text_repel(
  data = temp %>% dplyr::filter(
    Family == "Base"
  ) %>% dplyr::select(-Family),
  color = "red",
  size = 3.3, min.segment.length = 0.2, nudge_x = -0.05, nudge_y = 0.5
) + ggplot2::facet_grid(
  Family ~ PFrag
) + ggplot2::labs(
  x = "K", y = "MLE Exponent"
) + ggplot2::coord_cartesian(
  ylim = c(1, 5)
) + ggplotTheme()

summaryplot_KSMLE <- ggplot2::ggplot(
  temp %>% dplyr::filter(
    Family != "Base", Family != "Imperfect",
    !Shorthand %in% removed_shorthands
  ), ggplot2::aes(
    x = Avg_Signum_Xmax_over_simulation,
    ymin = ksmle_025, y = ksmle_500, ymax = ksmle_975, label = Shorthand
  )
) + ggplot2::geom_hline(
  yintercept = 2.5, linetype = "dashed"
) + ggplot2::geom_point(
  data = temp %>% dplyr::filter(
    Family != "Base", Family != "Imperfect",
    !Shorthand %in% removed_shorthands
  ) %>% dplyr::select(-PFrag),
  alpha = 0.2
) + ggplot2::geom_point(
  data = temp %>% dplyr::filter(
    Family == "Base"
  ) %>% dplyr::select(-PFrag, -Family),
  alpha = 0.2, color = "red"
) + ggplot2::geom_linerange(
) + ggplot2::geom_linerange(
  data = temp %>% dplyr::filter(
    Family == "Base"
  ) %>% dplyr::select(-Family),
  color = "red"
) + ggrepel::geom_text_repel(
  size = 3.3, min.segment.length = 0.2, nudge_x = 0.05, nudge_y = -0.5
) + ggrepel::geom_text_repel(
  data = temp %>% dplyr::filter(
    Family == "Base"
  ) %>% dplyr::select(-Family),
  color = "red",
  size = 3.3, min.segment.length = 0.2, nudge_x = -0.05, nudge_y = 0.5
) + ggplot2::facet_grid(
  Family ~ PFrag
) + ggplot2::labs(
  x = "K", y = "KS-MLE Exponent"
) + ggplot2::coord_cartesian(
  ylim = c(1, 5)
) + ggplotTheme()

ggplot2::ggsave(
  plot = summaryplot_MLE, filename = "SummaryPlot_MLE.png", path = directory,
  units = "cm", dpi = 320, width = 20.0, height = 18.0
)
ggplot2::ggsave(
  plot = summaryplot_KSMLE, filename = "SummaryPlot_KSMLE.png", path = directory,
  units = "cm", dpi = 320, width = 20.0, height = 18.0
)

