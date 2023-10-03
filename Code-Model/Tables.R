# Warning: If the libraries are not available, this script attempts to install!
# Settings: ###################################################################

tableNumber <- 2

# Location: Code-Model
# https://stackoverflow.com/a/36276269
directory <- getSrcDirectory(function(dummy) {dummy})

source(file.path(directory, "Figure_Functions.R"))

#directory <- '.'
# Location: Root.
directory <- file.path(directory, "..")

load_smoothing <- FALSE

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

Method_loaded <- (Methods <- c(
  "Basic Truncated MLE",
  "CSN + ZXX MLE and KS: Lower Bound"
))[1:2]

Method_conversion <- c(
  "Basic Truncated MLE" = 'MLE',
  "CSN + ZXX MLE and KS: Lower Bound" = 'KS-MLE'
)


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
files_analyses <- c(
  grep(
    x = dir(
      file.path(directory, "Analyses",
                if (tableNumber == 1) {
                  "Data20-10-10"
                } else if (tableNumber == 2) {
                  "Data20-08-25"
                } else {
                  stop("tableNumber not recognized.")
                }),
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

      # kernel = substring(file, 29, 30)
      #     n1 = n_1
      # ) %>% dplyr::select(
      #     -n_1
      # ) %>% dplyr::rename(
      #   N = clusters
    )

  if("clusters" %in% colnames(result_storage_subset)) {
    result_storage_subset <- result_storage_subset %>% dplyr::rename(
      N = clusters
    )
  }


  if(!"N" %in% colnames(result_storage_subset)) {
    # Now rare old standard from before PFrag exploration
    # Need to load the original files to extract.

    file_candidates <- dir(
      file.path(
        directory, if (tableNumber == 1) {
          "Data20-10-10"
        } else if (tableNumber == 2) {
          "Data20-08-25"
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
      by = "time"
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

    if (chunk == 1) {
      result_storage_subset_out <- rss_exponent %>% dplyr::full_join(
        rss_xmin, by = c('time', 'xmax', 'M', 'n_1', 'PFrag',
                         'N', 'Simulation',
                         'error', 'Kernel', "Variant", 'kernel', 'file'),
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
        rss_xmin, by = c('time', 'xmax', 'M', 'n_1', 'PFrag',
                         'N', 'Simulation',
                         'error', 'Kernel', "Variant", 'kernel', 'file'),
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

# Data for Tables: ############################################################

### Cyclicity: ################################################################

##### 1st Approx. Method: #####################################################
# This works well as a first approximation, but fairs poorly when the system
# should be spending a lot of time waiting for fragmentation to happen.
# It also neglects that our simulations have 10 time steps to a row (to save
# on time (write-frequency) and memory (storage)).
cyclicity_1st <- res_phy %>% dplyr::group_by(
  file
) %>% dplyr::mutate(
  signumXmax = sign(xmax - lag(xmax))
) %>% dplyr::ungroup(
) %>% dplyr::group_by(
  Kernel, Variant, PFrag, M
) %>% dplyr::summarise(
  meanK = mean(signumXmax, na.rm = TRUE),
  .groups = "drop"
)

##### 2nd Approx. Method: #####################################################
# Iterating from Summary_CalculateRunLengths.R.
# Note in particular that we are not always, or even necessarilly usually,
# in the shattering case here.

####### See Chunk 11: ########################################################
runlengths <- res_phy %>% dplyr::group_by(
  file, Kernel, Variant, kernel, PFrag, M, Simulation
) %>% dplyr::arrange(time) %>% dplyr::mutate(
  PFrag = as.numeric(gsub(pattern = "v0", replacement = "0.", PFrag)),
  change_xmax = xmax - dplyr::lag(xmax),
  Sgn_Max = sign(change_xmax),
  run = cumsum(ifelse(is.na(change_xmax), 1, change_xmax < 0))
) %>% dplyr::group_by(
  file, Kernel, Variant, kernel, PFrag, M, Simulation, run
) %>% dplyr::summarise(
  # Outside of the "shattering" context, it's hard to know what exactly got
  # fragmented (since the shards could end up anywhere).
  Last_xmax = dplyr::last(xmax),
  First_xmax = dplyr::first(xmax),
  Number_of_groups = dplyr::last(N),
  Sum_Signum_xmax = sum(Sgn_Max, na.rm = TRUE),
  Length = dplyr::n(),
  First_Time = dplyr::first(time),
  Last_Time = dplyr::last(time),
  .groups = "drop"
)

######## See Chunk 26: ########################################################
case_1 <- function(Timestep, Pop, PFr, Groups, Xmax) {
  # Note that this function is vectorized (hence lapply and indexing).
  event0 <- 10 * (Timestep - 1)
  # If timestep 10, then 100 events occurred by end of timestep, i.e. 91 : 100.
  events <- rep(10, length(event0))
  xmax <- ifelse(is.infinite(Xmax), Pop, Xmax)
  Probabilities <- lapply(
    seq_along(events),
    function(i, x, p, M, size) {
      # Form a vector of length events_before_geo
      retVal <- rep(0, x[i])
      # Populate it with probabilities of fragmentation before waiting.
      retVal <- dgeom(0:9, p[i] * size[i] / M[i])
      # Normalise
      retVal <- retVal / sum(retVal)
      retVal
    }, x = events, p = PFr, M = Pop, size = xmax)

  return(
    lapply(seq_along(Probabilities), function(ps, e0, p) {
      data.frame(
        Event_Number = e0[ps] + 1:length(p[[ps]]),
        Probabilities = p[[ps]]
      )}, e0 = event0, p = Probabilities)
  )
}

case_2 <- function(Timestep, Pop, PFr, Groups, Xmax) {
  # Note that this function is vectorized (hence lapply and indexing).
  event0 <- 10 * (Timestep - 1)
  # If timestep 10, then 100 events occurred by end of timestep, i.e. 91 : 100.
  events_before_geo <- Groups - 1
  xmax <- ifelse(is.infinite(Xmax), Pop, Xmax)
  Probabilities <- lapply(
    seq_along(events_before_geo),
    function(i, x, p, M, size) {
      # Form a vector of length events_before_geo
      retVal <- rep(0, x[i])
      # Populate it with probabilities of fragmentation before waiting.
      if (x[i] > 0) {
        retVal[1] <- p[i] * size[i] / M[i]
        if (x[i] > 1)
          for (j in 2:x[i]) {
            # Note the (i - 1) is the number of failures == coal. events.
            retVal[j] <- exp(sum(log(1 - retVal[1:(j - 1)]))) * p[i] * (
              size[i] + (M[i] - size[i]) / x[i] * (j - 1)
            ) / M[i]
          }
      }
      # Divide remaining probability mass  amongst a geometric distribution
      # with probability of success of PFr.
      p_remaining <- 1 - sum(retVal)
      # Determine how far to go by making sure we have --.-% coverage.
      mark_99 <- qgeom(0.9, p[i])
      # The remaining mass is divided up. Note dgeom 1st arg = number of failures.
      retVal <- c(retVal, p_remaining * dgeom(0:mark_99, p[i]))
      retVal
    }, x = events_before_geo, p = PFr, M = Pop, size = xmax)

  return(
    lapply(seq_along(Probabilities), function(ps, e0, p) {
      data.frame(
        Event_Number = e0[ps] + 1:length(p[[ps]]),
        Probabilities = p[[ps]]
      )}, e0 = event0, p = Probabilities)
  )
}

######## See Chunk 27: ########################################################
if (!load_smoothing) {
  runlengths_smoothed_chunks <- lapply(
    seq(from = 1, to = nrow(runlengths), by = 1000),
    function(i, r) {

      print(paste("i : ", i, ":", Sys.time()))

      tryCatch(
        r[i:(min(c(i + 999, nrow(r)))), ] %>% mutate(
          smoothed =
            # smoothed is an entry composed of a dataframe
            # with Run_Length and Probability_Mass entries
            case_when(
              Number_of_groups <= 10 ~
                case_2(Length, M, PFrag, Number_of_groups, Last_xmax),

              # In order execution, so we are in the > 10 case.
              TRUE ~
                case_1(Length, M, PFrag, Number_of_groups, Last_xmax)
            )
        ),
        error = function(e) {

          message(paste("error :", i, ":", e))

          lapply(
            seq(from = i, to = i + 1000 - 1, by = 100),
            function(i2, r, i) {
              tryCatch(
                r[i2:(min(c(i2 + 99, nrow(r)))), ] %>% mutate(
                  smoothed =
                    # smoothed is an entry composed of a dataframe
                    # with Run_Length and Probability_Mass entries
                    case_when(
                      Number_of_groups <= 10 ~
                        case_2(Length, M, PFrag, Number_of_groups, Last_xmax),

                      # In order execution, so we are in the > 10 case.
                      TRUE ~
                        case_1(Length, M, PFrag, Number_of_groups, Last_xmax)
                    )
                ),

                error = function(e2) {
                  message(paste("error check :", i, ":", i2, ": error :", e2))
                  NA
                },
                warning = function(w2) {
                  message(paste("error check :", i, ":", i2, ": warning :", w2))
                  NA
                }

              )
            }, r = r, i = i)
        },
        warning = function(w) {

          message(paste("warning :", i, ":", w))

          lapply(
            seq(from = i, to = i + 1000 - 1, by = 100),
            function(i2, r, i) {
              tryCatch(
                r[i2:(min(c(i2 + 99, nrow(r)))), ] %>% mutate(
                  smoothed =
                    # smoothed is an entry composed of a dataframe
                    # with Run_Length and Probability_Mass entries
                    case_when(
                      Number_of_groups <= 10 ~
                        case_2(Length, M, PFrag, Number_of_groups, Last_xmax),

                      # In order execution, so we are in the > 10 case.
                      TRUE ~
                        case_1(Length, M, PFrag, Number_of_groups, Last_xmax)
                    )
                ),
                error = function(e2) {
                  message(paste("warning check :", i, ":", i2, ": error :", e2))
                  NA
                },
                warning = function(w2) {
                  message(paste("warning check :", i, ":", i2, ": warning :", w2))
                  NA
                }
              )
            }, r = r, i = i)
        }
      )
    },
    r = runlengths
  )

  if (tableNumber == 1) {
    save.image(file = "Smoothed_Chunks_Finished_2023-09-11.RData")
  } else if (tableNumber == 2) {
    save.image(file = "Smoothed_Chunks_Finished_2023-09-18.RData")
  }
} else {
  if (tableNumber == 1) {
    load(file = "Smoothed_Chunks_Finished_2023-09-11.RData")
  } else if (tableNumber == 2) {
    load(file = "Smoothed_Chunks_Finished_2023-09-18.RData")
  }
}


######## See Chunk 28: ########################################################
runlengths_smoothed_subchunks <- lapply(
  runlengths_smoothed_chunks, function(d) {
    # Need to apply the appropriate descriptors to each sub-dataframe from main.
    dplyr::bind_rows(lapply(seq_along(d$smoothed), function(i, d) {
      d$smoothed[[i]] %>% dplyr::mutate(
        case = case_when(
          d$Number_of_groups[i] <= 10 ~ "2",
          # In order execution, so we are in the > 10 case.
          TRUE ~ "1"
        ),
        run = d$run[i],
        Simulation = d$Simulation[i],
        kernel = ifelse(!is.na(d$kernel[i]), as.character(d$kernel[i]),
                        paste("MixedKernel", d$Kernel[i], d$Variant[i],
                              sep = "_")),
        PFrag = d$PFrag[i],
        M = d$M[i],
        Signum = d$Sum_Signum_xmax[i] * 10 + case_when(
          case == "2" ~ c(
            rep(1, min(d$Number_of_groups[i] - 1, 9)),
            rep(0, max(nrow(d$smoothed[[i]]) - d$Number_of_groups[i] + 1, 1))
          ),
          case == "1" ~ 0,
          TRUE ~ NA_real_
        )
      )
    }, d = d))
  })

rm(runlengths_smoothed_chunks)

runlengths_smoothed <- dplyr::bind_rows(runlengths_smoothed_subchunks)
rm(runlengths_smoothed_subchunks)

######## See Chunk 29: ########################################################
# Same as above, but part 2.
# Note we also remove some invalid parts to make sure we have complete runs.
runlengths_smoothed_histvals <- runlengths_smoothed %>% group_by(
  PFrag, M, Simulation, kernel
) %>% filter(
  run != first(run),
  run != last(run)
) %>% ungroup() %>% group_by(
  PFrag, M, Event_Number, kernel
) %>% summarise(FragProb = sum(Probabilities), .groups = "drop")

# Fix for geom_area:
runlengths_smoothed_histvals_lead <- runlengths_smoothed_histvals %>% group_by(
  PFrag, M, kernel
) %>% filter(
  (lead(Event_Number) != Event_Number + 1) | is.na(lead(Event_Number))
) %>% mutate(Event_Number = Event_Number + 1,
             FragProb = 0)
runlengths_smoothed_histvals_lag <- runlengths_smoothed_histvals %>% group_by(
  PFrag, M, kernel
) %>% filter(
  (lag(Event_Number) != Event_Number - 1) | is.na(lag(Event_Number))
) %>% mutate(Event_Number = Event_Number - 1,
             FragProb = 0)
runlengths_smoothed_histvals <- runlengths_smoothed_histvals %>% rbind(
  runlengths_smoothed_histvals_lead,
  runlengths_smoothed_histvals_lag
)


######## See Chunk 31: ########################################################
summarystats2 <- runlengths_smoothed_histvals %>% group_by(
  PFrag, M, kernel
) %>% arrange(
  Event_Number
) %>% mutate(
  NormalisedFragProb = FragProb / sum(FragProb, na.rm = TRUE)#,
  # Cumul_Norm_FragProb = cumsum(ifelse(is.na(NormalisedFragProb), 0, NormalisedFragProb))
)

summarystats2_signum <- runlengths_smoothed %>% group_by(
  PFrag, M, kernel
) %>% select(
  -case
) %>% mutate(
  Event_Prob = Event_Number * Probabilities,
  Signum_Prob = Signum * Probabilities,
  # Previously, we had:
  # Avg_Signum_Xmax_over_simulation =
  # # ("Positive" steps - number of runs + number of unfinished runs) / number of steps
  # # This might technically over-penalise since each step is 10 long,
  # # which suggests dividing the lengths by 10 to normalise...
  # (sum(Sum_Signum_xmax) - length(n)/10 + length(unique(Simulation))/10)  / sum(n),
  # Now we have steps at event resolution, but events are probabilistic.
) %>% summarise(
  Avg_Signum_Xmax_over_simulation =
    (sum(Signum_Prob, na.rm = TRUE) # Probabilistic number of steps
     - (# Number of fragmentations of x_max = (runs - unfinished runs)
       length(unique(run))
       - length(unique(Simulation))
     )) / (# Divided by the total number of events.
       sum(Event_Prob, na.rm = TRUE)
     ), .groups = "drop"
) %>% mutate(
  r = M * PFrag / (1 - PFrag)
)

summarystats2 <- summarystats2 %>% group_by(
  PFrag, M, kernel
) %>% summarise(
  Mean = sum(NormalisedFragProb * Event_Number, na.rm = TRUE),
  Var = sum(NormalisedFragProb * (Event_Number - Mean)^2, na.rm = TRUE),
  Skew = sum(NormalisedFragProb * ((Event_Number - Mean) / sqrt(Var))^3, na.rm = TRUE),
  Kurt = sum(NormalisedFragProb * ((Event_Number - Mean) / sqrt(Var))^4, na.rm = TRUE),
  .groups = "drop"
) %>% mutate(
  r = M * PFrag / (1 - PFrag)
)

summarystats2 <- summarystats2 %>% left_join(
  summarystats2_signum,
  by = c("r", "M", "PFrag", "kernel")
)

save(summarystats2,
     file = paste0("summarystats2-",
                   if (tableNumber == 1) {
                     "2023-09-11"
                   } else if (tableNumber == 2) {
                     "2023-09-18"
                   },
                   #Sys.Date(),
                   ".RData"))

### Power-law exponent estimators: ############################################


