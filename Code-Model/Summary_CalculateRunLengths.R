# Calculate Run Lengths
# Borrowing from CoFr_Runlengths2_vik_20_06_03_rnotebook.Rmd

# Setup: #######################################################################

imagedpi <- 350#"retina"
imageunits <- "mm"
imageheight <- 297
imagewidth <- 420
imageheight_small <- 120
imagewidth_small <- 148

library(dplyr)
library(ggplot2)

load_smoothing <- FALSE

directory <- file.path('.')
data_folders <- dir(directory, pattern = "Data", full.names = TRUE)
analyses_folders <- dir(file.path(directory, "Analyses"),
                        pattern = "Data", full.names = TRUE)
imagepath <- file.path(directory, "Images2023-08-11")

data_folders_targets <- dir(data_folders[7])

Method_loaded <- (Methods <- c(
  "Basic Truncated MLE",
  "CSN + ZXX MLE and KS: Lower Bound"
))[1:2]

Method_conversion <- c(
  "Basic Truncated MLE" = 'MLE',
  "CSN + ZXX MLE and KS: Lower Bound" = 'KS-MLE'
)

chunking <- 5000

# > runif(1) * 10 ^ 8
# [1] 13692633
if (exists(".Random.seed")) {
  old_seed <- .Random.seed
}
set.seed(13692633)

# Functions: ###################################################################
# See also chunk 26 for functions to interpolate the probabilities.

ggplotTheme <- function(
    ..., base_size = 22
) {
  ggplot2::theme_bw(base_size) + ggplot2::theme(...)
}

# Load: ########################################################################
# From Figure_CreateSteadyStateExample6.R

files_analyses <- c(
  grep(
    x = dir(
      analyses_folders[7],
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

# Chunk 11: Init. Runlengths: ##################################################
# Create our first grouping of things into runlengths.
# The old chunk 11 expects res_phy to be a list of data frames.
# Here, we already have a singular data frame, so we need to make some changes.

runlengths <- res_phy %>% dplyr::group_by(
  file, kernel, PFrag, M, Simulation
) %>% dplyr::arrange(time) %>% dplyr::mutate(
  PFrag = as.numeric(gsub(pattern = "v0", replacement = "0.", PFrag)),
  change_xmax = xmax - dplyr::lag(xmax),
  Sgn_Max = sign(change_xmax),
  run = cumsum(ifelse(is.na(change_xmax), 1, change_xmax < 0))
) %>% dplyr::group_by(
  file, kernel, PFrag, M, Simulation, run
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

# Chunk 26: Interp. Functions: #################################################
# Create the interpolation functions.
# But we have a problem, the previous interpolation functions were "subtly"
# expecting us to be shattering, but we're not doing that anymore.
# This is most important in case 1, where we're not stuck with a single group.
# I don't see any obvious way of approximating the final size of the group
# immediately before fragmentation, so I don't think that we should.
# A possibility would be to go back and re-analyse the data. Two options:
#   look at the second largest group and use the sum
#   see if we have enough data to not have to thin by 10.
# At present though, I'm going to just go without the extra feature.

case_1 <- function(Timestep, Pop, PFr, Groups, Xmax) {
  # Note that this function is vectorized (hence lapply and indexing).
  event0 <- 10 * (Timestep - 1)
  # If timestep 10, then 100 events occurred by end of timestep, i.e. 91 : 100.
  events <- rep(10, length(event0))
  xmax <- ifelse(is.infinite(Xmax), Pop, Xmax)
  Probabilities <- lapply(seq_along(events),
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
  Probabilities <- lapply(seq_along(events_before_geo),
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

# Chunk 27: ####################################################################

# This calculation can be long, so we store it once we are done with it.
# It also is surprisingly unstable, hence the large numbers of error handling.
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
  save.image(file = "Smoothed_Chunks_Finished_2023-08-11.RData")
} else {
  load(file = "Smoothed_Chunks_Finished_2023-08-11.RData")
}

# Chunk 28: ####################################################################
# This chunk makes sure that signum (order param.) is kept updated, as well as
# creating a plottable data frame.
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
        kernel = d$kernel[i],
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

# runlengths_smoothed <- do.call(rbind, runlengths_smoothed_subchunks)
runlengths_smoothed <- dplyr::bind_rows(runlengths_smoothed_subchunks)
rm(runlengths_smoothed_subchunks)

# Chunk 29: ####################################################################
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

# Chunk 30: ####################################################################

ggplot2::ggplot(runlengths_smoothed_histvals,
                ggplot2::aes(
         x = Event_Number,
         weight = FragProb
       )
) + ggplot2::geom_histogram(
  binwidth = function(x) {max(c(floor(diff(range(x)) / 100), 1))}
) + ggplot2::facet_wrap(
  PFrag ~ kernel,
  nrow = length(unique(runlengths_smoothed_histvals$PFrag)),
  ncol = length(unique(runlengths_smoothed_histvals$kernel)),
  scales = "free"
) + labs(x = "Physical Time Run Lengths", y = "Weighted Counts") + ggplotTheme()


# Chunk 31: ####################################################################
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

save(summarystats2, file = paste0("summarystats2-",
                                  "2023-08-14"
                                  #Sys.Date(),
                                  ".RData"))

