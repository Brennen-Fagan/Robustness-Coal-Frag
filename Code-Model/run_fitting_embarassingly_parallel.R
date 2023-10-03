for (i in 394:1) {
  # Previously, we tried to use the parallel, doParallel, and foreach packages
# on YARCC. This did not work due to problems with ports.
# Now, we are going to try a more embarassingly parallel approach.
# Unfortunately, we will be sacrificing coordination between the nodes.
# In particular, if a node gets a particularly hard set of fits in a row,
# which is likely, e.g. high population, we will experience significant
# slow-down.

# command line Arguments:
#    File Number
#    Start Row
#    End Row
print("Reading Args")
cargs <- c(i, 1, 20000) #as.numeric(commandArgs(TRUE)) #c(3, 1, 20000)
print(cargs)

print("Reading Settings")
methods_to_use <- c(
  #1, # CSN's method
  #2, # Bootstrapped CSN
  3, # Truncated MLE
  #4, # MaxEnt on all Data
  #5, # MaxEnt with "optimised" boundaries
  6#, # KS to choose lower bound, MLE
  #7#, # KS to choose lower bound, Order Statistics
  #8 # KS to choose upper bound, Order Statistics
)

debug_skip_flux <- 1

methods_to_use_flux <- c(
  #1, # CSN's method
  #2, # Bootstrapped CSN
  3, # Truncated MLE
  #4, # MaxEnt on all Data
  #5, # MaxEnt with "optimised" boundaries
  6#, # KS to choose lower bound, MLE
  #7#, # KS to choose lower bound, Order Statistics
  #8 # KS to choose upper bound, Order Statistics

  #TODO Figure out what the differene is between method 1 and the
  #supplementary method to 2 (the former fits xmins by decreasing
  #xmax, while the latter just calls the default methods).
  #TODO Include methods for Lognormal
)

smoothing <- c(1, 10, 100, 1000)

estimate_rows_per_file <- cargs[3] - cargs[2] + 1
data_folder <- 'Analyses'
image_folder <- 'Images'

sparse_mode <- 1
overwrite <- 0

#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#D55E00", "#F0E442", "#0072B2", "#CC79A7")

# #https://stackoverflow.com/a/36276269
# directory <- getSrcDirectory(function(dummy) {dummy})
directory <- '.'

# Adjust libraries: ###########################################################
print("Loading libraries")
librarypath <- file.path(directory, "Rlibs")
if(!dir.exists(librarypath)){
  dir.create(librarypath, showWarnings = FALSE)
}
.libPaths(c(librarypath, .libPaths()))

packages <- c("ggplot2", "poweRlaw", "DEoptim", "R.matlab")

for(package in packages) {
  if(!require(package, character.only = TRUE)) {
    install.packages(package, lib = librarypath,
                     repos = 'https://cloud.r-project.org',
                     dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

#Methods:######################################################################
print("Reading Methods")
#Currently, we have a few of 'worthwhile' fitting methods.
#Unfortunately, none of our methods are perfect, each with
#different assumptions and sensitivities.

###1 CSN Classic, as implemented in Gillespie's poweRlaw:######################
#Pro: The standard approach, does well when needing to only
#look at the tail, gives good fits for most exponents.
#Con: The KS test can 'get lucky' and indicate a close, but
#not *good* fit to the data, leading to rejection of true
#power laws. This method also does not perform well if the
#power law's exponent is close to 1 and it is biased in the
#case of truncated distributions.
#Note: Extracted in full from PowerFitting_Experiment.R.
fit_poweRlaw_dis <- function(dataset, calculate_p = 0) {
  fit <- poweRlaw::displ$new(dataset)
  xmaxs <- c(Inf, sort(unique(dataset), decreasing = T))
  index <- 0
  errorflag <- 1
  while (errorflag) {
    index <- index + 1
    tryCatch(
      {
        fit$setXmin(poweRlaw::estimate_xmin(fit, xmax = xmaxs[index]))
        errorflag = 0
      }, error = function(c){
        #message(c) #Message issues?
      }
    )
    if (errorflag && index == length(xmaxs)){
      return(list(xmin = 0,
                  xmax = Inf,
                  exponent = Inf))
    }
  }

  if (calculate_p) {
    p = poweRlaw::bootstrap_p(fit, xmax = Inf)
    return(list(xmin = fit$getXmin(),
                exponent = fit$getPars(),
                bootstrap = p))
  }
  return(list(xmin = fit$getXmin(),
              xmax = Inf,
              exponent = fit$getPars()))
}

###2 CSN Bootstrap:############################################################
#Pro: Simple extension of CSN that in theory addresses the
#danger of 'getting lucky' with the KS test and removes any
#lack of independence amongst our data. Could also be used
#to extend some of the other methods mentioned here.
#Con: Bootstrapping potentially introduces bias in the data.
#Furthermore, it is not obvious how to extend some of the
#ideas of fit comparison in a bootstrapped context without
#suffering large hits to performance time.
#Note: extracted from R_Full_P_L_Testing.R.

#Supplementary Function:
GenerateManySamples <- function(Data,
                                TotalSamplesPerIteration,
                                TotalIterations) {
  #Assuming Data is a list of same size matrices,
  #where each row of a matrix is a sample. If
  #Data is instead a vector (as in the trace case)
  #then we take care of it in the main function...

  #First, determine the source of each data point.
  SampleGenerated <- #matrix(
    sample(x = length(Data),
           size = TotalSamplesPerIteration * TotalIterations,
           replace = TRUE)
  #nrow = TotalSamplesPerIteration,
  #ncol = TotalIterations
  #)
  RetVal <- rep(NA, length(SampleGenerated))

  #Determine the number of accesses to each source
  #and which sources are accessed.
  uniqueSources <- sort(unique(as.vector(SampleGenerated)))
  SourceAccesses <- table(SampleGenerated)

  for(Source in 1:length(uniqueSources)){
    #Open Source
    temp <- Data[[Source]]
    #Calculate Entries, in case of list.
    entries <- sample(length(temp),
                      replace = TRUE,
                      size = SourceAccesses[Source])
    #Calculate Rows of Entries that we will use per access.
    lengths <- sapply(entries, function(x){dim(temp[[x]])[1]})
    rows <- ceiling(runif(n = length(entries), min = 0, max = lengths))
    #Reweight and Select a column to return unif@rand group.
    #(Note: the reweighting is to allow the unif@rand selection.)
    SampleGenerated_part <- sapply(seq_along(entries), function(x){
      nonzeros <- which(temp[[entries[x]]][rows[x], ] > 0)
      return(
        nonzeros[
          sample(length(nonzeros), size = 1,
                 prob = temp[[entries[x]]][rows[x], nonzeros]/
                   sum(temp[[entries[x]]][rows[x], ]))
          ]
      )
    })
    #We then need to !CAREFULLY! put the draws back in the places
    #that they came from.
    RetVal[SampleGenerated == uniqueSources[Source]] <-
      SampleGenerated_part
  }

  return(matrix(RetVal,
                nrow = TotalSamplesPerIteration,
                ncol = TotalIterations))
}

#Supplementary Function:
FitDistribution <- function(DataSample,
                            poweRlawDist = poweRlaw::displ,
                            CalculatePValue = 0,
                            Censor = 0, NumThreads = 1) {
  Model <- poweRlawDist$new(DataSample)
  Model$setXmin(poweRlaw::estimate_xmin(Model, xmax = Inf))
  Fit <- if(CalculatePValue)
    tryCatch(
      list(poweRlaw::bootstrap_p(Model, threads = NumThreads,
                                 xmax = if(Censor)
                                   max(DataSample)
                                 else
                                   Inf
      )
      ),
      error = function(cond){
        #Fit not possible.
        temp <- list(NULL)
        temp[[1]]$p = NA
        return(temp)
      })
  else NULL
  return(c(Model$getPars(), Model$getXmin(),
           poweRlaw::get_ntail(Model) / poweRlaw::get_n(Model),
           Fit[[1]]$p)
  )
}

PowerLawFitSimple <- function(Data,
                              NumberOfSamples = 1000,
                              NumberOfBootstraps = 1000,
                              CensorTail = 0,
                              Print = 1,
                              NumberOfThreads = 1){
  #Setup storage
  PowerLawParams <- matrix(nrow = NumberOfBootstraps, ncol = 4)

  #Begin Testing: Parameter distributions
  if(Print) {
    print("Creating Parameter Distributions")
    StartTime <- Sys.time()
    print(c("Start time: ", StartTime))
    itertime <- Sys.time()
  }

  #Determine if we need to invoke a more complicated
  #bootstrapping procedure in order to generate the
  #samples or not. If the data is a whole matrix or
  #a collection of matrices, we expect it in a list.
  #We then use a supplementary function to generate
  #our data sets. Otherwise, we use a simpler
  #procedure.
  if (is.list(Data)) {
    Samples <- GenerateManySamples(Data = Data,
                                   SamplesPerIteration = NumberOfSamples,
                                   TotalIterations = NumberOfBootstraps)
  } else {
    Samples <- matrix(
      data = sample(Data,
                    size = NumberOfSamples * NumberOfBootstraps,
                    replace = T),
      nrow = NumberOfBootstraps,
      ncol = NumberOfSamples
    )
  }


  for(BootstrapNumber in 1:NumberOfBootstraps) {
    if(Print) {
      previtertime <- itertime
      itertime <- Sys.time()
      print(c(BootstrapNumber,
              itertime - previtertime,
              units(itertime - previtertime)))
    }
    #Take Sample
    Sample <- Samples[BootstrapNumber, ]

    #Returns in order of Parameters, Xmins, and TailProportion
    PowerLawParams[BootstrapNumber, ] <- tryCatch(
      FitDistribution(Sample, poweRlaw::displ,
                      CalculatePValue = TRUE, Censor = CensorTail,
                      NumThreads = NumberOfThreads),
      error = function(cond) {
        if(Print) {
          message(paste("poweRlaw failed to fit, skipping iteration.",
                        "Original message:"))
          message(cond)
        }
        return(c(NA, NA, NA, NA))
      })
    if(is.na({PowerLawParams[BootstrapNumber, 1]})) next

  }
  if(Print) {
    print(c("Finished at ", Sys.time()))
    print(Sys.time()-StartTime)
  }

  return( list(
    xmin = median(PowerLawParams[, 2]),
    xmax = Inf,
    exponent = median(PowerLawParams[, 1])
  ))
}

###3 Basic MLE:################################################################
#Pro: Simple ideas, easily adapts to the truncated environment
#or the infinite environment with decent results, comparable to 1.
#Con: It tries to fit a power law to the entire data set, and
#assumes that the range of the data is the range of the
#(truncated) power law. Hence, this does *not* adapt outside of
#its basic environment.
#Note: Extracted from PowerFitting_Experiment.R.
fit_displ_trunc_max_like <- function(dataset) {
  #If entire data set is to be power law fitted,
  #then estimator for xmin is just the minimum.
  xmax <- max(dataset)
  xmin <- min(dataset)

  #Solve equation 2 in ZXX.
  exponent <- uniroot(function(exp, d, xmin, xmax) {
    x <- xmin:xmax
    lnmn <- sum(log(d))/length(d)
    return(lnmn -
             sum(x ^ (-exp) * log(x)) /
             sum(x ^ (-exp))
    )
  }, interval = c(0, 100),
  d = dataset, xmin = xmin, xmax = xmax)$root

  return(list(xmin = xmin, xmax = xmax, exponent = exponent))
}

###4 Maximum Entropy:##########################################################
#Pro: The method assumes the minimum amount of information to
#create the given power law distribution, and should be able
#to determine by itself the upper and lower bounds of the
#power law data by simple maximization of the entropy.
#Con: I am not an expert in this field, so my interpretation
#may be incorrect. In addition, when allowed to vary the
#bounds, it appears that the resulting bounds are fairly
#conservative even when all the data *are* power law
#distributed. When the data are not all power law distributed,
#then it does not successfully retrieve the power law
#necessarily.
#Note: Extracted in full from PowerFitting_Experiment.R.
trunc_maximize_entropy <- function(data, fixed_edges = 0,
                                   bsval = 0,
                                   steptolval = 200,
                                   itermaxval = 200,
                                   traceval = 0) {
  unique_data <- sort(unique(data))

  solve_calced_chi_equal_empirical_chi <- function(data, N0, N) {
    echi <- mean(log(data[N0 <= data & data <= N]))

    cchi <- function(z, N0, N) {
      sum((N0 : N) ^ -z * log(N0 : N)) / sum((N0 : N) ^ -z)
    }

    #Minimize difference.
    optim(1, function(z) {
      abs(echi - cchi(z, N0, N))
    },
    lower = -10, upper = 10,
    method = "Brent")$par
  }

  if (fixed_edges) {
    return(list(
      exponent = solve_calced_chi_equal_empirical_chi(
        data, unique_data[1], unique_data[length(unique_data)]
      ),
      xmin = unique_data[1],
      xmax = unique_data[length(unique_data)]
    )
    )
  } else {
    pars <- DEoptim::DEoptim(
      function(x) {

        N0 <- x[1]
        N <- x[2]

        if (length(data[N0 <= data & data <= N]) == 0) {
          return(Inf)
        }

        z <- solve_calced_chi_equal_empirical_chi(data, N0, N)
        cchi <- (sum((N0 : N) ^ -z * log(N0 : N))
                 / sum((N0 : N) ^ -z))
        #Note: optim minimizes, so take the negative...
        res <- -(z * cchi + log(sum((N0 : N) ^ -z)))
        if (is.nan(res)) {
          print(paste(N0, N, z, cchi))
        }

        return(res)
      },
      #Enforce data as boundaries, else it will try to
      #include as much (non-existent) data as it can.
      lower = c(unique_data[1], unique_data[2]),
      upper = c(unique_data[length(unique_data) - 1],
                unique_data[length(unique_data)]),
      fnMap = function(x) c(floor(x[1]), floor(x[2])),
      control = list(trace = 10, NP = 20,
                     bs = bsval, steptol = steptolval,
                     itermax = itermaxval, trace = traceval)
    )

    z <-
      solve_calced_chi_equal_empirical_chi(
        data, pars$optim$bestmem[1], pars$optim$bestmem[2])

    return(list(
      exponent = z,
      xmin = pars$optim$bestmem[1],
      xmax = pars$optim$bestmem[2],
      entropy = pars$optim$bestval)
    )
  }
}

###5 Order Statistics:#########################################################
#Pro: When provided the appropriate tail, these methods appear
#successful.
#Con: Need to know which tail to provide to begin with, my best
#implementation still relies on the KS test, and the estimate
#for the opposite bound appears to be quite conservative.
#Additionally, relies on a control parameter that determines
#how much of the data is provided for the fitting. While we use
#1/10 or 9/10 (determined by tail), there is not an obvious
#choice, especially with mixed distributions. Furthermore, if
#both tails of the distribution are not power law, this, and the
#other methods all, seem to have problems.
#Note: Extracted in full from PowerFitting_Experiment.R.
fit_trunc_CSN_ZXX <- function(dataset, KS = "Below",
                              otherside = "Order", silent = 1) {
  #Mix CSN's idea of KS testing for best position with
  #ZXX's idea of MLEs, either through order statistics
  #or through the use of the min or max of the dataset.
  uniq = sort(unique(dataset))

  #Setup, replicating previous functions in this file.
  if ( KS == "Below" ) {
    uniq = uniq[1 : (length(uniq) - 1)]
    if (otherside == "Basic") {
      xmax = max(dataset)

      uniroot_func <- function(exp, d, xmin, xmax) {
        x <- xmin:xmax
        lnmn <- sum(log(d))/length(d)
        return(lnmn -
                 sum(x ^ (-exp) * log(x)) /
                 sum(x ^ (-exp))
        )
      }
    } else if (otherside == "Order") {
      r <- floor(9 / 10 * length(dataset))
      x_1_r <- sort(sort(dataset, partial = r)[1 : r], decreasing = F)

      #Derived.
      NegLogLikelihood <- function(
        x, #= c(exponent, maxim),
        minim = min(dataset), #<- and below will be replaced during KS search
        n_ = length(dataset),
        order_statistics = x_1_r,
        r_ = length(order_statistics),
        maxim = NULL #<- dummy argument, not used.
      ) {
        return(-(#Can ignore the constant
          - n_ * log(x[2] + sum((minim : order_statistics[r_]) ^ -x[1]))
          - (x[1]) * sum( log( order_statistics[1 : (r_ - 1)] ) )
          + (n_ - r_ + 1) * log(x[2])
          + log( (x[2] ^ -1 * order_statistics[r_] ^ -x[1] + 1) ^ (n_ - r_ + 1) - 1)
        ))
      }

      uniroot_func <- function(m, exp, int, xr){
        int - sum(((xr + 1):m) ^ (-exp))
      }
    } else {
      return("Unknown otherside argument")
    }
  } else if ( KS == "Above" ) {
    uniq = rev(uniq[2 : length(uniq)])
    if (otherside == "Basic") {
      xmin = min(dataset)

      uniroot_func <- function(exp, d, xmin, xmax) {
        x <- xmin:xmax
        lnmn <- sum(log(d))/length(d)
        return(lnmn -
                 sum(x ^ (-exp) * log(x)) /
                 sum(x ^ (-exp))
        )
      }
    } else if (otherside == "Order") {
      r <- floor(1 / 10 * length(dataset))
      x_1_r <- sort(-sort(-dataset, partial = r)[1 : r], decreasing = T)

      #Order statistic based likelihood:
      #Equation before Equation 3
      NegLogLikelihood <- function(
        x, #= c(exponent, minim),
        maxim = max(dataset), #<- and below will be replaced during KS search
        n_ = length(dataset),
        order_statistics = x_1_r,
        r_ = length(order_statistics),
        minim = NULL #<- dummy argument, not used.
      ) {
        return(-(
          - n_ * log(x[2] + sum((order_statistics[r_] : maxim) ^ -x[1]))
          - (x[1]) * sum( log( order_statistics[1 : (r_ - 1)] ) )
          + (n_ - r_ + 1) * log(x[2])
          + log((x[2] ^ -1 * order_statistics[r_] ^ -x[1] + 1) ^ (n_ - r_ + 1) - 1)
        ))
      }

      uniroot_func <- function(m, exp, int, xr){
        int - sum((m : (xr - 1)) ^ (-exp))
      }

    } else {
      return("Unknown otherside argument")
    }
  } else {
    return("Unknown KS Argument")
  }

  Best_KSd <- Inf
  Best_min <- 0
  Best_max <- Inf
  Best_exp <- Inf

  for (ksside in uniq) {
    #Truncate Dataset
    if (KS == "Below") {
      xmin <- ksside
    } else if (KS == "Above") {
      xmax <- ksside
    }

    #Perform MLE calculation(s)
    if (otherside == "Order") {
      if (KS == "Below") {
        temp_dataset <- dataset[xmin <= dataset]
        #Reduce r as appropriate:
        rprime <- max(floor(9 / 10 * length(temp_dataset)), 1)
        #Need to recalculate since we are reducing the appropriate side.
        #(Likely can shortcut somewhat by doing more than necessary
        # initially, and then adjusting "downward" repeatedly...)
        x_1_r <- sort(sort(temp_dataset, partial = rprime)[1 : rprime], decreasing = F)
        xmax <- NULL

      } else if (KS == "Above") {
        temp_dataset <- dataset[dataset <= xmax]
        #Reduce r as appropriate:
        rprime <- max(floor(1 / 10 * length(temp_dataset)), 1)
        #See above comment.
        x_1_r <- sort(-sort(-temp_dataset, partial = rprime)[1 : rprime], decreasing = T)
        xmin <- NULL

      }
      if (x_1_r[1] == x_1_r[rprime]) next
      if (max(temp_dataset) <= 1) next

      #Find acceptable starting value.
      start = 1
      arun <- NegLogLikelihood(
        c(start, max(temp_dataset)),
        maxim = xmax, minim = xmin,
        r_ = rprime, n_ = length(temp_dataset),
        order_statistics = x_1_r)
      while( (is.infinite(arun) || is.na(arun) || is.nan(arun)) && start < 10 ) {
        start = start + 1
        arun <- NegLogLikelihood(
          c(start, max(temp_dataset)),
          maxim = xmax, minim = xmin,
          r_ = rprime, n_ = length(temp_dataset),
          order_statistics = x_1_r)
      }
      if (start >= 10) next

      roots <- constrOptim(c(start, max(temp_dataset)), NegLogLikelihood, NULL,
                           matrix(c(1,0,0,1), nrow = 2), c(0,0),
                           maxim = xmax, minim = xmin,
                           r_ = rprime, n_ = length(temp_dataset),
                           order_statistics = x_1_r)$par

      xother <- tryCatch(
        uniroot(uniroot_func, interval = c(min(temp_dataset), max(temp_dataset)),
                exp = roots[1], int = roots[2], xr = x_1_r[rprime])$root,
        error = function(e) {
          if (attributes(e)$class[1] == "simpleError") {
            #Both values likely on same side.
            #Return the interval value that
            #yields a result closer to zero.
            if ( abs( uniroot_func(min(temp_dataset), exp = roots[1],
                                   int = roots[2], xr = x_1_r[rprime]) ) <
                 abs( uniroot_func(max(temp_dataset), exp = roots[1],
                                   int = roots[2], xr = x_1_r[rprime]) ) ) {
              return(min(temp_dataset))
            } else {
              return(max(temp_dataset))
            }
          } else {
            if (!silent) message(e)
            return(e)
          }
        }
      )

      exponent <- roots[1]

      if (KS == "Below") {
        xmax <- xother
        temp_dataset <- temp_dataset[temp_dataset <= xmax]
      } else if (KS == "Above") {
        xmin <- xother
        temp_dataset <- temp_dataset[temp_dataset >= xmin]
      }
      if (xmin + 1 > xmax) next

    } else if (otherside == "Basic") {
      temp_dataset <- dataset[xmin <= dataset & dataset <= xmax]
      lower <- 0
      upper <- 100
      arun <- uniroot_func(lower, d = temp_dataset, xmin = xmin, xmax = xmax)
      while( (is.infinite(arun) || is.na(arun) || is.nan(arun)) && lower < 10 ) {
        lower = lower + 1
        arun <- uniroot_func(lower, d = temp_dataset, xmin = xmin, xmax = xmax)
      }
      if (lower >= 10) next

      arun <- uniroot_func(upper, d = temp_dataset, xmin = xmin, xmax = xmax)
      while( (is.infinite(arun) || is.na(arun) || is.nan(arun)) && upper > 10 ) {
        upper = upper - upper/10
        arun <- uniroot_func(upper, d = temp_dataset, xmin = xmin, xmax = xmax)
      }
      if (upper <= lower) next

      exponent <- tryCatch(
        uniroot(uniroot_func, interval = c(lower, upper),
                d = temp_dataset, xmin = xmin, xmax = xmax)$root,
        error = function(e) {
          if (!silent) message(e)
          if (!silent) print(" ")
          if (attributes(e)$class[1] == "simpleError") {
            #Both values likely on same side.
            #Return the interval value that
            #yields a result closer to zero.
            if ( abs( uniroot_func(lower, temp_dataset, xmin, xmax) ) <
                 abs( uniroot_func(upper, temp_dataset, xmin, xmax) ) ) {
              return(lower)
            } else {
              return(upper)
            }
          } else {
            return(e)
          }
        }
      )


    }

    #Compare KS distances######################################################
    fitted <- (xmin : xmax) ^ (- exponent)
    fitted <- cumsum( fitted )
    if (fitted[length(fitted)] == 0) next
    fitted <- fitted / fitted[ length(fitted) ]

    empirical <- tabulate(temp_dataset)
    empirical <- cumsum( empirical )
    empirical <- empirical / empirical[ length(empirical) ]

    while(length(empirical) != length(fitted)) {
      if (length(empirical) < length(fitted)) {
        empirical <- c(empirical, rep(1, length(fitted) - length(empirical)))
      } else if (length(empirical) > length(fitted)) {
        empirical <- empirical[empirical > 0]
      }
    }

    KSdist <- max(abs(
      fitted - empirical
    ))

    if (KSdist < Best_KSd) {
      Best_min <- xmin
      Best_max <- xmax
      Best_exp <- exponent
      Best_KSd <- KSdist
      if (!silent) print(c(Best_min, Best_max, Best_exp, Best_KSd))
    }
  }
  return(list(xmin = Best_min, xmax = Best_max, exponent = Best_exp))
}

#Application: Step-wise Trace:#################################################
#Idea is to run the various methods on the various data sets
#in order to see how the extracted power law exponent behaves
#over time. We also seek to determine if our steady state
#may in truth be in some sense multimodal.

#Output plots: Data set by data set images with lines showing
#the xmin, xmax, and exponent values of each method.
#Additional outputs: Gelman - Rubin statistics.

###Create list of data locations:##############################################
print("Reading Locations")
#Blacklist the places we store our analyses (and their .mat conversions.)
folder_blacklist <- c('Analyses', 'Analyses_no_1s', data_folder)

#Add the data folder as a way to save progress.
if (!dir.exists(file.path(directory, data_folder))) {
  dir.create(file.path(directory, data_folder))
}

if (!dir.exists(file.path(directory, image_folder))) {
  dir.create(file.path(directory, image_folder))
}

folders <- dir(directory, full.names = TRUE)

#Remove inaccessible.
folders <- folders[dir.exists(folders)]

#Access Data Folder.
folders <- folders[grepl('Data', folders, fixed = TRUE)]

#Access subfolders.
folders <- dir(folders, full.names = TRUE)

#Remove inaccessible.
folders <- folders[dir.exists(folders)]

###Create list of methods, names, and arguments:###############################
print("Preparing Methods")
#TODO modify for additional methods
methods <- list(
  names = c(
    "CSN Classic",
    "CSN Bootstrap",
    "Basic Truncated MLE",
    "MaxEnt, All Data",
    "MaxEnt, Boundaries as Variables",
    "CSN + ZXX MLE and KS: Lower Bound",
    "CSN + ZXX's Order Statistics, KS: Lower Bound",
    "CSN + ZXX's Order Statistics, KS: Upper Bound"
  ),
  functions = c(
    fit_poweRlaw_dis,
    PowerLawFitSimple,
    fit_displ_trunc_max_like,
    trunc_maximize_entropy,
    trunc_maximize_entropy,
    fit_trunc_CSN_ZXX,
    fit_trunc_CSN_ZXX,
    fit_trunc_CSN_ZXX
  ),
  arguments = list(
    list(
      calculate_p = 0
    ),
    list(
      NumberOfSamples = 1000,
      NumberOfBootstraps = 1000,
      CensorTail = 0,
      Print = 0,
      NumberOfThreads = 1
    ),
    list(
      NULL
    ),
    list(
      fixed_edges = 1,
      bsval = 0, #<- Evidently does not work.
      steptolval = 50,
      itermaxval = 1000
    ),
    list(
      fixed_edges = 0,
      bsval = 0, #<- Evidently does not work.
      steptolval = 50,
      itermaxval = 1000
    ),
    list(
      KS = "Below",
      otherside = "Basic",
      silent = 1
    ),
    list(
      KS = "Below",
      otherside = "Order",
      silent = 1
    ),
    list(
      KS = "Above",
      otherside = "Order",
      silent = 1
    )
  ),
  bounded = c(#Specifically above. Below is assumed.
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1
  )
)

#Identify the used entries.
methods_flux <- lapply(methods, function(x) {
  lapply(seq_along(x), function(y) {
    if (any(methods_to_use_flux == y))
      x[[y]]
  })
})

#Identify the used entries.
methods <- lapply(methods, function(x) {
  lapply(seq_along(x), function(y) {
    if (any(methods_to_use == y))
      x[[y]]
  })
})

#Trim Unused Entries.
methods_flux <- lapply(methods_flux, function(x) {
  x[lapply(x, length) > 0]
})
methods <- lapply(methods, function(x) {
  x[lapply(x, length) > 0]
})

###Create storage for results.#################################################
print("Creating storage")
kern_names <- basename(folders)

file_names_patt <- if (sparse_mode) '.*_S[.]mat$' else '.*[^(_S)][.]mat$'
file_names_full <- unlist(sapply(folders, function(x) {
  dir(x, pattern = file_names_patt, full.names = TRUE)
  }))
file_names_base <- basename(file_names_full)

number_of_entries <- (length(methods$names) * estimate_rows_per_file)
result_storage_subset <- data.frame(
  time = rep(0, number_of_entries),
  xmin = rep(0, number_of_entries),
  xmax = rep(Inf, number_of_entries),
  clusters = rep(0, number_of_entries),
  exponent = rep(Inf, number_of_entries),
  method_type = factor(rep("Empty", number_of_entries),
                       levels = c("Empty", methods$names)),
  error = rep(NA, number_of_entries),
  kernel = factor(rep("Empty", number_of_entries),
                  levels = c("Empty", unique(unlist(kern_names)))),
  file = factor(rep("Empty", number_of_entries),
                levels = c("Empty", unique(unlist(file_names_base))))
)

flux_storage_subset <- data.frame( #TODO think about this and how to pack pars.
  time = rep(0, number_of_entries),
  xmin = rep(0, number_of_entries),
  xmax = rep(Inf, number_of_entries),
  clusters = rep(0, number_of_entries),
  exponent = rep(Inf, number_of_entries),
  pval = rep(Inf, number_of_entries),
  method_type = factor(rep("Empty", number_of_entries),
                       levels = c("Empty", methods$names)),
  error = rep(NA, number_of_entries),
  kernel = factor(rep("Empty", number_of_entries),
                  levels = c("Empty", unique(unlist(kern_names)))),
  file = factor(rep("Empty", number_of_entries),
                levels = c("Empty", unique(unlist(file_names_base))))
)

#Counter for accessing storage correctly.
rs_entry <- 1

###Begin Procedure:############################################################
print("Beginning Procedure")
# Choose File
file_this <- file_names_full[cargs[1]]
print(file_this)
# Open File
data <- if (substr(file_this,
                   nchar(file_this) - 3,
                   nchar(file_this)) == '.mat') {
  R.matlab::readMat(file_this)
} else NULL
print("Is null?")
print(is.null(data))
stopifnot(!is.null(data))

# Discard unused portion of File
# fluxC <- data[['stitchCFx']]
# fluxF <- data[['stitchFFx']]
# data <- data[['stitchPop']]
data <- data[["stitch"]]
data <- data[min(nrow(data), cargs[2]):min(nrow(data), cargs[3]), ]

directory_this <- file.path(
  directory, data_folder, dirname(file_this)
)

if (!dir.exists(directory_this)) {
  dir.create(directory_this, showWarnings = FALSE, recursive = TRUE)
}

file_save_name <- file.path(
  directory_this,
  paste0(substr(basename(file_this),
                1, nchar(basename(file_this)) - 4),
         '_', tolower(data_folder),
         '_', cargs[2], '_', cargs[3], '.RData')
)
print(file_save_name)

if(!overwrite && file.exists(file_save_name)) {
    print("File already exists. Skipping.")
} else {
	# Fit used portion, line by line
	for (data_row in 1:nrow(data)) {
	  for (f in 1 : length(methods_to_use)) {
		data_temp <- data[data_row, ]
		data_temp <- rep(which(data_temp > 0), data_temp[which(data_temp > 0)])

		fit_args <- formals(methods$functions[[f]])
		fit_args[[1]] <- data_temp
		fit_args <- modifyList(fit_args, methods$arguments[[f]])

		tryCatch({
		  retval <- tryCatch({
			do.call(what = methods$functions[[f]],
					args = as.list(fit_args))
		  },
		  error = function(e) {
			return(list(
			  xmin = 0,
			  xmax = Inf,
              clusters = 0,
			  exponent = Inf,
			  errorval = e
			))}
		  )
		},
		error = function(e) {
		  print(paste('Function Eval. Error:', cargs[1]))
		  print(e)
		})

		retval <-
		  if (is.null(retval$errorval)) {
			data.frame(
			  time = data_row,
			  xmin = retval$xmin,
			  xmax = retval$xmax,
              clusters = sum(data[data_row, ]),
			  exponent = retval$exponent,
			  method_type = methods$names[[f]],
			  error = FALSE
			)
		} else {
			data.frame(
			  time = data_row,
			  xmin = retval$xmin,
			  xmax = retval$xmax,
              clusters = sum(data[data_row, ]),
			  exponent = retval$exponent,
			  method_type = methods$names[[f]],
			  error = paste(format(retval$call), ':', retval$message)
			)
		}

		result_storage_subset[rs_entry, ] <- retval
		rs_entry <- rs_entry + 1
	  }
	}

  if(!debug_skip_flux) {
  #TODO
  # Perform a similar idea, but for the coalescence and fragmentation fluxes.
  # Need to perform the idea for four methods:
  # poweRlaw::displ,
  # poweRlaw::dislnorm,
  # CSN + ZXX truncated powerlaw
  # CSN + ZXX truncated powerlaw + KS
  # We record all fitted parameters, including probability of acceptance
  # where possible.
  #

  for (data_row in 1:nrow(fluxC)) {
    for (s in smoothing) {
      if(data_row + s > nrow(fluxC))
        next
      for (f in 1 : length(methods_to_use_flux)) {
        data_temp <- fluxC[data_row:data_row+s, ]
        # Same code for population, but accepts multiple rows at once.
        # Map values: order is preserved after modulo columns, then map back.
        data_temp <- ( rep(which(data_temp > 0),
                           data_temp[which(data_temp > 0)]) - 1 ) %% ncol(fluxC) + 1

        fit_args <- formals(methods_flux$functions[[f]])
        fit_args[[1]] <- data_temp
        fit_args <- modifyList(fit_args, methods_flux$arguments[[f]])

        #TODO modify returns.
        tryCatch({
          retval <- tryCatch({
            do.call(what = methods_flux$functions[[f]],
                    args = as.list(fit_args))
          },
          error = function(e) {
            return(list(
              xmin = 0,
              xmax = Inf,
              clusters = 0,
              exponent = Inf,
              errorval = e
            ))}
          )
        },
        error = function(e) {
          print(paste('Function Eval. Error:', cargs[1]))
          print(e)
        })

        retval <-
          if (is.null(retval$errorval)) {
            data.frame(
              time = data_row,
              xmin = retval$xmin,
              xmax = retval$xmax,
              clusters = sum(data[data_row, ]),
              exponent = retval$exponent,
              method_type = methods$names[[f]],
              error = FALSE
            )
          } else {
            data.frame(
              time = data_row,
              xmin = retval$xmin,
              xmax = retval$xmax,
              clusters = sum(data[data_row, ]),
              exponent = retval$exponent,
              method_type = methods$names[[f]],
              error = paste(format(retval$call), ':', retval$message)
            )
          }

        flux_storage_subset[rs_entry, ] <- retval
        rs_entry <- rs_entry + 1
      }
    }
  }
  } else {
    flux_storage_subset <- NULL
  }

	# Save the collection of fits to a new file.
	save(result_storage_subset, flux_storage_subset,
	     file = file_save_name, compress = TRUE)
}
}
