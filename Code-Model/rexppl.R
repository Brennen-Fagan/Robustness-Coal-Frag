runtrial <- 1
if (runtrial) {
    # > runif(1)*10^8
    # [1] 16764847
    set.seed(16764847)
}

exppl <- function(i, m, v) {
    retval <- ifelse(i > m | i < 1 | i %% 1 != 0, 0, i)

    index_nonzero <- which(retval != 0)

    n1 <- 1/(2 - v)

    rvnz <- retval[index_nonzero]

    retval[index_nonzero] <- ifelse(
        rvnz == 1, n1,
        exp(
            (rvnz - 1) * log(1 - v) - (2 * rvnz - 1) * log(2 - v) +
                lfactorial(2 * rvnz - 2) - 2 * lfactorial(rvnz)
        )
        )

    return(retval)
}

exppl_norm <- function(i, m, v) {
    temp <- exppl(1:m, m, v)
    return(temp[i] / sum(temp))
}

exppl_normcs <- function(i, m, v) {
    temp <- cumsum(exppl(1:m, m, v))
    return(temp[i] / tail(temp, n = 1))
}

exppl_normcs_comp <- function(i, m, v) {
    temp <- rev(cumsum(as.numeric(rev(exppl(1:m, m, v)))))
    return(temp[i] / head(temp, n = 1))
}

rexppl <- function(n, m, v, complementary = TRUE) {
    if (!complementary) {
        dist <- exppl_normcs(1:m, m, v)
        cut(runif(n), c(0, dist), include.lowest = TRUE, labels = FALSE)
    }
    else {
        dist <- exppl_normcs_comp(1:m, m, v)
        m - cut(runif(n), c(dist, 0), include.lowest = TRUE, labels = FALSE) + 1
    }
}

rexppl_partitions <- function(n, m, v) {
    # Approximately rexppl partitions, rather than exact.
    vals <- rexppl(n * m, m, v)
    parts <- list()
    for (i in 1:n) {
        valindex <- (which.max(cumsum(vals) > m) - 1)
        parts[[i]] <- vals[1:valindex]
        vals <- vals[(valindex + 1):length(vals)]
        if (sum(parts[[i]]) != m)
            parts[[i]][length(parts[[i]]) + 1] <- m - sum(parts[[i]])
    }
    return(parts)
}

if (runtrial) {
    methods_to_use <- c(
        1, # Truncated MLE
        2#, # KS to choose lower bound, MLE
        # 3,
        # 4
    )
    number_of_entries <- 1000

    #Methods:######################################################################
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

            while (length(empirical) != length(fitted)) {
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

    # A different method....:###################################################

    # if ('maxLik' %in% library()$results[, 1]) {
    #     expcutoffpl_LogLikelihood <- function(params, dataset) {
    #         # I.e. p(x) propto exp(gamma*x) * x^alpha
    #         mn <- min(dataset, na.rm = TRUE)
    #         mx <- max(dataset, na.rm = TRUE)
    #         n <- length(dataset)
    #         const <- sum(exp(params['gamma'] * mn:mx) * (mn:mx)^params['alpha'])
    #         return(- n * log(const)
    #                + params['gamma'] * sum(dataset)
    #                + params['alpha'] * sum(log(dataset)))
    #     }
    #
    #     expcutoffpl <- function(x, gamma, alpha, xmin, xmax, log.p = FALSE) {
    #         retval <- ifelse(x < xmin | x > xmax | x %% 1 != 0, 0, x)
    #         rvnz <- retval[retval != 0]
    #         inverseconstant <-  sum(exp(gamma * xmin:xmax) * (xmin:xmax)^alpha)
    #         logvals <- gamma * rvnz + alpha * log(rvnz) - log(inverseconstant)
    #         if (log.p) {
    #             retval[retval != 0] <- logvals
    #         } else {
    #             retval[retval != 0] <- exp(logvals)
    #         }
    #         return(retval)
    #     }
    #
    #     expcutoffpl_normcs <- function(x, gamma, alpha, xmin, xmax) {
    #         temp <- (cumsum(as.numeric((
    #             expcutoffpl(xmin:xmax, gamma, alpha, xmin, xmax)
    #         ))))
    #         return(temp[x] / tail(temp, n = 1))
    #     }
    #
    #     expcutoffpl_normcs_comp <- function(x, gamma, alpha, xmin, xmax) {
    #         temp <- rev(cumsum(as.numeric(rev(
    #             expcutoffpl(xmin:xmax, gamma, alpha, xmin, xmax)
    #         ))))
    #         return(temp[x] / head(temp, n = 1))
    #     }
    #
    #     rexpcutoffpl <- function(n, gamma, alpha, xmin, xmax, complementary = TRUE) {
    #         if (!complementary) {
    #             dist <- expcutoffpl_normcs(xmin:xmax, gamma, alpha, xmin, xmax)
    #             cut(runif(n), c(0, dist), include.lowest = TRUE, labels = FALSE)
    #         }
    #         else {
    #             dist <- expcutoffpl_normcs_comp(xmin:xmax, gamma, alpha, xmin, xmax)
    #             xmax - cut(runif(n), c(dist, 0), include.lowest = TRUE, labels = FALSE) + 1
    #         }
    #     }
    #
    #     fit_expcutoffpl <- function(dataset) {
    #         tempwrapper <- function(params) {
    #             expcutoffpl_LogLikelihood(params, dataset)
    #         }
    #         maxLik::maxLik(tempwrapper,
    #                        start = c('gamma' = -10,
    #                                  'alpha' = -2.5),
    #                        constraints = list(
    #                            ineqA = matrix(c(-1,0, 0, -1),2,2),
    #                            ineqB = c(0,0)
    #                            )
    #                        ) # ineqA %*% params + ineqB > 0
    #     }
    #
    #     fit_expcutoffpl_mle <- function(dataset) {
    #         retval <- fit_expcutoffpl(dataset)
    #         return(list(
    #             xmin = min(dataset),
    #             xmax = max(dataset),
    #             exponent = retval$estimate[2],
    #             errorval = ifelse(retval$code == 0, NA, retval$code),
    #             adtl1 = retval$estimate[1]
    #         ))
    #     }
    #
    #     fit_expcutoffpl_ksmle <- function(dataset) {
    #         # For each candidate xmin, calculate params
    #         # Realistically only need to do two at a time: best and current.
    #         candidates <- sort(unique(dataset))
    #         candidates <- candidates[1:(length(candidates) - 1)]
    #         xmax <- max(dataset)
    #
    #         best <- list(
    #             ksdist = NA,
    #             xmin = NA,
    #             xmax = xmax,
    #             adtl1 = NA,
    #             exponent = NA
    #         )
    #
    #
    #         for (xmin in candidates) {
    #             tempdata <- dataset[dataset >= xmin]
    #             if (length(tempdata) < 10) {
    #                 break
    #             }
    #             retval <- fit_expcutoffpl(tempdata)
    #             # Calculate ksdist
    #             thispmf <- expcutoffpl(xmin:xmax,
    #                                     retval$estimate[1],
    #                                     retval$estimate[2],
    #                                     xmin, xmax)
    #             thiscdf <- cumsum(thispmf)
    #
    #             empipmf <- tabulate(tempdata)
    #             empicdf <- cumsum(empipmf)
    #             empicdf <- empicdf/empicdf[length(empicdf)]
    #
    #             while (length(empicdf) != length(thiscdf)) {
    #                 if (length(empicdf) < length(thiscdf)) {
    #                     empicdf <- c(empicdf, rep(1, length(thiscdf) - length(empicdf)))
    #                 } else if (length(empicdf) > length(thiscdf)) {
    #                     empicdf <- empicdf[empicdf > 0]
    #                 }
    #             }
    #
    #             thisksdist <- max(abs(
    #                 thiscdf - empicdf
    #             ))
    #
    #             # end
    #             if (is.na(best$ksdist) || best$ksdist > thisksdist) {
    #                 best$ksdist <- thisksdist
    #                 best$xmin <- xmin
    #                 best$xmax <- xmax
    #                 best$adtl1 <- retval$estimate[1]
    #                 best$exponent <- retval$estimate[2]
    #             }
    #         }
    #
    #         return(best)
    #
    #     }
    # }

    ############################################################################

    methods <- list(
        names = c(
            "Basic Truncated MLE",
            "CSN + ZXX MLE and KS: Lower Bound"#,
            # "Exponential Cutoff PL MLE",
            # "Exponential Cutoff PL KSMLE"
        ),
        functions = c(
            fit_displ_trunc_max_like,
            fit_trunc_CSN_ZXX#,
            # fit_expcutoffpl_mle,
            # fit_expcutoffpl_ksmle
        ),
        arguments = list(
            list(
                NULL
            ),
            list(
                KS = "Below",
                otherside = "Basic",
                silent = 1
            )#,
            # list(
            #     NULL
            # ),
            # list(
            #     NULL
            # )
        ),
        bounded = c(#Specifically above. Below is assumed.
            1,
            1#,
            # 1,
            # 1
        )
    )

    #Identify the used entries.
    methods <- lapply(methods, function(x) {
        lapply(seq_along(x), function(y) {
            if (any(methods_to_use == y))
                x[[y]]
        })
    })

    #Trim Unused Entries.
    methods <- lapply(methods, function(x) {
        x[lapply(x, length) > 0]
    })

    result_storage_subset <- data.frame(
        time = rep(0, number_of_entries),
        xmin = rep(0, number_of_entries),
        xmax = rep(Inf, number_of_entries),
        exponent = rep(Inf, number_of_entries),
        adtl1 = rep(NA, number_of_entries),
        method_type = factor(rep("Empty", number_of_entries),
                             levels = c("Empty", methods$names)),
        N = rep(0, number_of_entries),
        n_1 = rep(0, number_of_entries),
        error = rep(NA, number_of_entries),
        stringsAsFactors = FALSE
    )


    #Counter for accessing storage correctly.
    rs_entry <- 1

    data <- rexppl_partitions(number_of_entries, 10000, 0.2)

    for (data_row in 1:length(data)) {
        print(data_row)
        for (f in 1 : length(methods_to_use)) {
            print(f)
            data_temp <- data[[data_row]]

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
                        exponent = Inf,
                        adtl1 = NA,
                        errorval = e
                    ))}
                )
            },
            error = function(e) {
                print(paste('Function Eval. Error:', data_row))
                print(e)
            })

            retval <-
                if (is.null(retval$errorval) || is.na(retval$errorval)) {
                    data.frame(
                        time = data_row,
                        xmin = retval$xmin,
                        xmax = retval$xmax,
                        exponent = retval$exponent,
                        adtl1 = ifelse(is.null(retval$adtl1), NA, retval$adtl1),
                        method_type = methods$names[[f]],
                        N = sum(data[[data_row]]),
                        n_1 = sum(data[[data_row]] == 1),
                        error = FALSE,
                        stringsAsFactors = FALSE
                    )
                } else {
                    data.frame(
                        time = data_row,
                        xmin = retval$xmin,
                        xmax = retval$xmax,
                        exponent = retval$exponent,
                        adtl1 = retval$adtl1,
                        method_type = methods$names[[f]],
                        N = sum(data[[data_row]]),
                        n_1 = sum(data[[data_row]] == 1),
                        error = paste(format(retval$call), ':', retval$message),
                        stringsAsFactors = FALSE
                    )
                }

            result_storage_subset[rs_entry, ] <- retval
            rs_entry <- rs_entry + 1
        }
    }

    library(ggplot2)
    library(dplyr)

    Method_conversion <- c(
        "Basic Truncated MLE" = 'MLE',
        "CSN + ZXX MLE and KS: Lower Bound" = 'KS-MLE'#,
        # "Exponential Cutoff PL MLE" = 'EPL MLE',
        # "Exponential Cutoff PL KSMLE" = 'EPL KS-MLE'
    )

    chunking <- 500
    ggplotTheme <- function(
        ..., base_size = 22
    ) {
        ggplot2::theme_bw(base_size) + ggplot2::theme(...)
    }

    result_storage_subset <-
        result_storage_subset %>% dplyr::mutate(
            method_type = Method_conversion[as.character(method_type)]
        ) %>% dplyr::mutate(
            M = 10000,
            time = time * 10, # Sample Time to Event Time.
            Simulation = "Dummy",
            kernel = "Dummy",
            file = "Dummy",
            n1 = n_1
        ) %>% dplyr::select(
            -n_1
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
                rss_xmin, by = c('time', 'xmax', 'M', 'n1', 'N', 'Simulation',
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
                    rss_xmin, by = c('time', 'xmax', 'M', 'n1', 'N', 'Simulation',
                                     'error', 'kernel', 'file'),
                    suffix = c('exponent', 'xmin')
                ) %>% dplyr::arrange(time)
        }
    }

    result_storage_subset <- result_storage_subset_out

    # # Detect appropriate file from files_data to load.
    # fname_run <- strsplit(fname, split = "_")
    # print(fname_run)
    # if (fname_run[[1]][7] == "Run") {
    #     fname_run = fname_run[[1]][8]
    # } else {
    #     stop("Error in file loaded. Does not conform to expected convention.")
    # }
    #
    # file_data <- grep(pattern = paste0("10000_Run_", fname_run),
    #                   x = files_data,
    #                   value = TRUE, fixed = TRUE)
    #
    # file_data_matrix <- R.matlab::readMat(
    #     file.path(file_data)
    # )$stitch # Name of the Matrix in MATLAB format.
    #
    # # Add relevant values to result_storage_subset.
    # result_storage_subset$N <- rowSums(file_data_matrix)
    # result_storage_subset$n1 <- file_data_matrix[, 1]

    # Save.
    res_phy <- result_storage_subset

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
        ggplot2::aes(x = `KS-MLE:xmin`,
                     y = `KS-MLE:exponent`,
                     color = `KS-MLE:xmin`)
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
    ) + ggplot2::scale_x_discrete(name = "Method") + ggplot2::coord_cartesian(
        ylim = c(2.45, 3)
    ) # 900 x 600 method_violinsSimulated
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

}



