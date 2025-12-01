

#' Summarize data frame
#'
#' Helper function to summarize data for input to [SLDAssay::get.mle]
#'
#' @param counts Column name of counts data
#' @param dilution Column name of dilutions tested
#' @param cutoff Threshold to score positive wells
#' @param groupby Optional sample identifier
#' @returns A data.frame object summarizing results of LDA
#' @import dplyr
#' @importFrom rlang sym
format_for_slda <- function(df, counts, dilutions, cutoff, groupby = NULL){
  # Group by multiple samples if provided
  if (!is.null(groupby)){
    df <- df %>% group_by( !!sym(groupby), .add = TRUE)
  }
  df <- df %>%
    group_by( !!sym(dilutions) , .add = TRUE) %>%
    mutate(pos_wells = !!sym(counts) > cutoff) %>%
    mutate(pos = sum(pos_wells == TRUE),
           neg = sum(pos_wells == FALSE),
           replicates = n()) %>%
    distinct(pos,neg,replicates) %>%
    ungroup() %>%
    mutate(pct_neg = 100*neg/replicates)
  return (df)
}

#' Convert list to data frame
#'
#' Helper function to format list output from [SLDAssay::get.mle] into a data frame
#'
#' @param results List returned from [SLDAssay::get.mle]
#' @returns A data.frame of results
format_results_as_df <- function(results){
  # Split multiple values into separate entries (min,max)
  for (name in names(results)){
    value <- unlist(results[name])
    if (length(value)==2){
      results[paste0(name,"_min")] = value[[1]]
      results[paste0(name,"_max")] = value[[2]]
      results[name] <- NULL
    }
  }
  return(as.data.frame(results))
}

#' Get Maximum Likelihood Estimate (MLE)
#'
#' Calculates the MLE for LDA experiment(s). Uses [SLDAssay::get.mle] under the hood.
#'
#' Returned statistics include:
#'    - `MLE`: MLE
#'    - `BC_MLE`: Bias-corrected MLE
#'    - `Exact_PGOF`: P-value for goodness of fit. The probability of an experimental result as rare as or rarer than that obtained, assuming that the model is correct. Low values (e.g. < 0.01) indicate rare or implausible experimental results. Samples with a very low values might be considered for retesting
#'    - `Asymp_PGOF`: P-value calculated using an asymptotic Chi-Squared distribution with D-1 degrees of freedom, where D is the number of dilution levels
#'    - `Exact_CI`: Exact confidence interval, computed from the likelihood ratio test (recommended)
#'    - `Asymp_CI`: Wald asymptotic confidence interval, based on the normal approximation to the binomial distribution
#'
#' @param df A data.frame object containing columns 'counts' and 'dilutions' or a path to a .csv file containing those columns
#' @param counts Character name of column in df where cell counts are located
#' @param dilutions Character name of column in df where cell dilutions are located
#' @param groupby Optional character name of column in df containing sample or group names. The MLE will be calculated for each sample in this column. If not provided, will treat the entire data frame as one sample (the default)
#' @param cutoff Threshold for scoring counts as positive or negative (default: 25 cells)
#' @param conf_level The confidence level of the interval (default: 0.95)
#' @param mc Number of Monte Carlo samples. Default is exact (no MC sampling), unless more than 15,000 possible positive well outcomes exist, in which case 15,000 MC samples are taken. Use mc=FALSE for exact computation
#' @returns Results as a list, where `results$summary` is the data frame summarizing the LDA experiment, and `results$statistics` is the maximum likelihood estimate, confidence interval, and goodness of fit p-value
#' @importFrom SLDAssay get.mle
#' @importFrom utils read.csv
#' @importFrom rlang sym
#' @importFrom dplyr bind_rows select everything
#' @export
getMLE <- function(df,
                    counts,
                    dilutions,
                    groupby = NULL,
                    cutoff = 25,
                    conf_level=0.95,
                    mc=15000){
  # Load df from csv or validate is data.frame
  if (inherits(df,"character")){
    if (file.exists(df)){
      df <- read.csv(df)
    }else{
      stop(sprintf("File does not exist: %s", df))
    }
  }else if (!inherits(df,"data.frame")){
    stop("Must provide a data.frame for df")
  }
  # Check columns exist
  if (!(counts %in% colnames(df)) | !(dilutions %in% colnames(df))){
    stop(sprintf("%s and %s must be column names in df", counts, dilutions))
  }else if (!is.null(groupby)){
    if (!(groupby %in% colnames(df))) stop(sprintf("%s must be a column name in df", groupby))
  }
  # Format for SLDAssay
  df <- format_for_slda(df,
                        counts = counts,
                        dilutions = dilutions,
                        cutoff = cutoff,
                        groupby = groupby)

  # Get MLE for one group
  if (is.null(groupby)){
    message("Calculating MLE...")
    res <- get.mle(df[["pos"]],
                   df[["replicates"]],
                   df[[dilutions]],
                   conf.level = conf_level,
                   monte = mc,
                   iupm = FALSE,
                   na.rm = FALSE)
    res <- format_results_as_df(res)
  # Get MLE per group
  }else{
    res <- lapply(unique(df[[groupby]]), function(group){
      message(sprintf("Calculating MLE for %s...", group))
      # Select group for MLE
      df_sub <- df[df[[groupby]] == group,]
      res_sub <- get.mle(df_sub[["pos"]],
                         df_sub[["replicates"]],
                         df_sub[[dilutions]],
                         conf.level = conf_level,
                         monte = mc,
                         iupm = FALSE,
                         na.rm = FALSE)
      # Add group name
      res_sub[[groupby]] <- group
      # Format as data frame and return
      format_results_as_df(res_sub)
    })
    # Make into single data frame
    res <- bind_rows(res)
    # Make groupby column the first
    res <- select(res, !!sym(groupby), everything())
  }

  # Add to results data frame
  res["CONF_LEVEL"] <- conf_level
  res["CUTOFF"] <- cutoff
  res["MONTE_CARLO"] <- mc
  message("Done")
  return(list(summary = df, statistics = res))
}
