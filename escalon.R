require(e1071)
require(doParallel)
require(signal)
library(ggplot2)
library(gridExtra)
require(RColorBrewer)
require(pracma)
require(scales)


##########################Funciones######################

############# SVM
#retardos
retardos_multi <- function( signalData, lags){
  signal.uni <- signalData
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag)
  col.names <- list("PAMn","VFSCn")
  columns <- NULL
  lagged.columns.names <- c()
  for(colname in col.names){
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname]
    if(!is.null(lag.order) && lag.order > 0)
      for(i in 1:lag.order){
        new.colname <- paste(colname, paste0("lag", i), sep = ".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- signal.uni[lag.mat[, i+1], colname]
      }
  }
  folded.signal <- data.frame(columns)
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}

normalise.signal <- function(
  signal,
  signal.baseline.value,
  signal.min.value
)
{
  (signal - signal.min.value) / abs(signal.baseline.value - signal.min.value)
}

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


############ mfARI 
are.tolerably.equal <- function(x, y, tol, ...) UseMethod("tolerably.equal")

are.tolerably.equal <- function(x, y, tol = .Machine$double.eps^0.5, ...)
{
  abs(x - y) < tol
}

estimate.mfARI.CBFV.parameters.type.2 <- function(
  time.instants,
  normalised.CBFV.signal,
  min.CBFV.sample,
  min.CBFV.time.instant,
  transient.CBFV.duration,
  steady.CBFV.duration,
  time.tol = min(diff(time.instants)) / 100,
  ...
)
{
  ans <- list()
  
  ans[["Delta.tau"]] <- transient.CBFV.duration
  ans[["tau"]] <- min.CBFV.time.instant + transient.CBFV.duration
  ans[["transient.CBFV.samples"]] <- which(
    time.instants >= min.CBFV.time.instant - time.tol &
      time.instants <= ans[["tau"]] + time.tol
  )
  ans[["transient.CBFV.time.instants"]] <- 
    time.instants[ans[["transient.CBFV.samples"]]]
  ans[["transient.normalised.CBFV.signal"]] <-
    normalised.CBFV.signal[ans[["transient.CBFV.samples"]]]
  ans[["transient.CBFV.slope"]] <- 
    mean(ans[["transient.normalised.CBFV.signal"]]) /
    mean(ans[["transient.CBFV.time.instants"]] - min.CBFV.time.instant)
  ans[["transient.CBFV.line"]] <- 
    ans[["transient.CBFV.slope"]] *
    (time.instants[ans[["transient.CBFV.samples"]]] -
       min.CBFV.time.instant)
  if(length(ans[["transient.normalised.CBFV.signal"]]) > 0)
    ans[["transient.CBFV.MSE"]] <- get.MSE(
      sim = ans[["transient.CBFV.line"]],
      obs = ans[["transient.normalised.CBFV.signal"]]
    )
  else
    ans[["transient.CBFV.MSE"]] <- Inf
  ans[["Ks"]] <- tail(ans[["transient.normalised.CBFV.signal"]], 1)
  #ans[["Ks"]] <- tail(normalize(ans[["transient.normalised.CBFV.signal"]]), 1)
  
  
  ans[["nominal.steady.CBFV.duration"]] <- steady.CBFV.duration
  ans[["steady.CBFV.samples"]] <- which(
    time.instants >= ans[["tau"]] - time.tol &
      time.instants <= ans[["tau"]] + steady.CBFV.duration + time.tol
  )
  ans[["steady.CBFV.time.instants"]] <-
    time.instants[ans[["steady.CBFV.samples"]]]
  steady.CBFV.length <- length(ans[["steady.CBFV.samples"]])
  ans[["steady.CBFV.duration"]] <-
    ans[["steady.CBFV.time.instants"]][length(ans[["steady.CBFV.time.instants"]])] -
    ans[["tau"]]
  ans[["steady.normalised.CBFV.signal"]] <-
    normalised.CBFV.signal[ans[["steady.CBFV.samples"]]]
  ans[["steady.CBFV.line"]] <- rep(ans[["Ks"]], steady.CBFV.length)
  if(steady.CBFV.length > 0)
    ans[["steady.CBFV.MSE"]] <- get.MSE(
      sim = ans[["steady.CBFV.line"]],
      obs = ans[["steady.normalised.CBFV.signal"]]
    )
  else
    ans[["steady.CBFV.MSE"]] <- Inf
  
  ans[["CBFV.response.duration"]] <-
    ans[["Delta.tau"]] +
    ans[["steady.CBFV.duration"]]
  if(is.finite(ans[["transient.CBFV.MSE"]]) &&
     is.finite(ans[["steady.CBFV.MSE"]]))
    ans[["CBFV.response.MSE"]] <-
    (ans[["Delta.tau"]] / ans[["CBFV.response.duration"]]) *
    ans[["transient.CBFV.MSE"]] +
    (ans[["steady.CBFV.duration"]] / ans[["CBFV.response.duration"]]) *
    ans[["steady.CBFV.MSE"]]
  else
    ans[["CBFV.response.MSE"]] <- Inf
  
  return(ans)
}


search.mfARI.CBFV.parameters <- function(
  time.instants,
  normalised.CBFV.signal,
  min.CBFV.sample,
  min.CBFV.time.instant,
  min.transient.CBFV.duration,
  max.transient.CBFV.duration,
  steady.CBFV.duration,
  time.tol = min(diff(time.instants)) / 100,
  estimation.function = estimate.mfARI.CBFV.parameters.type.2,
  keep.search.results = FALSE,
  add.search.plots = FALSE,
  id = NULL,
  ...
)
{
  min.delta.tau <- min.CBFV.time.instant + min.transient.CBFV.duration
  max.delta.tau <- min.CBFV.time.instant + max.transient.CBFV.duration
  last.time <- max.delta.tau + steady.CBFV.duration
  last.instant <-time.instants[length(time.instants)]
  
  if(last.time > last.instant + time.tol)
  {
    last.time <- last.instant
    max.delta.tau <- last.time - steady.CBFV.duration
    max.transient.CBFV.duration <- max.delta.tau - min.CBFV.time.instant
    if(is.null(id))
      msg <- "Warning: signal is too short for specified parameters;"
    else
      msg <- sprintf(
        "Warning: signal is too short for specified parameters (%s);",
        id
      ) 
    msg <- paste(
      msg,
      "maximum transient CBFV signal duration was reduced to")
    msg <- paste(
      msg,
      sprintf('%.1f seconds', max.transient.CBFV.duration)
    )
    warning(msg)
  }
  i <- which(
    time.instants >= min.delta.tau - time.tol &
      time.instants <= max.delta.tau + time.tol
  )
  
  ans <- list()
  ans[["Delta.tau.values"]] <- time.instants[i] - min.CBFV.time.instant
  
  plots <- list()
  results <- list()
  best.result <- list(CBFV.response.MSE = Inf)
  
  for(delta.tau in ans[["Delta.tau.values"]])
  {
    current.estimation <- estimation.function(
      time.instants = time.instants,
      normalised.CBFV.signal = normalised.CBFV.signal,
      min.CBFV.sample = min.CBFV.sample,
      min.CBFV.time.instant = min.CBFV.time.instant,
      transient.CBFV.duration = delta.tau,
      steady.CBFV.duration = steady.CBFV.duration,
      time.tol = time.tol,
      ...
    )
    
    
    if(current.estimation[["CBFV.response.MSE"]] <
       best.result[["CBFV.response.MSE"]])
      best.result <- current.estimation
    
    
    str.delta.tau <- as.character(delta.tau)
    if(keep.search.results)
    {
      lce <- list(current.estimation)
      names(lce) <- str.delta.tau
      results <- c(results, lce)
    }
    
    if(add.search.plots)
    {
      p <- get.plot.CBFV.parameters.search(
        time.instants,
        normalised.CBFV.signal,
        min.CBFV.sample,
        min.CBFV.time.instant,
        current.estimation,
        id = id,
        ...
      )
      lp <- list(p)
      names(lp) <- str.delta.tau
      plots <- c(plots, lp)
    }
  }
  
  if(keep.search.results)
    ans[["search.results"]] <- results
  
  if(add.search.plots)
    ans[["search.plots"]] <- plots
  
  ans <- c(ans, best.result)
  return(ans)
}


get.MSE <-function(sim, obs, ...) UseMethod("get.MSE")

get.MSE.default <- function(sim, obs, ...)
{
  valid.classes <- c("integer", "numeric", "ts")
  if(!all(inherits(sim, valid.classes), inherits(obs, valid.classes)))
    stop("Invalid argument type: 'sim' & 'obs' have to be of class ",
         valid.classes)
  
  mse <- mean((sim - obs)^2, ...)
  
  return(mse)
}

get.MSE.matrix <- function(sim, obs, ...)
{
  # Check that 'sim' and 'obs' have the same dimensions
  if(!all.equal(dim(sim), dim(obs)))
    stop(paste0("Invalid argument: dim(sim) != dim(obs) ",
                "(", "[", paste(dim(sim), collapse = " "), "]", " != ",
                "[", paste(dim(obs), collapse = " "), "]", ")"))
  
  mse <- colMeans((sim - obs)^2, ...)
  
  return(mse)
}

get.MSE.data.frame <- function(sim, obs, ...)
{
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
  
  get.MSE.matrix(sim = sim, obs = obs, ...)
  
}


normalise.ABP.signal <- function(
  time.instants,
  ABP.signal,
  sample.release,
  sampling.time,
  baseline.initial.time = time.instants[1],
  baseline.final.time = time.instants[sample.release],
  min.ABP.max.delta.time = 5 * 0.8,
  time.tol = sampling.time / 100
)
{
  ans <- list()
  valid <- !is.na(ABP.signal)
  time.release <- time.instants[sample.release]
  
  # Adds baseline data to answer
  if(baseline.final.time < baseline.initial.time)
    stop("final time for the ABP baseline cannot be ",
         "before its initial time")
  if(baseline.final.time > time.release)
    stop("final time for the ABP baseline cannot be ",
         "after the release of the thigh cuffs")
  
  i <- which(valid & time.instants >= baseline.initial.time - time.tol)
  if(length(i) == 0)
    stop("no time instant for the beginning of the ABP baseline could be determined")
  i <- i[1]
  ans[["ABP.baseline.initial.time"]] <- time.instants[i]
  ans[["ABP.baseline.initial.sample"]] <- i
  
  i <- which(valid & time.instants <= baseline.final.time + time.tol)
  if(length(i) == 0)
    stop("no time instant for the end of the ABP baseline could be determined")
  i <- tail(i, 1)
  ans[["ABP.baseline.final.time"]] <- time.instants[i]
  ans[["ABP.baseline.final.sample"]] <- i
  
  ans[["ABP.baseline.samples"]] <-
    ans[["ABP.baseline.initial.sample"]]:ans[["ABP.baseline.final.sample"]]
  ans[["ABP.baseline.duration"]] <-
    ans[["ABP.baseline.final.time"]] - ans[["ABP.baseline.initial.time"]]
  ans[["ABP.baseline.value"]] <- mean(
    ABP.signal[ans[["ABP.baseline.samples"]]],
    na.rm = TRUE
  )
  
  # Sets min ABP window
  ans[["min.ABP.max.delta.time"]] <- min.ABP.max.delta.time
  i <- time.instants > time.release + time.tol
  j <- time.instants < time.release + ans[["min.ABP.max.delta.time"]] +
    time.tol
  ans[["min.ABP.samples"]] <- which(valid & i & j)
  if(length(ans[["min.ABP.samples"]]) == 0)
    stop("no candidates for min ABP")
  
  # Gets minimum ABP
  ans[["min.ABP.window"]] <- ABP.signal[ans[["min.ABP.samples"]]]
  min.sample.in.window <- which.min(ans[["min.ABP.window"]])
  min.ABP.sample <- min.sample.in.window + ans[["min.ABP.samples"]][1] - 1
  ans[["min.ABP.time.instant"]] <- time.instants[min.ABP.sample]
  ans[["min.ABP.sample"]] <- min.ABP.sample
  ans[["min.ABP.value"]] <- ABP.signal[min.ABP.sample]
  
  # Gets min ABP type
  #min.info <- findpeaks(
  #  -ans[["min.ABP.window"]],
  #  zero = "-",
  #  sortstr = TRUE
  #)
  #edit: modificacion de min.info... funcionando
  min.info <- pracma::findpeaks(
    x = -ans[["min.ABP.window"]],
    zero = "-",
    sortstr = TRUE
  )
  ans[["min.ABP.type"]] <- "minimum value in window"
  if(!is.null(min.info) && min.info[1, 2] == min.sample.in.window)
    ans[["min.ABP.type"]] <- "local minimum"
  
  # Measures min ABP drop
  ans[["min.ABP.drop"]] <-
    ans[["ABP.baseline.value"]] - ans[["min.ABP.value"]]
  ans[["min.ABP.drop.pc"]] <-
    ans[["min.ABP.drop"]] / ans[["ABP.baseline.value"]]
  ans[["min.ABP.drop.pc"]] <- round(ans[["min.ABP.drop.pc"]] * 100, 2)
  
  # Normalises the signal
  # ans[["normalised.ABP.signal"]] <- normalise.signal(
  #   signal = ABP.signal,
  #   signal.baseline.value = ans[["ABP.baseline.value"]],
  #   signal.min.value = ans[["min.ABP.value"]]
  # )
  ans[["normalised.ABP.signal"]] <- normalize(ABP.signal)
  invisible(ans)
}


normalise.CBFV.signal <- function(
  time.instants,
  CBFV.signal,
  sample.release,
  sampling.time,
  baseline.initial.time = time.instants[1],
  baseline.final.time = time.instants[sample.release],
  min.CBFV.max.delta.time = 8 * 0.8,
  time.tol = sampling.time / 100
)
{
  ans <- list()
  valid <- !is.na(CBFV.signal)
  time.release <- time.instants[sample.release]
  # Adds baseline data to answer
  if(baseline.final.time < baseline.initial.time)
    stop("final time for the CBFV baseline cannot be before its initial time")
  if(baseline.final.time > time.release)
    stop("final time for the CBFV baseline cannot be after the release of the thigh cuffs")
  
  i <- which(valid & time.instants >= baseline.initial.time - time.tol)
  if(length(i) == 0)
    stop("no time instant for the beginning of the CBFV baseline could be determined")
  i <- i[1]
  ans[["CBFV.baseline.initial.time"]] <- time.instants[i]
  ans[["CBFV.baseline.initial.sample"]] <- i
  
  i <- which(valid & time.instants <= baseline.final.time + time.tol)
  if(length(i) == 0)
    stop("no time instant for the end of the CBFV baseline could be determined")
  i <- tail(i, 1)
  ans[["CBFV.baseline.final.time"]] <- time.instants[i]
  ans[["CBFV.baseline.final.sample"]] <- i
  
  ans[["CBFV.baseline.samples"]] <- ans[["CBFV.baseline.initial.sample"]]:ans[["CBFV.baseline.final.sample"]]
  ans[["CBFV.baseline.duration"]] <- ans[["CBFV.baseline.final.time"]] - ans[["CBFV.baseline.initial.time"]]
  
  # Gets baseline value
  ans[["CBFV.baseline.value"]] <- mean(CBFV.signal[ans[["CBFV.baseline.samples"]]], na.rm = TRUE)
  
  # Sets min CBFV window
  ans[["min.CBFV.max.delta.time"]] <- min.CBFV.max.delta.time
  i <- time.instants > time.release + time.tol
  j <- time.instants < time.release + ans[["min.CBFV.max.delta.time"]] +
    time.tol
  ans[["min.CBFV.samples"]] <- which(valid & i & j)
  if(length(ans[["min.CBFV.samples"]]) == 0)
    stop("no candidates for min CBFV")
  
  # Gets minimum CBFV
  ans[["min.CBFV.window"]] <- CBFV.signal[ans[["min.CBFV.samples"]]]
  min.sample.in.window <- which.min(ans[["min.CBFV.window"]])
  min.CBFV.sample <- min.sample.in.window + ans[["min.CBFV.samples"]][1] - 1
  ans[["min.CBFV.time.instant"]] <- time.instants[min.CBFV.sample]
  ans[["min.CBFV.sample"]] <- min.CBFV.sample
  ans[["min.CBFV.value"]] <- CBFV.signal[min.CBFV.sample]
  
  # Gets min CBFV type
  min.info <- findpeaks(
    -ans[["min.CBFV.window"]],
    zero = "-",
    sortstr = TRUE
  )
  ans[["min.CBFV.type"]] <- "minimum value in window"
  if(!is.null(min.info) && min.info[1, 2] == min.sample.in.window)
    ans[["min.CBFV.type"]] <- "local minimum"
  
  # Normalises the signal
  # ans[["normalised.CBFV.signal"]] <- normalise.signal(
  #   signal = CBFV.signal,
  #   signal.baseline.value = ans[["CBFV.baseline.value"]],
  #   signal.min.value = ans[["min.CBFV.value"]]
  # )
  ans[["normalised.CBFV.signal"]] <- normalize(CBFV.signal)
  
  invisible(ans)
}


get.mfARI.parameters <- function(
  time.instants,
  ABP.signal,
  CBFV.signal,
  sampling.time,
  time.release,
  baseline.initial.time = time.instants[1],
  baseline.final.time = time.release,
  min.ABP.max.delta.time = 5 * 0.8,
  transient.ABP.duration = 6,
  min.CBFV.max.delta.time = 8 * 0.8, 
  min.transient.CBFV.duration = 1.5,
  max.transient.CBFV.duration = 10,
  steady.CBFV.duration = 6,
  min.Ks = 0.0,
  max.Ks = 1.022095,
  min.ABP.angle = 0.0,
  max.ABP.angle = 90.0,
  min.CBFV.angle = 0.0,
  max.CBFV.angle = 90.0,
  min.Phi = 0.0,
  max.Phi = 35.238932,
  keep.details = TRUE,
  ABP.rounding.digits = 2,
  normalised.ABP.rounding.digits = ABP.rounding.digits * 2,
  CBFV.rounding.digits = 2,
  normalised.CBFV.rounding.digits = CBFV.rounding.digits * 2,
  Ks.rounding.digits = 4,
  Phi.rounding.digits = 2,
  time.tol = sampling.time / 100,
  ... # Passed to search.mfARI.CBFV.parameters
)
{
  # Simple validations
  #valores dados en los parametros de la funcion, excepto el sampling time
  if(transient.ABP.duration < sampling.time)
    stop("duration of the transient ABP signal is too short")
  if(steady.CBFV.duration < sampling.time)
    stop("duration of the steady CBFV signal is too short")
  if(max.transient.CBFV.duration < min.transient.CBFV.duration)
    stop("maximum duration of the transient CBFV signal must be equal or higher than the specified minimum duration")
  if(max.Ks < min.Ks)
    stop("maximum value for Ks must be equal or higher than the specified minimum value")
  if(max.ABP.angle < min.ABP.angle)
    stop("maximum value for the ABP recovery angle must be equal or higher than the specified minimum value")
  if(max.CBFV.angle < min.CBFV.angle)
    stop("maximum value for the CBFV recovery angle must be equal or higher than the specified minimum value")
  if(max.Phi < min.Phi)
    stop("maximum value for the difference between ABP and CBFV recovery angles must be equal or higher than the specified minimum value")
  
  # Validates the lengths of the signals
  if(length(ABP.signal) != length(CBFV.signal))
    stop("ABP and CBFV signals must be of the same length")
  if(length(ABP.signal) != length(time.instants))
    stop("number of time instants must coincide with the signals")
  
  # Validates sampling time
  if(sampling.time < 0.1 - time.tol)
    stop("sampling time must be equal or higher than 0.1")
  
  # Defines the detailed answer list
  dans <- list()
  
  # Copies the original signals
  dans[["time.instants"]] <- round(time.instants, 1)
  dans[["ABP.signal"]] <- round(ABP.signal, ABP.rounding.digits)
  dans[["CBFV.signal"]] <- round(CBFV.signal, CBFV.rounding.digits)
  
  # Sets frequency
  dans[["sampling.time"]] <- round(sampling.time, 1)
  dans[["frequency"]] <- 1 / sampling.time
  
  # Sets release sample
  dans[["time.release"]] <- time.release
  i <- which(are.tolerably.equal(time.instants, time.release, time.tol))
  if(length(i) != 1)
    stop("a unique time instant for the thigh-cuff release could not be determined")
  dans[["sample.release"]] <- i
  if(!is.finite(ABP.signal[i]))
    stop("invalid ABP signal value at the time of thigh-cuff release")
  if(!is.finite(CBFV.signal[i]))
    stop("invalid CBFV signal value at the time of thigh-cuff release")
  
  # Normalises ABP 
  nabp <- normalise.ABP.signal(
    time.instants = dans[["time.instants"]],
    ABP.signal = dans[["ABP.signal"]],
    sample.release = dans[["sample.release"]],
    sampling.time = dans[["sampling.time"]],
    baseline.initial.time = baseline.initial.time,
    baseline.final.time = baseline.final.time,
    min.ABP.max.delta.time = min.ABP.max.delta.time,
    time.tol = time.tol
  )
  #cat(indent, "-- normalizando --\n", sep = "")
  nabp[["normalised.ABP.signal"]] <- round(
    nabp[["normalised.ABP.signal"]],
    normalised.ABP.rounding.digits
  )
  dans <- c(dans, nabp)
  
  # Normalises CBFV
  ncbfv <- normalise.CBFV.signal(
    time.instants = dans[["time.instants"]],
    CBFV.signal = dans[["CBFV.signal"]],
    sample.release = dans[["sample.release"]],
    sampling.time = dans[["sampling.time"]],
    baseline.initial.time = baseline.initial.time,
    baseline.final.time = baseline.final.time,
    min.CBFV.max.delta.time = min.CBFV.max.delta.time,
    time.tol = time.tol
  )
  ncbfv[["normalised.CBFV.signal"]] <- round(
    ncbfv[["normalised.CBFV.signal"]],
    normalised.CBFV.rounding.digits
  )

  dans <- c(dans, ncbfv)
  
  if(any(nabp[["ABP.baseline.samples"]] != ncbfv[["CBFV.baseline.samples"]]))
    warning("baseline samples are different for ABP and CBFV")
  
  # Search for best Ks and delta-tau
  search.results <- search.mfARI.CBFV.parameters(
    time.instants = time.instants,
    normalised.CBFV.signal = dans[["normalised.CBFV.signal"]],
    min.CBFV.sample = dans[["min.CBFV.sample"]],
    min.CBFV.time.instant = dans[["min.CBFV.time.instant"]],
    min.transient.CBFV.duration = min.transient.CBFV.duration,
    max.transient.CBFV.duration = max.transient.CBFV.duration,
    steady.CBFV.duration = steady.CBFV.duration,
    time.tol = time.tol,
    ...)
  dans <- c(dans, search.results)
  
  dans[["Delta.tau"]] <- round(dans[["Delta.tau"]], 1)
  
  # Bounds Ks parameter
  if(dans[["Ks"]] < min.Ks)
    dans[["Ks"]] <- min.Ks
  if(dans[["Ks"]] > max.Ks)
    dans[["Ks"]] <- max.Ks
  dans[["Ks"]] <- round(dans[["Ks"]], Ks.rounding.digits)
  
  # Gets transient ABP representation
  dans[["nominal.transient.ABP.duration"]] <- transient.ABP.duration
  nominal.transient.ABP.end.time <- dans[["min.ABP.time.instant"]] + transient.ABP.duration
  dans[["transient.ABP.samples"]] <- which(time.instants >= dans[["min.ABP.time.instant"]] & time.instants <= nominal.transient.ABP.end.time)
  dans[["transient.ABP.time.instants"]] <- time.instants[dans[["transient.ABP.samples"]]]
  dans[["transient.ABP.duration"]] <- dans[["transient.ABP.time.instants"]][length(dans[["transient.ABP.time.instants"]])] - dans[["min.ABP.time.instant"]]
  dans[["transient.normalised.ABP.signal"]] <- dans[["normalised.ABP.signal"]][dans[["transient.ABP.samples"]]]
  dans[["transient.ABP.slope"]] <- mean(dans[["transient.normalised.ABP.signal"]]) / mean(dans[["transient.ABP.time.instants"]] - dans[["min.ABP.time.instant"]])
  dans[["transient.ABP.line"]] <- dans[["transient.ABP.slope"]] * (time.instants[dans[["transient.ABP.samples"]]] - dans[["min.ABP.time.instant"]])
  
  # Gets transient ABP angle
  dans[["transient.ABP.angle"]] <- atan(dans[["transient.ABP.slope"]]) * 180 / pi
  dans[["bounded.transient.ABP.angle"]] <- dans[["transient.ABP.angle"]]
  if(dans[["bounded.transient.ABP.angle"]] < min.ABP.angle)
    dans[["bounded.transient.ABP.angle"]] <- min.ABP.angle
  if(dans[["bounded.transient.ABP.angle"]] > max.ABP.angle)
    dans[["bounded.transient.ABP.angle"]] <- max.ABP.angle
  
  # Gets transient CBFV angle
  dans[["transient.CBFV.angle"]] <- atan(dans[["transient.CBFV.slope"]]) * 180 / pi
  dans[["bounded.transient.CBFV.angle"]] <- dans[["transient.CBFV.angle"]]
  if(dans[["bounded.transient.CBFV.angle"]] < min.CBFV.angle)
    dans[["bounded.transient.CBFV.angle"]] <- min.CBFV.angle
  if(dans[["bounded.transient.CBFV.angle"]] > max.CBFV.angle)
    dans[["bounded.transient.CBFV.angle"]] <- max.CBFV.angle
  
  # Gets Phi parameter
  
  dans[["transient.angles.difference"]] <- dans[["transient.CBFV.angle"]] - dans[["transient.ABP.angle"]]
  dans[["bounded.transient.angles.difference"]] <- dans[["bounded.transient.CBFV.angle"]] - dans[["bounded.transient.ABP.angle"]]
  dans[["Phi"]] <- dans[["bounded.transient.angles.difference"]]
  if(dans[["Phi"]] < min.Phi)
    dans[["Phi"]] <- min.Phi
  if(dans[["Phi"]] > max.Phi)
    dans[["Phi"]] <- max.Phi
  dans[["Phi"]] <- round(dans[["Phi"]], Phi.rounding.digits)
  
  if(!keep.details)
  {
    ans <- list()
    
    ans[["ABP.baseline.samples"]] <- dans[["ABP.baseline.samples"]]
    ans[["ABP.baseline.value"]] <- dans[["ABP.baseline.value"]]
    
    ans[["min.ABP.samples"]] <- dans[["min.ABP.samples"]]
    ans[["min.ABP.sample"]] <- dans[["min.ABP.sample"]]
    ans[["min.ABP.value"]] <- dans[["min.ABP.value"]]
    ans[["min.ABP.type"]] <- dans[["min.ABP.type"]]
    ans[["min.ABP.drop.pc"]] <- dans[["min.ABP.drop.pc"]]
    
    ans[["CBFV.baseline.samples"]] <- dans[["CBFV.baseline.samples"]]
    ans[["CBFV.baseline.value"]] <- dans[["CBFV.baseline.value"]]
    
    ans[["min.CBFV.samples"]] <- dans[["min.CBFV.samples"]]
    ans[["min.CBFV.sample"]] <- dans[["min.CBFV.sample"]]
    ans[["min.CBFV.value"]] <- dans[["min.CBFV.value"]]
    ans[["min.CBFV.type"]] <- dans[["min.CBFV.type"]]
    
    ans[["normalised.ABP.signal"]] <- dans[["normalised.ABP.signal"]]
    ans[["normalised.CBFV.signal"]] <- dans[["normalised.CBFV.signal"]]
    
    ans[["transient.ABP.samples"]] <- dans[["transient.ABP.samples"]]
    ans[["transient.ABP.slope"]] <- dans[["transient.ABP.slope"]]
    ans[["transient.ABP.line"]] <- dans[["transient.ABP.line"]]
    
    ans[["transient.CBFV.samples"]] <- dans[["transient.CBFV.samples"]]
    ans[["transient.CBFV.slope"]] <- dans[["transient.CBFV.slope"]]
    ans[["transient.CBFV.line"]] <- dans[["transient.CBFV.line"]]
    ans[["steady.CBFV.samples"]] <- dans[["steady.CBFV.samples"]]
    ans[["steady.CBFV.line"]] <- dans[["steady.CBFV.line"]]
    ans[["steady.CBFV.duration"]] <- dans[["steady.CBFV.duration"]]
    
    ans[["Delta.tau"]] <- dans[["Delta.tau"]]
    ans[["Ks"]] <- dans[["Ks"]]
    ans[["Phi"]] <- dans[["Phi"]]
  }
  else
    ans <- dans
  
  invisible(ans)
}

#Calculo de mfari para VE
#Calculo de mfari sin phi porque solo se utiliza para comparar en caso de maniobras no con VE
get.mfARI <- function(Ks, Delta.tau, Phi)
{
  #calculo de mfari con phi
  invisible(1.48740 + 3.85688 * Ks - 0.12420 * Delta.tau + 0.10999 * Phi)
  #Calculo de mfari sin phi
  #invisible(6.75148 + 2.51501 * Ks - 0.61895 * Delta.tau)
}



############# ARI ##############
get.theoretical.CBFV.response <- function(
  T = 1.9,
  D = 0.75,
  K = 0.9,
  time.instants,
  ABP.normalised,
  sampling.time = min(round(diff(time.instants), 3)),
  stabilisation.time = 1
)
{
  # Initialises the answer
  ans <- list()
  ans[["T"]] <- T
  ans[["D"]] <- D
  ans[["K"]] <- K
  ans[["time.instants"]] <- time.instants
  ans[["sampling.time"]] <- sampling.time
  ans[["ABP.normalised"]] <- ABP.normalised
  
  frequency <- 1 / ans[["sampling.time"]]
  nsamples <- length(ans[["time.instants"]])
  nsamples.stabilisation <-
    round(stabilisation.time / ans[["sampling.time"]])
  
  P <- c(
    rep(ans[["ABP.normalised"]][1], nsamples.stabilisation),
    ans[["ABP.normalised"]],
    rep(ans[["ABP.normalised"]][nsamples], nsamples.stabilisation)
  )
  
  # Gets dP
  dP <- P - 1
  
  # Applies Tiecks' equations to obtain the CBFV signal
  X1 <- vector(mode = "numeric", length = length(P))
  X2 <- vector(mode = "numeric", length = length(P))
  CBFV <- vector(mode = "numeric", length = length(P))
  
  divisor <- frequency * ans[["T"]]
  X1[1] <- 0
  X2[1] <- 0
  CBFV[1] <- 1
  for(t in 2:length(P))
  {
    X1[t] <- X1[t-1] + (dP[t] - X2[t-1]) / divisor
    X2[t] <- X2[t-1] + (X1[t] - 2 * ans[["D"]] * X2[t-1]) / divisor
    CBFV[t] <- 1 + dP[t] - ans[["K"]] * X2[t]
  }
  
  if(nsamples.stabilisation > 0)
    CBFV <- CBFV[-(1:nsamples.stabilisation)]
  CBFV <- CBFV[1:nsamples]
  ans[["CBFV.theoretical.response"]] <- CBFV
  
  invisible(ans)
}

get.AT.templates.parameters <- function()
{
  K <- c(0.00, 0.20, 0.40, 0.60, 0.80, 0.90, 0.94, 0.96, 0.97, 0.98)
  D <- c(1.60, 1.60, 1.50, 1.15, 0.90, 0.75, 0.65, 0.55, 0.52, 0.50)
  T <- c(2.00, 2.00, 2.00, 2.00, 2.00, 1.90, 1.60, 1.20, 0.87, 0.65)
  ARI <- 0:9
  
  data.frame(T, D, K, ARI)
}

get.AT.decimal.templates.parameters <- function(rounding.digits = 6)
{
  orig <- get.AT.templates.parameters()
  
  ARI.decimal <- round(seq(0, 9, 0.1), 1)
  K.decimal <-
    pracma::interp1(orig[["ARI"]], orig[["K"]], ARI.decimal, 'spline')
  K.decimal <- round(K.decimal, rounding.digits)
  D.decimal <- 
    pracma::interp1(orig[["ARI"]], orig[["D"]], ARI.decimal, 'spline')
  D.decimal <- round(D.decimal, rounding.digits)
  T.decimal <-
    pracma::interp1(orig[["ARI"]], orig[["T"]], ARI.decimal, 'spline')
  T.decimal <- round(T.decimal, rounding.digits)
  
  data.frame(
    T = T.decimal,
    D = D.decimal,
    K = K.decimal,
    ARI = ARI.decimal
  )
  print(data.frame(
    T = T.decimal,
    D = D.decimal,
    K = K.decimal,
    ARI = ARI.decimal
  ))
}

get.normalised.dari.templates <- function(
  time.instants,
  normalised.abp.signal,
  sampling.time,
  sample.release,
  baseline.initial.time = time.instants[1],
  baseline.final.time = time.instants[sample.release],
  min.cbfv.max.delta.time = 8 * 0.8,
  stabilisation.time = 30,
  at.param.rounding.digits = 6,
  cbfv.rounding.digits = 2,
  time.tol = sampling.time / 100,
  ... # Not used
)
{
  at.params <- get.AT.decimal.templates.parameters(
    rounding.digits = at.param.rounding.digits
  )
  
  .tmp.fun <- function(i)
    get.theoretical.CBFV.response(
      T = at.params[i, 1],
      D = at.params[i, 2],
      K = at.params[i, 3],
      time.instants = time.instants,
      ABP.normalised = normalised.abp.signal,
      sampling.time = sampling.time,
      stabilisation.time = stabilisation.time
    )
  r <- 1:nrow(at.params)
  templates <- lapply(r, .tmp.fun)
  
  .tmp.fun <- function(t)
    normalise.CBFV.signal(
      time.instants = time.instants,
      CBFV.signal = t[["CBFV.theoretical.response"]],
      sample.release = sample.release,
      sampling.time = sampling.time,
      baseline.initial.time = baseline.initial.time,
      baseline.final.time = baseline.final.time,
      min.CBFV.max.delta.time = min.cbfv.max.delta.time,
      time.tol = time.tol
    )
  normalised.templates <- lapply(templates, .tmp.fun)
  
  .tmp.fun <- function(t) t[["normalised.CBFV.signal"]]
  normalised.templates <- lapply(normalised.templates, .tmp.fun)
  
  names(normalised.templates) <- 
    sapply(r, function(i) sprintf("%.1f", at.params[i, 4]))
  
  normalised.templates
}


get.best.templates <- function(
  time.instants,
  signal,
  templates,
  referential.time.instant = 0,
  delta.time.before.ref = 0,
  delta.time.after.ref = 20 * 0.8,
  comparison.function = get.MSE,
  keep.details = TRUE,
  time.tol = min(diff(time.instants)) / 100,
  ... # Pass over to comparison.function()
)
{
  # Validates lengths
  lengths <- c(length(time.instants), length(signal))
  lengths <- c(lengths, sapply(templates, length))
  lengths <- unique(lengths)
  if(length(lengths) != 1)
    stop("time, signal and templates must have the same length")
  
  # Initialises detailed answer
  ans <- list()
  ans[["time.instants"]] <- time.instants
  ans[["signal"]] <- signal
  ans[["templates"]] <- templates
  ans[["referential.time.instant"]] <- referential.time.instant
  
  # Finds referential sample
  i <- which(are.tolerably.equal(
    time.instants,
    referential.time.instant,
    time.tol
  ))
  if(length(i) != 1)
    stop("a unique referential time instant could not be determined")
  ans[["referential.sample"]] <- i
  
  # Finds initial sample
  ans[["delta.time.before.ref"]] <- delta.time.before.ref
  ans[["initial.time.instant"]] <- referential.time.instant - delta.time.before.ref
  i <- which(are.tolerably.equal(
    time.instants,
    ans[["initial.time.instant"]],
    time.tol
  ))
  if(length(i) != 1)
    stop("initial time instant could not be found in specified time instants")
  ans[["initial.sample"]] <- i
  
  # Finds final sample
  ans[["delta.time.after.ref"]] <- delta.time.after.ref
  ans[["final.time.instant"]] <-
    referential.time.instant + delta.time.after.ref
  i <- which(are.tolerably.equal(
    time.instants,
    ans[["final.time.instant"]],
    time.tol
  ))
  if(length(i) != 1)
    stop("final time instant could not be found in specified time instants")
  ans[["final.sample"]] <- i
  
  # Sets relevant interval
  ans[["relevant.samples"]] <-
    ans[["initial.sample"]]:ans[["final.sample"]]
  
  # Gets relevant segments
  .tmp.fun <- function(s) s[ans[["relevant.samples"]]]
  ans[["relevant.signal.segment"]] <- .tmp.fun(ans[["signal"]])
  ans[["relevant.template.segments"]] <-
    lapply(ans[["templates"]], .tmp.fun)
  
  # Gets fit values
  .tmp.fun <- function(t) 
    comparison.function(ans[["relevant.signal.segment"]], t, ...)
  ans[["fit.values"]] <-
    sapply(ans[["relevant.template.segments"]], .tmp.fun)
  
  # Orders templates and determines the best one
  ans[["ranking"]] <- order(ans[["fit.values"]])
  
  # Deletes details if corresponds
  if(!keep.details)
  {
    i <- ans[["ranking"]][1]
    e <- ans[["fit.values"]][i]
    ans <- list(best.template.index = i, best.fit.value = e)
  }
  
  ans
}




################################### Programa principal #################################
run<-function()
{
  ### escalon
  filepath<-("C:/Users/Usuario/Google Drive/USACH/TESIS/programas")
  setwd(filepath)

  # #crear archivo escalon
  # x <- seq(pi/2,pi,length.out=100)
  # y <- sin(x)
  # nombrearch <- paste ("seno", "txt", sep = ".", collapse = NULL)
  # write.table(y, file = nombrearch, sep = "\t", row.names= FALSE)
  
  ## leer escalon
  filename <- paste ("Esc_neg", "txt", sep = ".", collapse = NULL)
  escalon<- read.table(filename, header=TRUE)
  inverseStep<-escalon[1:350,]
  w<-0.2
  filbut<- butter(2,w)
  inverseStep<-filtfilt(filbut, inverseStep)
  inverseStep<-inverseStep[20:200]
  print("lectura archivo")
  filepath <- ("C:/Users/Usuario/Google Drive/USACH/TESIS/programas/cluster_VE5mins_NARX_p1/cluster/Results/Univariado/NARX_grilla1/")
  setwd(filepath)
  subjects = c( #"P046-VE_",
                # 'ANGL-VE_','PAAR-VE_',
                # 'P040-VE_',
                # 'P041-VE_','P061-VE_',
                # 'P044-VE_',
                # 'P060-VE_',
                # 'P045-VE_','P052-VE_',
                # 'P048-VE_','P062-VE_',
                # 'P049-VE_',
                'P051-VE_'#,
                # 'P050-VE_','P057-VE_',
                # 'P056-VE_','P066-VE_',
                # 'P063-VE_',
                # 'P064-VE_',
                #'P067-VE_'
              )
  
  for (subject in subjects){
    namesubject <- subject
    filename <- paste (namesubject, "txt", sep = ".", collapse = NULL)
    archivo <- read.table(filename, header=TRUE)
    mejoresModelos1<-archivo[order(archivo[,7], decreasing = TRUE),]
    
    ##parametros
    gamma<-mejoresModelos1$gamma
    nu <- mejoresModelos1$nu
    cost <- mejoresModelos1$cost
    fold<-mejoresModelos1$fold
    PAMn<-0
    VFSCn<-0
    Ts=0.2
    
    mejoresModelos<-mejoresModelos1
    
    ### evaluacion de cada modelo
    for (i in 1:length(mejoresModelos[,1])){
      filepath<-("C:/Users/Usuario/Google Drive/USACH/TESIS/programas/cluster_VE5mins_NARX_p1/cluster/Data")
      setwd(filepath)
      if(fold[i]==1){
        namearch <- paste (namesubject, "02izq", sep = "_", collapse = NULL)
        filename <- paste (namearch, "txt", sep = ".", collapse = NULL)
        fold1 <- read.table(filename, header=TRUE)
        names(fold1)<-c("PAMn","VFSCn")
        attach(fold1)
        PAMn <- fold1[1]
        VFSCn <- fold1[2]
      } else{
        namearch <- paste (namesubject, "01izq", sep = "_", collapse = NULL)
        filename <- paste (namearch, "txt", sep = ".", collapse = NULL)
        fold2 <- read.table(filename, header=TRUE)
        names(fold2)<-c("PAMn","VFSCn")
        attach(fold2)
        PAMn <- fold2[1]
        VFSCn <- fold2[2]
      }
      
      #modelo
      PAMn<-normalize(PAMn)
      VFSCn<-normalize(VFSCn)
      data <- data.frame(PAMn,VFSCn)
      lag<-list(PAMn = mejoresModelos[i,1],VFSCn = 0)
      signal.train <- retardos_multi(data, lag)
      retDatos=signal.train$folded.signal
      x=subset(retDatos, select = -VFSCn)
      x=subset(x, select = -PAMn)
      y=retDatos$VFSCn
      mejorModelo <- svm(x, y, kernel = "radial",type = "nu-regression",
                         cost = cost[i], nu = nu[i], gamma=gamma[i])
      
      ### respuesta al escalon
      PAMn=normalize(inverseStep)
      VFSCn=normalize(inverseStep)
      data <- data.frame(PAMn,VFSCn)
      lag<-list(PAMn = mejoresModelos[i,1],VFSCn = 0)
      signal.train <- retardos_multi(data, lag)
      retDatos=signal.train$folded.signal
      x=subset(retDatos, select = -VFSCn)
      x=subset(x, select = -PAMn)
      y=retDatos$VFSCn
      stepTime=seq(Ts,(length(retDatos$PAMn))*Ts,Ts)
      stepResponse <- predict(mejorModelo, x )
      
      pamnormalizado<-normalize(retDatos$PAMn)
      respuestanormalizada<-normalize(stepResponse)
      
      tiempoManiobra<- 0
      j<-0
      for (variable in stepTime) {
        tiempoManiobra[j]<- variable-6.2
        j<-j+1
      }
      
      # Criterios de evaluacion del modelo
      # si cumple con los criterios entonces graficar y calcular mfari
      beforeDrop=5;
      iniCaida=beforeDrop/Ts
      finCaida=(beforeDrop+6)/Ts

      caida = min(respuestanormalizada[iniCaida:finCaida])
      plateau=max(respuestanormalizada[(beforeDrop-2)/Ts:(beforeDrop+2)/Ts])
      #calculo de la varianza
      varPlateau=var(respuestanormalizada[1:(beforeDrop)/Ts])

      #grafico del escalon y la respuesta
      largo.mvre<-length(tiempoManiobra)
      pamnormalizado<-pamnormalizado[1:largo.mvre]
      respuestanormalizada<-respuestanormalizada[1:largo.mvre]
      
      # plot(tiempoManiobra,pamnormalizado,type="l", col="red")
      # lines(tiempoManiobra,respuestanormalizada, col = "blue")
      
      #calculo de mfari y sus parametros
      time.instants<-tiempoManiobra
      abp.signal<-pamnormalizado
      cbfv.signal<-respuestanormalizada
      sampling.time<-0.2
      time.release<-0
      
      mfariParams <- get.mfARI.parameters(
        time.instants = time.instants,
        ABP.signal = abp.signal,
        CBFV.signal = cbfv.signal,
        sampling.time = sampling.time,
        time.release = time.release,
        keep.details = FALSE
      )
      
      mfari<- get.mfARI(
        Ks = mfariParams[["Ks"]],
        Delta.tau = mfariParams[["Delta.tau"]],
        Phi = mfariParams[["Phi"]]
      )
      cbfvmin<-mfariParams[["min.CBFV.value"]]
      ks<- mfariParams[["Ks"]]
      delta.tau<-mfariParams[["Delta.tau"]]
      phi<- mfariParams[["Phi"]]
      
      #calculo del ARI
      time.tol<-sampling.time/100
      sample.release <- which(are.tolerably.equal(
        time.instants,
        time.release,
        time.tol
      ))
      normalised.cbfv.templates<-get.normalised.dari.templates(
        time.instants = time.instants,
        normalised.abp.signal = abp.signal,
        sampling.time = sampling.time,
        sample.release = sample.release,
        time.tol = time.tol
      )
      
      ari<-get.best.templates(
        time.instants = time.instants,
        signal = cbfv.signal,
        templates = normalised.cbfv.templates,
        keep.details = TRUE
      )
      # i <- seq(11, 91, 20)
      # temps <- ari[["relevant.template.segments"]][i]
      ibest <- ari[["ranking"]][1]
      best <- ari[["relevant.template.segments"]][[ibest]]
      
      #que solo me grafique cuando encuentre un ks>0
      if(delta.tau>=3 && ks>=cbfvmin){
      
      plot(tiempoManiobra,pamnormalizado,type="l", col="red",xlab = "Tiempo (seg)", ylab = "VFSC Estimado Normalizado (cm/seg)", main = "Respuesta al escalón inverso de presión")
      lines(tiempoManiobra,respuestanormalizada, col = "blue")
    
      legend("topright", 
             c("Escalon de presión", "respuesta al escalon"), 
             title = "Autorregulacion", 
             pch =1, 
             col=c("red","blue"),
             lty=c(1,1),
             inset = 0.01)
      legend("bottomright", 
             c(paste("corr=",round(mejoresModelos[i,7],2)),paste("Ks=",ks),paste("Delta.tau=",delta.tau), paste("Phi=",phi),paste("mfARI=",mfari)), 
             title = "Parametros mfARI")
      
      #print(ari)
      print(best)
      readline(prompt="Press [enter] to continue")
      }
    }
  }
}
