
#' strata_HTcondbiasest
#' The function to estimate the conditional biases for the HT estimator of a total, for the stratified sampling designs. For one or several given variables of interest as input, this function returns a vector with the estimates of the conditional bias of the HT estimator where each row corresponds to a sample unit.
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Population size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param di  Inverse of the first order inclusion probabilities
#' @param id  Primary key to keep as an identifier of every unit
#' @param remerge  True/False to remerge the resulting conditional biases with the original data set
#'
#' @return Returns a dataframe with the conditional bias of each sampled unit
#' @export


"strata_HTcondbiasest" <- function(data,
                                   strataname = NULL,
                                   varname = NULL,
                                   gnh,
                                   method = c("si","poisson","rejective"),
                                   pii = NULL,
                                   di = NULL,
                                   id = NULL,
                                   remerge = T) {
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
  } else if (!missing(pii) & !missing(di)) {
    warning("Warning: di is redundant. Only pii is being used\n")
  }
  if (missing(gnh)) {
    stop("the population size vector is missing\n")
  }
  if (missing(strataname) | is.null(strataname)) {
    stop("no variable name to use for stratification has been specified\n")
  }
  if (any(table(data[,strataname]) <= 1)) {
    stop("each stratum must contain at least 2 units\n")
  }
  if (remerge==F) {
    if (missing(id)) {
      warning("Warning: no column is specified as an identifier, one is added automatically (id)\n")
      id <- "id"
      identifier <- c(1:nrow(data))
    } else {
      if (!(id %in% colnames(data))) {
        stop("the specified identifier is not a column name\n")
      }
      if (sum(duplicated(id))>1) {
        stop("there are dupkeys on the id variable")
      }
      identifier <- data[[id]]
    }
  }

  data = data.frame(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m)))
    stop("the name of the variable is wrong\n")
  ms = match(strataname, colnames(data))
  if (any(is.na(ms)))
    stop("the name of the strata is wrong\n")
  x1 = data.frame(unique(data[,ms]))
  bc = matrix(0, nrow=nrow(data), ncol=length(m))
  cgn = 0

  for (i in 1:nrow(x1)) {
    print(paste("Stratum ", i, ":", sep=""))
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    resint = HTcondbiasest(data=datastr, varname=varname, gn=gnh[i], method=method, pii=piisrt, id="none", remerge=F)
    bc[(cgn+1):(cgn+nh),] = as.matrix(resint)
    # bc[(cgn+1):(cgn+nh)]=nh[i]/(nh[i]-1)*(gnh[i]/nh-1)*(y-mean(y))
    cgn = cgn+nh
  }
  if (remerge==T) {
    result = cbind.data.frame(data, bc)
    colnames(result) = c(colnames(data), colnames(resint))
  } else {
    result = data.frame(bc)
    if (length(m) == 1) {
      colnames(result) = c("condbias")
    } else {
      colnames(result) = c(paste0("condbias", varname))
    }
    result <- cbind.data.frame(identifier, result)
    colnames(result)[1] <- id
  }
  print("Done!")
  result
}

