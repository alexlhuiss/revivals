
#' strata_robustest
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Population size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#'
#' @return Computes the robust estimator of Beaumont et al.(2013) by stratum using the conditional bias and the minmax criterion to compute the tuning constant
#' @export

"strata_robustest" <- function(data,
                               strataname = NULL,
                               varname = NULL,
                               gnh,
                               method = c("si","poisson","rejective"),
                               pii) {
  if (missing(gnh)) {
    stop("the population size vector is missing\n")
  }
  if (missing(strataname) | is.null(strataname)) {
    stop("no variable name to use for stratification has been specified\n")
  }
  data = data.frame(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  ms = match(strataname, colnames(data))
  if (any(is.na(ms))) {
    stop("the name of the strata is wrong\n")
  }
  x1 = data.frame(unique(data[,ms]))
  rht = matrix(0, nrow=nrow(x1), ncol=length(m))

  for (i in 1:nrow(x1)) {
    print(paste("Stratum ", i, ":", sep=""))
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    rht[i,] = robustest(data=datastr, varname=varname, gn=gnh[i], method=method, pii=piisrt)
  }
  result = cbind.data.frame(rht)
  colnames(result) = c(paste0("RHT_", colnames(datastr)))
  rownames(result) = c(paste("Stratum", 1:nrow(x1)))
  print("Done!")
  result
}
