
#' robustest
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#'
#' @return Computes the robust estimator of Beaumont et al.(2013) using the conditional bias and the minmax criterion to compute the tuning constant
#' @export
#'

"robustest" <- function (data,
                         varname = NULL,
                         gn,
                         method = c("si","poisson","rejective"),
                         pii) {
  data = data.frame(data)
  pn = nrow(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  data2 = cbind.data.frame(data[, m], index)
  colnames(data2) = c(varname, "index")
  if (any(is.na(data[,m]))) {
    stop("missing values for some y-variables\n")
  }
  if (length(m) == 1) {
    htestim = crossprod(data[,m], 1/pii)
    htbc = HTcondbiasest(data=data, varname=varname, gn=gn, method=method, pii=pii, id="none")[,c(seq(1+length(data),length(data)+length(varname)))]
    result = htestim - (min(htbc)+max(htbc))/2
  } else {
    htestim = (1/pii) %*% as.matrix(data[,m])
    htbc = HTcondbiasest(data=data, varname=varname, gn=gn, method=method, pii=pii, id="none")[,c(seq(1+length(data),length(data)+length(varname)))]
    result = htestim-(apply(X=htbc,MARGIN = 2,min )+apply(X =htbc,MARGIN = 2,max ))/2
  }
  colnames(result) = c(paste0("RHT_", colnames(htbc)))
  result
}

