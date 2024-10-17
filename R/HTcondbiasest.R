

#' HTcondbiasest
#' The function to estimate the conditional biases for the HT estimator of a total, for the non-stratiﬁed sampling designs. For one or several given variables of interest as input, this function returns a vector with the estimates of the conditional bias of the HT estimator where each row corresponds to a sample unit.
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param di  Inverse of the first order inclusion probabilities
#' @param id  Primary key to keep as an identifier of every unit
#' @param remerge  True/False to remerge the resulting conditional biases with the original data set
#'
#' @return Returns a data set with the conditional bias of each sampled unit
#' @export
#'

"HTcondbiasest" <- function (data,
                             varname = NULL,
                             gn,
                             method = c("si", "poisson", "rejective"),
                             pii = NULL,
                             di = NULL,
                             id = NULL,
                             remerge = T) {
  if (nrow(data) <= 1) {
    stop("your sample must contain at least 2 units\n")
  }
  if (sum(is.na(data[,varname])) != 0) {
    stop("missing values for some y-variables\n")
  }
  if (any(data[,varname] < 0)) {
    warning("Warning: your variable(s) of interest contain negative values\n")
  }
  if (missing(method) | length(method) > 1) {
    warning("Warning: the method is not specified or multiple methods are specified.\nBy default, the method is 'si'\n")
    method = "si"
  }
  if (!(method %in% c("si", "poisson", "rejective"))) {
    stop("the name of the method is wrong\n")
  }
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
  } else if (!missing(pii) & !missing(di)) {
    warning("Warning: di is redundant. Only pii is being used\n")
  }
  if (missing(gn)) {
    stop("the population size is missing\n")
  }
  if (method == "si" & !isTRUE(all.equal(gn, sum(1/pii)))) {
    warning("Warning: the sum of the inclusion probabilities is not equal to N\n")
    print(paste("N:", gn))
    print(paste("sum(1/pii):", sum(1/pii)))
  }
  if (remerge==F) {
    if (missing(id)) {
      warning("Warning: no column is specified as an identifier, one is added automatically (id)\n")
      id <- "id"
      identifier <- c(1:nrow(data))
    } else {
      if (id != "none") {
        if (!(id %in% colnames(data))) {
          stop("the specified identifier is not a column name\n")
        }
        if (sum(duplicated(data[,id]))>1) {
          stop("there are dupkeys on the id variable")
        }
        identifier <- data[[id]]
      }
    }
  }

  data = data.frame(data)
  pn = nrow(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  data2 = cbind.data.frame(data[, m], index)
  colnames(data2) = c(varname, "index")
  if (method == "si") {
    if (length(m)==1){
      bc = (pn/(pn-1))*(gn/pn-1)*(data[,m]-mean(data[,m]))
    } else {
      bc = (pn/(pn-1))*(gn/pn-1)*(data[,m]-matrix(data=colMeans(data[,m]),nrow=pn,ncol=length(m), byrow =T))
    }
  }
  if (method == "poisson") {
    bc = (1/pii-1) * data[,m]
  }
  if (method == "rejective") {
    gd = sum(1-pii)
    gb = t(1/pii-1) %*% as.matrix(data[,m]) / gd
    print(paste0("D: ", round(gd, 3)))
    print(paste0("N/D: ", round(gn/gd, 3)))
    warning("Warning: please make sure that D is large enough and N/D is bounded\n")
    bc = (1/pii - 1) * (data[,m] - t(t(gb) %*% c(pii)))
  }
  if (remerge == T) {
    result = cbind.data.frame(data, bc)
    colnames(result) = c(colnames(data), paste0("condbias", colnames(bc)))
  } else {
    result = data.frame(bc)
    if (length(m) == 1) {
      colnames(result) = c("condbias")
    } else {
      colnames(result) = c(paste0("condbias", colnames(bc)))
    }
    if (id != "none") {
      result <- cbind.data.frame(identifier, result)
      colnames(result)[1] <- id
    }
  }
  result
}


