#' Tidy method for \code{clm} objects.
#'
#' This method tidies the coefficients of a \code{clm} object into a data frame.
#'
#' If you have missing values in your model data, you may need to refit the model with
#' \code{na.action = na.exclude}.
#'
#' @param x clm object
#' @param conf.int Logical: whether to include a confidence interval
#' @param conf.level confidence level of the interval, used only if \code{conf.int = TRUE}
#' @param exponentiate Logical: whether to exponentiate the coefficient estimates and confidence intervals
#' @param quick Logical: whether to compute a smaller and faster version
#' @param cluster Logical: whether to use clustered standard errors
#' @param vcov an alternative variance/covariance matrix, used only if \code{cluster = TRUE}
#' @param ... extra arguments
#'
#' @return tidy data frame
#' @export
#'
#' @importFrom broom tidy
#'
#' @examples
#' \dontrun{
#' require(ordinal)
#' fm1 <- clm(rating ~ temp * contact, data = wine)
#' tidy(fm1)
#' }
tidy.clm <- function (x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE,
                      quick = FALSE, cluster = FALSE, vcov = NULL, ...) {
  if (quick) {
    co <- coef(x)
    ret <- data.frame(term = names(co), estimate = unname(co))
    return(process_clm(ret, x, conf.int = FALSE, exponentiate = exponentiate))
  }
  s <- summary(x)
  ret <- tidy.summary.clm(s, cluster = cluster, vcov = vcov)
  process_clm(ret, x, conf.int = conf.int, conf.level = conf.level,
              exponentiate = exponentiate)
}

process_clm <- function(ret, x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE){
  if (exponentiate){
    trans <- exp
  } else {
    trans <- identity
  }
  if (conf.int) {
    CI <- suppressMessages(confint(x, level = conf.level, type = 'Wald'))
    colnames(CI) <- c('conf.low', 'conf.high')
    ret <- cbind(ret, trans(unrowname(CI)))
  }
  ret$estimate <- trans(ret$estimate)
  ret
}

unrowname <- function(x){
  rownames(x) <- NULL
  x
}

#' Tidy method for summaries of \code{clm} objects
#'
#' @param x clm object
#' @param cluster Logical: whether to calculate clustered standard errors
#' @param vcov alternative variance/covariance matrix, only used if \code{cluster = TRUE}
#' @param digits  digits to print
#' @param ... additional arguments
#'
#' @return tidy data frame
#' @export
#'
#' @examples
#' \dontrun{
#' require(ordinal)
#' fm1 <- clm(rating ~ temp * contact, data = wine)
#' tidy(summary(fm1))
#' }
tidy.summary.clm <- function(x, cluster = FALSE, vcov = NULL, digits = 2, ...){
  requireNamespace("plyr", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)

  co <- coef(x)
  nn <- c("estimate", "std.error", "statistic", "p.value")
  if (inherits(co, "listof")) {
    ret <- plyr::ldply(co, fix_data_frame, nn[1:ncol(co[[1]])],
                       .id = "response")
    ret$response <- stringr::str_replace(ret$response, "Response ",
                                         "")
  }
  else {
    ret <- fix_data_frame(co, nn[1:ncol(co)])
  }
  if (cluster){
    if (is.matrix(vcov)) {
      vcov <- sqrt(diag(vcov))
    }
    ret[, 'std.error'] <- vcov
    ret[, 'statistic'] <- ret[, 'estimate'] / vcov
    ret[, 'p.value'] <- pnorm(abs(ret[, 'statistic', drop = TRUE]), lower.tail = FALSE) * 2
  }
  ret
}

#' Glance method for \code{clm} objects
#'
#' @param x clm object
#' @param ... extra arguments
#'
#' @return data frame
#' @export
#'
#' @importFrom broom glance
#'
#' @examples
#' \dontrun{
#' require(ordinal)
#' fm1 <- clm(rating ~ temp * contact, data = wine)
#' glance(fm1)
#' }

glance.clm <- function(x, ...){
  ret <- as.data.frame(x[c("logLik")])
  rownames(ret) <- NULL
  finish_glance(ret, x)
}

#' Clustered standard errors for \code{clm}
#'
#' @param model clm object
#' @param cluster vector identifying clusters
#'
#' @return variance/covariance matrix
#' @export
get_CL_vcov<-function(model, cluster){

  if (length(model$na.action) > 0) {
    cluster <- cluster[-model$na.action]
  }

  #calculate degree of freedom adjustment
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- length(coef(model))
  dfc <- (M / (M - 1)) * ((N - 1) / (N - K))

  #calculate the uj's
  uj  <- apply(sandwich::estfun(model), 2, function(x) tapply(x, cluster, sum))

  #use sandwich to get the var-covar matrix
  vcovCL <- dfc * sandwich::sandwich(model, meat = crossprod(uj) / N)
  return(vcovCL)
}

get_confint<-function(model, vcovCL){
  t <- qt(.975,  model$df.residual)
  ct <- coeftest(model, vcovCL)
  est <- cbind(ct[, 1], ct[, 1] - t * ct[, 2], ct[, 1] + t * ct[, 2])
  colnames(est) <- c("Estimate","LowerCI","UpperCI")
  return(est)
}
