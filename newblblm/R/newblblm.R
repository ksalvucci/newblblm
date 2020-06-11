
## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' @aliases NULL
#' @import purrr
#' @import stats
#' @import furrr
#' @import tidyverse
#' @import ggplot2
#' @import future
#' @importFrom magrittr %>%
#' @importFrom graphics boxplot abline plot
#' @importFrom utils capture.output
#' @importFrom grDevices rgb
"_PACKAGE"
#' Linear Regression with Little Bag of Bootstraps
#' @details This model uses blb to find a lm fit.
#' @param formula a formula: a symbolic description of the model to be fitted.
#' @param data a data.frame or a vector list of file names.
#' @param m an integer for the number of subsets you want your data to be split into.
#' @param B an integer for the number of bootstraps you would like to use on each subset.
#' @return returns an object of class "newblblm".
#' The applicable generic functions include boxplot, coef, confint, plot, predict, print, and sigma.
#' @export
#' @examples
#' fit <- newblblm(mpg ~ wt, data = mtcars, m = 3, B = 1000)

newblblm <- function(formula, data, m = 10, B = 5000) {
  estimates <- future_estimate(formula = formula, data = data, m = m, B = B)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "newblblm"
  invisible(res)
}


#' find the estimates
#' @param formula a symbol description of the model to be fitted
#' @param data a data.frame of a vector list of file names
#' @param m integer. the number of subsamples
#' @param B integer. the number of bootstraps for each subsample
future_estimate <- function(formula, data, m, B) {
  if (class(data) == "data.frame") {
    data_list <- split_data(data, m)
    future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
    )
  } else if (class(data) == "character") {
    data %>% future_map(~ {
      df <- read_csv(., col_types = cols())
      lm_each_subsample(formula = formula, data = df, n = nrow(df), B = B)
    })
  }
}



#' split data into m parts of approximated equal sizes
#' @param data data.frame to split
#' @param m integer. the number of subsamples to make
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#' @param formula a symbolic description of the model to be fitted
#' @param data the data the model will be fitted on
#' @param n the size of the original data
#' @param B the number of bootstraps to be created
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#' @param formula a symbolic description of the model to be fitted
#' @param data the data the model will be fitted on
#' @param n the size of the original data
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @param formula a symbolic description of the model to be fitted
#' @param data the data the model will be fitted on
#' @param freqs the frequency of each observation in the dataset
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @param fit a class newblblm object
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit a class newblblm object
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print newblblm
print.newblblm <- function(x, ...) {
  cat("newblblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma newblblm
sigma.newblblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - level
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef newblblm
coef.newblblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint newblblm
confint.newblblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}



#' @export
#' @method predict newblblm
predict.newblblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}


#' Generic function fitplot
#' @param object the model
#' @param data data.frame to be plotted
#' @param boots logical
#'
#' @export
#'
fitplot <- function(object, data, boots = FALSE) {
  UseMethod("fitplot")
}

#' Plot a newblblm Linear Model
#' @description Plot the data alongside the linear regression line computed from the newblblm procdure.
#' @param object a newblblm class object
#' @param data data.frame. The original data that the model is fitted on
#' @param boots logical. The option to show the linear regression line for each bootstrap
#'
#' @export
#' @method fitplot newblblm
#'
fitplot.newblblm <- function(object, data, boots = FALSE) {
x <- attr(terms(object$formula), "term.labels")
if (length(x) > 1) {
  print("Must be a two-dimensional model")
} else {
  y <- as.character(object$formula)[2]
  datax <- data[x] %>% pluck(1)
  datay <- data[y] %>% pluck(1)
  intercept <- coef(object) %>% pluck(1)
  slope <- coef(object) %>% pluck(2)
  plot(datax, datay, xlab = x, ylab = y, main = paste(y, "~", x), pch = 20)
  abline(intercept, slope, col = "blue", lwd = 2)
  if (boots) {
    plotboots(object, data)
  }
}
}


#' plot the regression line for each bootstrap
#' @param object a class newblblm object
#' @param data the original data to be plotted
plotboots <- function(object, data) {
  i <- 0
  j <- 0
  while (i < length(object$estimates)) {
    subsetest <- object$estimates %>% pluck(i)
    i <- i + 1
    while (j < length(subsetest)) {
      bootest <- subsetest %>% pluck(j)
      int <- bootest$coef %>% pluck(1)
      slo <- bootest$coef %>% pluck(2)
      abline(int, slo, col = alpha(rgb(0, 0, 0), 0.15), lwd = 0.75)
      j <- j + 1
    }
  }
}



#' Generic function coefplot
#' @param object the model
#' @param violin logical. Option to show a violin plot.
#' @param boots logical. Option to show a scatter plot.
#' @param box logical. Option to show a boxplot.
#'
#' @export
#'
coefplot <- function(object, violin = FALSE, boots = FALSE, box = FALSE) {
  UseMethod("coefplot")
}

#' Plot Boostrap Coefficient Estimates
#' @description Provides a variety of plot to visualize the bootstrap estimates.
#' @param object a newblblm class object
#' @param violin logical. Option to show a violin plot of the distribution of the bootstrap estimates.
#' @param boots logical. Option to show the estimation from each bootstrap.
#' @param box logical. Option to show a boxplot for the distribution of the bootstrap estimates.
#'
#' @export
#' @method coefplot newblblm
#'
coefplot.newblblm <- function(object, violin = FALSE, boots = FALSE, box = FALSE) {
  i <- 1
  x <- NULL
  while (i <= length(object$estimates)) {
    subset <- object$estimates %>% pluck(i)
    i <- i + 1
    j <- 1
    while (j <= length(subset)) {
      boot <- subset %>%
        pluck(j) %>%
        as.data.frame()
      bootest <- as.numeric(boot$coef)[-1]
      x <- c(x, bootest)
      j <- j + 1
    }
  }
  predictors <- attr(terms(object$formula), "term.labels")
  totalboots <- (i - 1) * (j - 1)
  y <- rep(predictors, totalboots)
  dat <- data.frame(x, y)
  x1 <- as.numeric(coef(object)[-1])
  y1 <- predictors
  p <- ggplot(subset(dat, y %in% predictors), aes(x = x, y = y, colour = y)) +
    facet_wrap(~y, ncol = 1, scales = "free")
  if (violin) {
    p <- p + geom_violin(trim = FALSE)
  }
  if (boots) {
    p <- p + geom_jitter(width = 0, alpha = 0.4)
  }
  if (box) {
    p <- p + geom_boxplot(alpha = 0, colour = "black")
  }
  p
}



