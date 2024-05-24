# https://chrissl789.github.io/modern-regression/notes/Lec03.SimpleLinearRegression.html

library(knitr)
library(pander)

# updated htest supports CI
pander.htest <- function (x, caption = attr(x, "caption"), ...) {
  if (is.null(caption)) {
    if (is.null(pander:::storage$caption)) {
      caption <- paste0(x$method, ": `", gsub("( and | by )", 
                                              "`\\1`", paste(x$data.name, collapse = "")), 
                        "`")
    }
    else {
      caption <- get.caption()
    }
  }
  res <- data.frame(placeholder = "FOO")
  if (!is.null(x$statistic)) {
    res$"Test statistic" <- as.numeric(x$statistic)
  }
  if (!is.null(x$parameter)) {
    res[names(x$parameter)] <- x$parameter
  }
  if (!is.null(x$p.value)) {
    res$"P value" <- paste(format(round(x$p.value, panderOptions("round")), 
                                  trim = TRUE, digits = panderOptions("digits"), decimal.mark = panderOptions("decimal.mark")), 
                           add.significance.stars(x$p.value))
  }
  if(!is.null(x$conf.int)) {
    cilevel <- sprintf('%s%% CI', round(100 * attr(x$conf.int, 'conf.level')))
    civalue <- format(x$conf.int[1:2], digits = panderOptions("digits"))
    res[cilevel] <- sprintf('[%s, %s]', civalue[1], civalue[2])
  }
  if (!is.null(x$alternative)) {
    res["Alternative hypothesis"] <- x$alternative
  }
  if (!is.null(x$estimate)) {
    if (!is.null(names(x$estimate))) {
      res[names(x$estimate)] <- x$estimate
    }
    else {
      res["Estimate"] <- x$estimate
    }
  }
  res$placeholder <- NULL
  pandoc.table(res, caption = caption, ...)
}

pander.summary.stanreg <- function(x, ...) {
  pd <- panderOptions('digits')
  on.exit(panderOptions('digits', pd))
  atts <- attributes(x)
  digits <- max(pd, atts$print.digits)
  panderOptions('digits', digits)
  res <- data.frame(placeholder = "FOO")
  res[,'function'] <- atts$stan_function
  res[,'family'] <- atts$family
  res[,'formula'] <- rstanarm:::formula_string(atts$formula)
  res[,'algorithm'] <- atts$algorithm
  if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling") {
    res[,'sample'] <- atts$posterior_sample_size
  }
  res[,'observations'] <- atts$nobs
  if (!is.null(atts$npreds)) {
    res[,'predictors'] <- atts$npreds
  }
  if (!is.null(atts$call$subset)) {
    res[,'subset'] <- deparse(atts$call$subset)
  }
  if (!is.null(atts$ngrps)) {
    res[,'groups'] <- sprintf('%s (%s)', names(atts$ngrps), unname(atts$ngrps))
  }
  res$placeholder <- NULL
  pandoc.table(res, caption = 'Model Info', ...)

  if (rstanarm:::used.optimizing(atts) || rstanarm:::used.variational(atts)) {
    hat <- "khat"
    str_diag <- "Monte Carlo diagnostics"
    str1 <- "and khat is the Pareto k diagnostic for importance sampling"
    str2 <- " (perfomance is usually good when khat < 0.7).\n"
  } else {
    hat <- "Rhat"
    str_diag <- "MCMC diagnostics"
    str1 <- "and Rhat is the potential scale reduction factor on split chains"
    str2 <- " (at convergence Rhat=1).\n"
  }
  sel <- which(colnames(x) %in% c("mcse", "n_eff", hat))
  has_mc_diagnostic <- length(sel) > 0
  if (has_mc_diagnostic) {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    mcse_hat <- format(round(x[, c("mcse", hat), drop = FALSE], digits), nsmall = digits)
    n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
    mcdx <- cbind(mcse_hat, n_eff)
  }
  else {
    xtemp <- x
  }
  ppd_nms <- grep("^mean_PPD", rownames(x), value = TRUE)
  has_ppd_diagnostic <- !atts$no_ppd_diagnostic && length(ppd_nms) > 0
  if (has_ppd_diagnostic) {
    ppd_estimates <- xtemp[rownames(xtemp) %in% ppd_nms, , drop = FALSE]
  } else {
    ppd_estimates <- NULL
  }
  xtemp <- xtemp[!rownames(xtemp) %in% c(ppd_nms, "log-posterior"), , drop = FALSE]
  pandoc.table(xtemp, caption = 'Estimates', ...)

  if (has_ppd_diagnostic) {
    pandoc.table(ppd_estimates, caption = 'Fit Diagnostics', ...)
    cat("\nThe mean_ppd is the sample average posterior predictive distribution of the outcome variable.\n")
  }

  if (has_mc_diagnostic) {
    pandoc.table(mcdx, caption = str_diag, ...)
    dx_txt <- sprintf(
'\nFor each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, %s %s\n',
str1, str2
)
    cat(dx_txt)
  }
}

pander.prior_summary.stanreg <- function(x, ...) {
  pd <- panderOptions('digits')
  on.exit(panderOptions('digits', pd))
  atts <- attributes(x)
  .dig <- max(pd, attr(x, "print_digits"))
  panderOptions('digits', .dig)
  .fr2 <- function(y, .digits = .dig, ...) format(y, digits = .digits, ...)
  .fr3 <- function(y, .nsmall = .dig) .fr2(y, nsmall = .nsmall)

  QR <- attr(x, "QR")
  sparse <- attr(x, "sparse")
  model_name <- attr(x, "model_name")
  msg <- sprintf("Priors for model '%s'", model_name)
  stan_function <- attr(x, "stan_function")
  if (stan_function == "stan_mvmer") {
    return(rstanarm:::print.prior_summary.stanreg(x, ...))
  }

  pi_txt <- paste0('Intercept', if(!sparse) " (after predictors centered)")
  p_txt <- paste0('Coefficients', if(QR) " (in Q-space)")
  pa_txt <- ''
  rdf <- data.frame(matrix(NA, nrow = 3, ncol = 2))
  pfun <- function(p) {
    d <- p$dist
    a <- p$adjusted_scale
    if(d == 'exponential') {
      sprintf('%s(rate = %s)', d, .fr3(c(p$rate, 1 / a)))
    } else {
      sprintf('%s(location = %s, scale = %s)', d, .fr2(p$location), .fr3(c(p$scale, a)))
    }
  }
  if (!is.null(x[["prior_intercept"]])) {
    rdf[1,] <- pfun(x[["prior_intercept"]])
  }
  if (!is.null(x[["prior"]])) {
    rdf[2,] <- pfun(x[["prior"]])
  }
  if (!is.null(x[["prior_aux"]])) {
    p <- x[["prior_aux"]]
    aux_dist <- p$dist
    if (aux_dist %in% c("normal", "student_t", "cauchy")) {
      p[['dist']] <- paste0('half-', aux_dist)
    }
    rdf[3,] <- pfun(p)
    pa_txt <- sprintf('Auxiliary (%s)', p$aux_name)
  }
  names(rdf) <- c('Specified prior', 'Adjusted prior')
  rownames(rdf) <- c(pi_txt, p_txt, pa_txt)
  # remove rows with no data
  rdf <- rdf[rowSums(!is.na(rdf)) > 0,]
  pandoc.table(rdf, caption = msg, ...)
}

# knit_print.data.frame = function(x, ...) {
#   # res = paste(c("", "", knitr::kable(x)), collapse = "\n")
#   # knitr::asis_output(res)
#   pander::pander(x, ...)
# }
# registerS3method("knit_print", "data.frame", knit_print.data.frame, envir = asNamespace("knitr"))
gen_pander <- function(x, ...) pander::pander(x, ..., split.table = 120, split.cells = 20)

preglist <- c(
  'data.frame',
  'ols',
  'finalfit', # consider split.cells and split.table
  'summary.lm', # missing residuals and fstatistic
  'matrix',
  'htest', # updated function
  'summary.stanreg', # new function
  'prior_summary.stanreg' # new function - very limited
)

for(i in seq_along(preglist)) {
  registerS3method("knit_print", preglist[i], gen_pander, envir = asNamespace("knitr"))
}
