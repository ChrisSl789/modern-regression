# https://chrissl789.github.io/modern-regression/notes/Lec03.SimpleLinearRegression.html

library(knitr)
library(pander)

panderOptions("digits", 5)
panderOptions("round", 10)
panderOptions('keep.trailing.zeros', TRUE)
panderOptions('missing', '')
panderOptions("table.split.table", 120)
panderOptions("table.split.cells", 20)

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
  # `see rstanarm:::.print_scalar_prior`
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
      out <- sprintf('~ %s(rate = %s)', d, .fr3(c(p$rate, 1 / a)))
    } else {
      out <- sprintf('~ %s(location = %s, scale = %s)', d, .fr2(p$location), .fr3(c(p$scale, a)))
    }
    if(length(out) == 1) out <- c(out, '')
    ### CHECK
    # when length=2, labels are specified/adjusted
    # when length=1, use specified/NA or other?
    ###
    out
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

pander.summary.glm <- function(x, ...) {
  pd <- max(4L, panderOptions('digits')+1)
  pander:::pander.summary.lm(x, ...)
  if (!is.null(x$na.action) && nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n\n", sep = "")
  }
  if(!is.null(x$aic)) {
    cat(sprintf('AIC: %s\n\n', format(x$aic, digits = pd)))
  }
  if(!is.null(x$iter)) {
    cat(sprintf('Number of Fisher Scoring Iterations: %s\n', x$iter))
  }
}

pander.glm <- function(x, ...) {
  pander:::pander.glm(x, ...)
  pd <- max(4L, panderOptions('digits')+1)
  res <- matrix(NA, 5, 1)
  rownames(res) <- c('Degrees of Freedom', 'Degrees of Freedom', 'Null Deviance', 'Residual Deviance', 'AIC')
  if(!is.null(x$df.null)) {
    res[1,1] <- sprintf('%s Total (i.e. Null)', x$df.null)
  }
  if(!is.null(x$df.residual)) {
    res[2,1] <- sprintf('%s Residual', x$df.residual)
  }
  if(!is.null(x$null.deviance)) {
    res[3,1] <- format(x$null.deviance, digits = pd)
  }
  if(!is.null(x$deviance)) {
    res[4,1] <- format(x$deviance, digits = pd)
  }
  if(!is.null(x$aic)) {
    res[5,1] <- format(x$aic, digits = pd)
  }
  res <- res[!is.na(res[,1]),,drop = FALSE]
  pandoc.table(res, ...)
  if (!is.null(x$na.action) && nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n\n", sep = "")
  }
}

pander.data.frame.ff <- function(x, ...) {
  rownames(x) <- NULL
  pander:::pander.data.frame(x, ...)
}

pander.summary.lm <- function(x, ...) {
  title <- pandoc.formula.return(x$call$formula, text = "Fitting linear model:")
  cat(sprintf('\n\n%s\n\n', title))
  pd <- max(4L, panderOptions('digits')+1)
  resid <- x$residuals
  rdf <- x$df[2L]
  res_cap <- 'Residuals'
  if(!is.null(x$weights) && diff(range(x$weights))) res_cap <- paste('Weighted', res_cap)
  if (rdf > 5L) {
    stopifnot(length(dim(resid)) == 0)
    rq <- zapsmall(quantile(resid), pd + 1)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    pandoc.table(rq, caption = res_cap, ...)
  } else if (rdf > 0L) {
    pandoc.table(resid, ...)
  }
  pander:::pander.summary.lm(x, caption = '', ...)

  if (!is.null(x$fstatistic)) {
    fstat <- formatC(x$fstatistic[1L], digits = pd)
    pval <- format.pval(pf(x$fstatistic[1L], x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), digits = pd)
    cat(sprintf('F-statistic: %s on %s and %s DF,  pvalue: %s\n', fstat, x$fstatistic[2L], x$fstatistic[3L], pval))
  }
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
}

# knit_print.data.frame = function(x, ...) {
#   # res = paste(c("", "", knitr::kable(x)), collapse = "\n")
#   # knitr::asis_output(res)
#   pander::pander(x, ...)
# }
# registerS3method("knit_print", "data.frame", knit_print.data.frame, envir = asNamespace("knitr"))
gen_pander <- function(x, ...) pander::pander(x, ...)

preglist <- c(
  'data.frame',
  'ols',
  'data.frame.ff', # finalfit()
  'summary.lm', # missing residuals and fstatistic
  'matrix',
  'htest', # updated function
  'summary.stanreg', # new function
  'prior_summary.stanreg', # new function - very limited
  'glm',
  'summary.glm', # updated function
  'summary.lm', # updated function
  'table',
  'lrm',
  'summary.lrm',
  'robcov', # WANT: format pvalue, formatNP?
  'summary.robcov',
  'lrtest'
  # epi.2by2 - in epiR, new method
)

# highlight STATA?

for(i in seq_along(preglist)) {
  registerS3method("knit_print", preglist[i], gen_pander, envir = asNamespace("knitr"))
}
