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

pander.summary.survreg <- function(x, summary = TRUE,
                                   digits = panderOptions("digits"),
                                   round = panderOptions("round"),
                                   keep.trailing.zeros = panderOptions("keep.trailing.zeros"),
                                   ...) {
  if (!is.null(cl <- x$call)) {
    cat("\nCall:", pandoc.formula.return(cl), "", sep = "\n\n")
  }
  if (!is.null(x$fail)) {
    cat(" Survreg failed.", x$fail, "\n\n")
    return(invisible())
  }
  if (summary) {
    pandoc.table(x$table, caption = "Model statistics", digits = digits, 
                 round = round, keep.trailing.zeros = keep.trailing.zeros, 
                 ...)
  }
  else {
    coef <- x$coef
    if (any(nas <- is.na(coef))) {
      if (is.null(names(coef))) {
        names(coef) <- paste("b", 1:length(coef), sep = "")
      }
      cat("\nCoefficients: (", sum(nas), " not defined because of singularities)\n", 
          sep = "")
    }
    pandoc.table(coef, caption = "Coefficients", digits = digits, 
                 round = round, keep.trailing.zeros = keep.trailing.zeros, 
                 ...)
  }
  if (nrow(x$var) == length(coef)) {
    cat("\nScale fixed at", format(x$scale), "\n")
  }
  else if (length(x$scale) == 1) {
    cat("\nScale=", format(x$scale), "\n")
  }
  else {
    pandoc.table(x$scale, caption = "Scale", ...)
  }
  nobs <- length(x$linear)
  if(is.null(x$linear)) nobs <- x$n
  chi <- 2 * diff(x$loglik)
  df <- sum(x$df) - x$idf
  pandoc.table(data.frame(`Loglik(model)` = x$loglik[2], `Loglik(intercept only)` = x$loglik[1]), ...)
  if (df > 0) {
    cat("Chisq=", p(chi, wrap = ""), "on", p(df, wrap = ""), 
        "degrees of freedom, p=", p(signif(1 - pchisq(chi, 
                                                      df), 2), wrap = ""), "\n\n")
  }
  else {
    cat("\n")
  }
  if (summary) {
    if (x$robust) {
      cat("(Loglikelihood assumes independent observations)\n\n")
    }
    cat("Number of Newton-Raphson Iterations:", p(trunc(x$iter), 
                                                  wrap = ), "\n\n")
  }
  omit <- x$na.action
  if (length(omit)) {
    cat("n=", nobs, " (", naprint(omit), ")\n", sep = "")
  }
  else {
    cat("n=", nobs, "\n")
  }
  if (summary) {
    if (!is.null(correl <- x$correlation)) {
      p <- dim(correl)[2]
      if (p > 1) {
        ll <- lower.tri(correl)
        correl <- apply(correl, c(1, 2), p, wrap = "", 
                        digits = digits, round = round, keep.trailing.zeros = keep.trailing.zeros)
        correl[!ll] <- ""
        pander(correl[-1L, -ncol(correl)], digits = digits, 
               round = round, keep.trailing.zeros = keep.trailing.zeros, 
               caption = "Correlation of Coefficients", ...)
      }
    }
  }
  invisible()
}

pander.summary.survfit <- function (x, digits = panderOptions("digits"), ...) {
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  if (!is.null(cl <- x$call)) {
    cat("\nCall:", pandoc.formula.return(cl), "", sep = "\n\n")
  }
  omit <- x$na.action
  if (length(omit)) 
    cat(naprint(omit), "\n")
  if (x$type == "right" || is.null(x$n.enter)) {
    mat <- cbind(x$time, x$n.risk, x$n.event, x$surv)
    cnames <- c("time", "n.risk", "n.event")
  }
  else if (x$type == "counting") {
    mat <- cbind(x$time, x$n.risk, x$n.event, x$n.censor, 
                 x$surv)
    cnames <- c("time", "n.risk", "n.event", "censored")
  }
  if (is.matrix(x$surv)) 
    ncurve <- ncol(x$surv)
  else ncurve <- 1
  if (ncurve == 1) {
    cnames <- c(cnames, "survival")
    if (!is.null(x$std.err)) {
      if (is.null(x$lower)) {
        mat <- cbind(mat, x$std.err)
        cnames <- c(cnames, "std.err")
      }
      else {
        mat <- cbind(mat, x$std.err, x$lower, x$upper)
        cipct <- x$conf.int * 100
        cnames <- c(cnames, "std.err", sprintf('lower %s%% CI', cipct), sprintf('upper %s%% CI', cipct))
      }
    }
  }
  else cnames <- c(cnames, paste("survival", seq(ncurve), sep = ""))
  if (!is.null(x$start.time)) {
    mat.keep <- mat[, 1] >= x$start.time
    mat <- mat[mat.keep, , drop = FALSE]
    if (is.null(dim(mat))) 
      stop(paste("No information available using start.time =", x$start.time, "."))
  }
  if (!is.matrix(mat)) 
    mat <- matrix(mat, nrow = 1)
  if (!is.null(mat)) {
    dimnames(mat) <- list(rep("", nrow(mat)), cnames)
    if (is.null(x$strata)) 
      pandoc.table(mat, ...)
    else {
      strata <- x$strata
      if (!is.null(x$start.time)) 
        strata <- strata[mat.keep]
      for (i in levels(strata)) {
        who <- (strata == i)
        cat("               ", i, "\n")
        pandoc.table(mat[who, ], ...)
        cat("\n")
      }
    }
  }
  else stop("There are no events to print.  Please use the option ", 
            "censored=TRUE with the summary function to see the censored ", 
            "observations.")
  invisible(x)
}

pander.summary.coxph <- function (x, digits = panderOptions("digits"), ...) {
  signif.stars = getOption("show.signif.stars")
  expand = FALSE
  if (!is.null(cl <- x$call)) {
    cat("\nCall:", pandoc.formula.return(cl), "", sep = "\n\n")
  }
  if (!is.null(x$fail)) {
    cat(" Coxreg failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  omit <- x$na.action
  m <- matrix(c(x$n), dimnames = list('n', NULL))
  if (!is.null(x$nevent)) 
    m <- rbind(m, 'number of events' = x$nevent)
  pandoc.table(t(m))
  if (length(omit)) 
    cat("   (", naprint(omit), ")\n", sep = "")
  if (nrow(x$coef) == 0) {
    cat("   Null model\n")
    return()
  }
  if (expand && !is.null(x$cmap)) {
    signif.stars <- FALSE
    tmap <- x$cmap
    cname <- colnames(tmap)
    printed <- rep(FALSE, length(cname))
    for (i in 1:length(cname)) {
      if (!printed[i]) {
        j <- apply(tmap, 2, function(x) all(x == tmap[, 
                                                      i]))
        printed[j] <- TRUE
        tmp2 <- x$coefficients[tmap[, i], , drop = FALSE]
        names(dimnames(tmp2)) <- c(paste(cname[j], collapse = ", "), 
                                   "")
        rownames(tmp2) <- rownames(tmap)[tmap[, i] > 
                                           0]
        printCoefmat(tmp2, digits = digits, P.values = TRUE, 
                     has.Pvalue = TRUE, signif.legend = FALSE, signif.stars = signif.stars, 
                     ...)
        if (!is.null(x$conf.int)) {
          tmp2 <- x$conf.int[tmap[, i], , drop = FALSE]
          rownames(tmp2) <- rownames(tmap)[tmap[, i] > 
                                             0]
          names(dimnames(tmp2)) <- c(paste(cname[j], 
                                           collapse = ", "), "")
          print(tmp2, digits = digits, ...)
        }
      }
    }
    cat("\n States:", paste(paste(seq(along.with = x$states), 
                                  x$states, sep = "= "), collapse = ", "), "\n")
  } else {
    if (!is.null(x$coefficients)) {
      pandoc.table(x$coefficients, ...)
    }
    if (!is.null(x$conf.int)) {
      pandoc.table(x$conf.int, ...)
    }
  }
  if (!is.null(x$concordance)) {
    v_conc <- format(round(x$concordance, digits))
    vdf <- data.frame(Concordance = v_conc[1], se = v_conc[2], row.names = NULL)
    pander(vdf, ...)
  }
  pdig <- max(1, digits - 4)
  v_test <- c(x$logtest['test'], x$waldtest['test'], x$sctest['test'])
  v_df <- c(x$logtest['df'], x$waldtest['df'], x$sctest['df'])
  v_pval <- c(x$logtest['pvalue'], x$waldtest['pvalue'], x$sctest['pvalue'])
  f_rn <- c('Likelihood ratio test', 'Wald test', 'Score (logrank) test')
  if (!is.null(x$robscore)) {
    v_test <- c(v_test, x$robscore['test'])
    v_df <- c(v_df, x$robscore['df'])
    v_pval <- c(v_pval, x$sctest['pvalue'])
    f_rn <- c(f_rn, 'Robust score test')
  }
  f_test <- format(round(v_test, 2))
  f_df <- v_df
  f_pval <- format.pval(v_pval, digits = pdig)
  rdf <- data.frame(test = f_test, df = f_df, p = f_pval, row.names = f_rn)
  pander(rdf, ...)
  if (x$used.robust) 
    cat("  (Note: the likelihood ratio and score tests", 
        "assume independence of\n     observations within a cluster,", 
        "the Wald and robust score tests do not).\n")
  invisible()
}

pander.aftreg <- function (x, digits = panderOptions("digits"), ...) {
  if (!is.null(cl <- x$call)) {
    cat("\nCall:", pandoc.formula.return(cl), "", sep = "\n\n")
  }
  if (!is.null(x$fail)) {
    cat(" aftreg failed.\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  if (x$pfixed) {
    n.slsh <- 1
  } else {
    n.slsh <- 2 * x$n.strata
  }
  coef <- x$coefficients
  se <- sqrt(diag(x$var))
  wald.p <- 1 - pchisq((coef/se)^2, 1)
  if (is.null(coef) || is.null(se)) 
    stop("Input is not valid")
  cnames <- c('Covariate', 'W.mean', 'Coef', 'Time-Expn', 'se(Coef)', 'Wald p')
  if (x$param == "lifeAcc") {
    cnames[4] <- 'Time-Accn'
  }
  e.coef <- exp(coef)
  ett <- formatC(1, width = 9, digits = 0, format = "f")
  noll <- formatC(0, width = 5, digits = 0, format = "f")
  factors <- attr(x$terms, "factors")
  resp <- attr(x$terms, "response")
  row.strata <- attr(x$terms, "specials")$strata
  if (!is.null(row.strata)) 
    col.strata <- which(factors[row.strata, ] == 1)
  else col.strata <- NULL
  if (!is.null(x$covars)) {
    if (!is.null(col.strata)) {
      factors <- attr(x$terms, "factors")[-c(resp, row.strata), 
                                          -col.strata, drop = FALSE]
    }
    else {
      factors <- attr(x$terms, "factors")[-c(resp, row.strata), 
                                          , drop = FALSE]
    }
    covar.names <- c(x$covars, names(coef)[(length(coef) - n.slsh + 1):length(coef)])
    term.names <- colnames(factors)
    isF <- x$isF
  }
  ord <- attr(x$terms, "order")
  if (!is.null(col.strata)) 
    ord <- ord[-col.strata]
  index <- 0
  xdf <- as.data.frame(matrix(NA, 10, 6))
  bdf <- xdf
  x_new <- xdf[1,]
  if (!is.null(x$covars)) {
    n.rows <- length(term.names)
    for (term.no in 1:n.rows) {
      if (ord[term.no] == 1) {
        covar.no <- which(factors[, term.no] == 1)
        if (isF[covar.no]) {
          cat(covar.names[covar.no], "\n")
          no.lev <- length(x$levels[[covar.no]])
          x_new[1,1] <- x$levels[[covar.no]][1]
          x_new[1,2] <- x$w.means[[covar.no]][1]
          x_new[1,3] <- noll
          x_new[1,4] <- ett
          x_new[1,5] <- '(reference)'
          for (lev in 2:no.lev) {
            index <- index + 1
            xdf[index,1] <- x$levels[[covar.no]][lev]
            xdf[index,2] <- x$w.means[[covar.no]][lev]
            xdf[index,3] <- coef[index]
            xdf[index,4] <- e.coef[index]
            xdf[index,5] <- se[index]
            xdf[index,6] <- wald.p[index]
          }
        } else {
          index <- index + 1
          xdf[index,1] <- covar.names[covar.no]
          xdf[index,2] <- x$w.means[[covar.no]]
          xdf[index,3] <- coef[index]
          xdf[index,4] <- e.coef[index]
          xdf[index,5] <- se[index]
          xdf[index,6] <- wald.p[index]
        }
      } else if (ord[term.no] > 1) {
        cat(format(term.names[term.no], width = 16), "\n")
        niv <- numeric(ord[term.no])
        covar.no <- which(factors[, term.no] == 1)
        for (i in 1:ord[term.no]) {
          if (isF[covar.no[i]]) {
            niv[i] <- length(x$levels[[covar.no[i]]]) - 1
          } else {
            niv[i] <- 1
          }
        }
        stt <- index + 1
        for (index in stt:(stt + prod(niv) - 1)) {
          vn <- sub(covar.names[covar.no[1]], "", names(coef)[index])
          for (i in 1:ord[term.no]) {
            vn <- sub(covar.names[covar.no[i]], "", vn)
          }
          xdf[index,1] <- ''
          xdf[index,2] <- substring(vn, 1, 22)
          xdf[index,3] <- coef[index]
          xdf[index,4] <- e.coef[index]
          xdf[index,5] <- se[index]
          xdf[index,6] <- wald.p[index]
        }
      }
    }
    if(rowSums(!is.na(x_new)) > 0) xdf <- rbind(x_new, xdf)
    names(xdf) <- cnames
    xdf <- xdf[rowSums(!is.na(xdf)) > 0,]
    pandoc.table(xdf, ...)
  }
  base_index <- index
  for (i in 1:n.slsh) {
    jup <- length(coef)
    ss.names <- names(coef[(jup - n.slsh + 1):jup])
    index <- index + 1
    bdf[index-base_index, 1] <- ss.names[i]
    bdf[index-base_index, 3] <- coef[index]
    bdf[index-base_index, 5] <- se[index]
    bdf[index-base_index, 6] <- wald.p[index]
  }
  bdf <- bdf[rowSums(!is.na(bdf)) > 0,]
  names(bdf) <- cnames
  pandoc.table(bdf, caption = 'Baseline parameters', ...)
  cat("Baseline life expectancy: ", x$baselineMean)
  cat('\n\n')
  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  if (is.null(x$df)) 
    df <- sum(!is.na(coef)) - n.slsh
  else df <- round(sum(x$df), 2)
  ldf <- data.frame()
  if (x$pfixed) {
    cat(" Shape is fixed at ", x$shape, "\n\n")
  }
  str1 <- paste(formatC("Events", width = 25, flag = "-"), x$n.events)
  str2 <- paste(formatC("Total time at risk", width = 25, flag = "-"),
                formatC(x$ttr, digits = 5, format = "fg"))
  str3 <- paste(formatC("Max. log. likelihood", width = 25, flag = "-"),
                formatC(x$loglik[2], digits = 5, format = "fg"))
  bigstr <- c(str1, str2, str3)
  if (df > 0.5) {
    val_lr <- format(round(logtest, 2))
    val_df <- formatC(df, digits = 0, format = "f")
    val_p <- format.pval(1 - pchisq(logtest, df), digits = 6)
    str4 <- paste(formatC("LR test statistic", width = 25, flag = "-"), val_lr)
    str5 <- paste(formatC("Degrees of freedom", width = 25, flag = "-"), val_df)
    str6 <- paste(formatC("Overall p-value", width = 25, flag = "-"), val_p)
    bigstr <- c(bigstr, str4, str5, str6)
  }
  pander(paste(bigstr, collapse = '\n\n'), ...)
  if (length(x$icc)) {
    cat("\n")
    cat("   number of clusters=", x$icc[1], "    ICC=", format(x$icc[2:3]), "\n")
  }
  invisible(x)
}

pander.coeftest <- function (x, digits = panderOptions("digits"), ...) {
  if (is.null(d <- dim(x)) || length(d) != 2L) 
    stop("'x' must be coefficient matrix/data frame")
  xm <- data.matrix(x)
  xmethod <- attr(xm, 'method')
  attr(xm, 'method') <- NULL
  attr(xm, 'df') <- NULL
  attr(xm, 'nobs') <- NULL
  attr(xm, 'logLik') <- NULL
  class(xm) <- 'numeric'
  pandoc.table(xm, caption = xmethod, ...)
  invisible(x)
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
  'lrtest',
  'coxph',
  'survdiff',
  'survfit',
  'summary.coxph', # new function
  'summary.survfit', # new function
  'summary.survreg', # updated function
  'aftreg', # - package "eha" - new function/fragile
  'coeftest', # package "lmtest" - new function
  'lm',
  'grouped_df' # "tibble"
  # epi.2by2 - in epiR, new method
)

# need consistent pvalue format

for(i in seq_along(preglist)) {
  registerS3method("knit_print", preglist[i], gen_pander, envir = asNamespace("knitr"))
}
