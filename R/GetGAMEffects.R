#' Get GAM Effects
#'
#' @description Use this version if you want the confidence intervals to be based on the cv runs.
# This function has to be very complicated to account for possibility of missing data in various places
# note, needs additional work to fit gaussian models
#' @param model a model fit using the "gam" function from mgcv
#' @param data data frame; data used to fit the models and cv models
#' @param cv.model.list list;cv models to be used to determine confidence bounds, usually from CrossValidateModel function
#' @param vars character; a vector of model terms to be plotted, or "all" to return all terms
#' @param scale character; Either "log" or "abundance", what should the outputs be in
#' @param scale.factor a factor to multiply results
#'
#' @return a list of data frames containing effect estimates and optionally confidence intervals for each term
#' @export
#'
#' @examples
GetGAMEffects <- function(model,
                          data,
                          cv.model.list = NULL,
                          vars = "all",
                          scale = "log",
                          scale.factor = 1) {

  # detect and remove any models that are NULL or NA, and the accompanying types
  if (is.null(cv.model.list) == F) {
    list.index <- 1
    cv.model.list1 <- list()
    for (m in 1:length(cv.model.list)) {
      if (is.list(cv.model.list[[m]])) {
        cv.model.list1[[list.index]] <- cv.model.list[[m]]
        list.index <- list.index + 1
      }
    }
  }
  # get a list of the actual terms in the model;Need a check in case the p model has variables the main one doesn't
  var.table <- AutodetectGAMTerms(model)

  if (model$family$family == "ziplss") {
    p.only <- var.table[[2]]$term[var.table[[2]]$term %in% var.table[[1]]$term == F]
    d.only <- var.table[[1]]$term[var.table[[1]]$term %in% var.table[[2]]$term == F]
    var.table <- unique(rbind(var.table[[1]], var.table[[2]]))
  } else {
    p.only <- NA
  }

  average.vars <- stats::na.omit(c(var.table$term[var.table$type %in% c("factor", "offset") == F], var.table$term2))
  fac.vars <- var.table$term[var.table$type == "factor"]
  off.var <- var.table$term[var.table$type == "offset"]

  if (tolower(vars) == "all") {
    do.rows <- which(var.table$type != "offset")
  } else {
    var.match <- unlist(strsplit(vars, split = "[*]"))
    do.rows <- vector(length = length(var.match))
    for (i in 1:length(var.match)) {
      do.rows[i] <- which(var.table$term == var.match[i] | var.table$term2 == var.match[i])
    }
    do.rows <- unique(do.rows)
  }

  # take the averages to set things up for making the predictions
  average.dat <- apply(data.frame(data[, average.vars]), FUN = "mean", MARGIN = 2)
  fac.dat <- rep(0, times = length(fac.vars))
  off.dat <- apply(data.frame(data[, off.var]), FUN = "mean", MARGIN = 2)

  average.dat <- c(average.dat, fac.dat, off.dat)
  average.dat <- data.frame(t(average.dat))
  colnames(average.dat) <- c(average.vars, fac.vars, off.var)

  # fix the factors
  if (length(fac.vars) > 0) {
    for (f in 1:length(fac.vars)) {
      average.dat[, fac.vars[f]] <- factor(average.dat[, fac.vars[f]], levels = levels(as.factor(data[, fac.vars[f]])))
    }
  }

  effects.list <- list()
  list.index <- 1
  var.table1 <- var.table[do.rows, ]

  # get down to business
  vars2d <- var.table1[var.table1$dims == 2, ]
  vars1d <- var.table1[var.table1$type == "smooth" & var.table1$dims == 1, ]
  varsf <- var.table1[var.table1$type == "factor", ]

  if (nrow(vars2d) > 0) {
    for (v in 1:nrow(vars2d)) {
      term2d <- c(vars2d$term[v], vars2d$term2[v])
      x.seq <- seq(from = min(data[, term2d[1]]), to = max(data[, term2d[1]]), length.out = 40)
      y.seq <- seq(from = min(data[, term2d[2]]), to = max(data[, term2d[2]]), length.out = 40)

      v2.data <- data.frame(
        x = rep(x.seq, times = 40), y = rep(y.seq, each = 40),
        average.dat[, names(average.dat) %in% term2d == F]
      )
      names(v2.data)[1:2] <- term2d

      v2.preds <- as.numeric(mgcv::predict.gam(model, type = "terms", newdata = v2.data, terms = paste0("s(", term2d[1], ",", term2d[2], ")"))[, 1])

      # special handling for ziplss models
      if (model$family$family == "ziplss") {
        v2.probs <- as.numeric(mgcv::predict.gam(model, type = "terms", newdata = v2.data, terms = paste0("s.1(", term2d[1], ",", term2d[2], ")"))[, 1])
        v2.preds <- exp(v2.preds) * (1 - exp(-exp(v2.probs))) / (1 - stats::dpois(0, exp(v2.preds)))
        v2.preds <- log(v2.preds)
      }

      # adjust the scale as necessary; in log scale, use a small number instead of zero
      if (scale == "abund") {
        v2.preds <- exp(v2.preds) * scale.factor
      } else {
        v2.preds <- v2.preds + log(scale.factor)
      }
      # we are going to piggy back off the mgcv functions to figure out which ones are NAs
      # that way the plot won't extrapolate out of the sample area
      # best way to do this is via the built in plot function, but need a dummy png so that it doesn't output
      grDevices::png("trashme.png")
      x <- plot(model, scale = 0, se = F, pages = 1)
      grDevices::dev.off()
      file.remove("trashme.png")

      xvec <- vector(length = length(x))

      for (n in 1:length(x)) {
        xvec[n] <- x[[n]]$xlab == term2d[1]
      }
      v2nas <- which(is.na(x[[which(xvec)[1]]]$fit))

      v2.preds[v2nas] <- NA
      v2.preds[is.infinite(v2.preds)] <- NA

      effects.list[[list.index]] <- data.frame(x.seq, y.seq, effect = v2.preds)
      names(effects.list[[list.index]])[1:2] <- c("x", "y")
      names(effects.list)[list.index] <- paste(term2d, collapse = "*")

      list.index <- list.index + 1
    }
  }

  # now start the smooths
  if (nrow(vars1d) > 0) {
    for (v in 1:nrow(vars1d)) {
      term1d <- vars1d$term[v]
      v.data <- data.frame(
        seq(from = min(data[, term1d]), to = max(data[term1d]), length.out = 100),
        average.dat[names(average.dat) != term1d]
      )
      names(v.data)[1] <- term1d

      v.name <- paste0("s(", term1d, ")")


      # first get the main effect
      main.pred <- as.numeric(mgcv::predict.gam(model, newdata = v.data, type = "terms", terms = v.name))
      se.pred <- as.numeric(mgcv::predict.gam(model, newdata = v.data, type = "terms", terms = v.name, se.fit = T)[[2]])
      # special handling for ziplss models
      if (term1d %in% p.only) {
        main.pred <- main.pred[1:nrow(v.data)]
        se.pred <- se.pred[1:nrow(v.data)]
      }
      if (model$family$family == "ziplss") {
        main.pred[-10 > main.pred] <- log(.01)
        main.prob <- as.numeric(mgcv::predict.gam(model, type = "terms", newdata = v.data, terms = paste0("s.1(", term1d, ")"))[, 1])
        main.pred <- exp(main.pred) * (1 - exp(-exp(main.prob))) / (1 - stats::dpois(0, exp(main.pred)))
        main.pred <- log(main.pred)
      }

      if (scale == "abund") {
        main.pred <- exp(main.pred) * scale.factor
      } else {
        main.pred <- main.pred + log(scale.factor)
      }
      main.pred[is.infinite(main.pred)] <- NA
      # now loop through to get the cv effects
      if (is.null(cv.model.list) == F) {
        effect.dat <- matrix(nrow = 100, ncol = length(cv.model.list))
        for (f in 1:length(cv.model.list)) {
          if (is.list(cv.model.list[[f]])) {
            cv.pred <- as.numeric(mgcv::predict.gam(cv.model.list[[f]],
                                                    newdata = v.data, type = "terms",
                                                    terms = v.name
            ))
            if (term1d %in% p.only) {
              cv.pred <- cv.pred[1:nrow(v.data)]
            }
            if (model$family$family == "ziplss") {
              cv.pred[-10 > cv.pred] <- log(.01)
              cv.prob <- as.numeric(mgcv::predict.gam(cv.model.list[[f]], type = "terms", newdata = v.data, terms = paste0("s.1(", term1d, ")"))[, 1])
              cv.pred <- exp(cv.pred) * (1 - exp(-exp(cv.prob))) / (1 - stats::dpois(0, exp(cv.pred)))
              cv.pred <- log(cv.pred)
            }
            if (scale == "abund") {
              cv.pred <- exp(cv.pred) * scale.factor
            } else {
              cv.pred <- cv.pred + log(scale.factor)
            }
          } else {
            cv.pred <- NA
          }
          effect.dat[, f] <- cv.pred
        }
        colnames(effect.dat) <- paste0("CV", 1:length(cv.model.list))

        effect.dat[is.infinite(effect.dat)] <- NA

        uppers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", probs = .95, na.rm = T)
        lowers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", probs = .05, na.rm = T)
        variance <- apply(X = effect.dat, MARGIN = 1, FUN = var, na.rm = T)

        out.dat <- data.frame(x = v.data[, 1], effect = main.pred, var = se.pred^2, cvvar = variance, upper = uppers, lower = lowers, effect.dat)
      } else {
        out.dat <- data.frame(x = v.data[, 1], effect = main.pred, var = se.pred^2)
      }

      effects.list[[list.index]] <- out.dat
      names(effects.list)[list.index] <- term1d
      list.index <- list.index + 1
    }
  }

  # Now do the factors
  if (nrow(varsf) > 0) {
    for (v in 1:nrow(varsf)) {
      termf <- varsf$term[v]

      # a little bit to make sure the names match
      if (termf %in% names(model$model)) {
        fac.name <- termf
      } else {
        fac.name <- paste0("as.factor(", termf, ")")
      }

      # set up the data
      f.data <- data.frame(
        unique(as.character(data[, termf])),
        average.dat[names(average.dat) != termf]
      )
      names(f.data)[1] <- termf

      # first get the main effect
      main.pred <- as.numeric(mgcv::predict.gam(model, newdata = f.data, type = "terms", terms = fac.name))
      se.pred <- as.numeric(mgcv::predict.gam(model, newdata = f.data, type = "terms", terms = fac.name, se.fit = T)[[2]])

      # special handling for ziplss models
      if (termf %in% p.only) {
        main.pred <- main.pred[1:nrow(f.data)]
        se.pred <- se.pred[1:nrow(f.data)]
      }
      if (model$family$family == "ziplss") {
        main.pred[-10 > main.pred] <- log(.01)
        main.prob <- as.numeric(mgcv::predict.gam(model, type = "terms", newdata = f.data, terms = paste0("as.factor(", termf, ").1"))[, 1])
        main.pred <- exp(main.pred) * (1 - exp(-exp(main.prob))) / (1 - stats::dpois(0, exp(main.pred)))
        main.pred <- log(main.pred)
      }

      if (scale == "abund") {
        main.pred <- exp(main.pred) * scale.factor
      } else {
        main.pred <- main.pred + log(scale.factor)
      }

      if (is.null(cv.model.list) == F) {
        effect.dat <- matrix(nrow = length(main.pred), ncol = length(cv.model.list))
        for (f in 1:length(cv.model.list)) {
          if (is.list(cv.model.list[[f]])) {
            cv.pred <- as.numeric(mgcv::predict.gam(cv.model.list[[f]],
                                                    newdata = f.data, type = "terms",
                                                    terms = fac.name
            ))
            if (termf %in% p.only) {
              cv.pred <- cv.pred[1:nrow(f.data)]
            }
            if (model$family$family == "ziplss") {
              cv.pred[-10 > cv.pred] <- log(.01)
              cv.prob <- as.numeric(mgcv::predict.gam(cv.model.list[[f]], type = "terms", newdata = f.data, terms = paste0("as.factor(", termf, ").1"))[, 1])
              cv.pred <- exp(cv.pred) * (1 - exp(-exp(cv.prob))) / (1 - stats::dpois(0, exp(cv.pred)))
              cv.pred <- log(cv.pred)
            }
          } else {
            cv.pred <- NA
          }

          if (scale == "abund") {
            cv.pred <- exp(cv.pred) * scale.factor
          } else {
            cv.pred <- cv.pred + log(scale.factor)
          }
          effect.dat[, f] <- cv.pred
        }
        colnames(effect.dat) <- paste0("CV", 1:length(cv.model.list))

        # Going to have to homebrew the factor plots
        uppers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", probs = .95, na.rm = T)
        lowers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", probs = .05, na.rm = T)
        variance <- apply(X = effect.dat, MARGIN = 1, FUN = var, na.rm = T)

        out.dat <- data.frame(x = f.data[, 1], effect = main.pred, var = se.pred^2, cvvar = variance, upper = uppers, lower = lowers, effect.dat)
      } else {
        out.dat <- data.frame(x = f.data[, 1], effect = main.pred, var = se.pred^2)
      }

      effects.list[[list.index]] <- out.dat
      names(effects.list)[list.index] <- termf
      list.index <- list.index + 1
    }
  }
  return(effects.list)
}
