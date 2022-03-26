# This script includes functions for fitting GAM and HGAM models

# Setting lists of required packages & installing it
# rpackages <- c("mgcv", "raster", "PresenceAbsence")
#
# # rJava, rgeos, maps
# which_not_installed <- which(rpackages %in% rownames(installed.packages()) == FALSE)
#
# if (length(which_not_installed) > 1) {
#   install.packages(rpackages[which_not_installed], dep = TRUE)
# }
# rm(rpackages, which_not_installed)
#
# require(mgcv)
# require(raster)
# require(PresenceAbsence)


#' Make a GAM formula
#'
#' @description Improved version designed to use the tables produced by the Autodetect functions
# allows for easily dropping terms for things like term selections and deviance explained.
#' @param yvar Name of dependent variable for gam models
#' @param gam.table Data frame of parameters for GAM formula
#' @param hgam Logical; do you want an hgam formula
#'
#' @return Returns a formula object, or list of formulas for hgam
#' @export
#'
#'
#' @examples
AssembleGAMFormula <- function(yvar, gam.table, hgam = F) {

  # logic to handle different possibilities to supply for hgam
  if (hgam) {
    if (is.data.frame(gam.table[[1]])) {
      gam.table1 <- gam.table[[1]]
    } else {
      gam.table1 <- gam.table
    }
    if (is.data.frame(gam.table[[2]])) {
      gam.table2 <- gam.table[[2]]
    } else {
      gam.table2 <- gam.table1
    }
    form.list <- list(gam.table1, gam.table2)
  } else {
    gam.table1 <- gam.table
    form.list <- list(gam.table1)
  }

  for (f in 1:length(form.list)) {
    g.table <- form.list[[f]]
    out.terms <- rep(NA, nrow(g.table))

    for (t in 1:nrow(g.table)) {
      if (g.table$type[t] == "smooth") {
        out.term0 <- paste0("s(", g.table$term[t])
      }
      if (g.table$type[t] == "factor") {
        out.term0 <- paste0("as.factor(", g.table$term[t])
      }
      if (g.table$type[t] == "offset") {
        out.term0 <- paste0("offset(", g.table$term[t])
      }

      if (is.na(g.table$term2[t]) == F) {
        out.term0 <- paste0(out.term0, ",", g.table$term2[t])
      }
      if (is.na(g.table$bs[t]) == F) {
        out.term0 <- paste0(out.term0, ",bs='", g.table$bs[t], "'")
      }
      if (is.na(g.table$k[t]) == F) {
        out.term0 <- paste0(out.term0, ",k=", g.table$k[t])
      }
      if (is.na(g.table$m[t]) == F & is.na(g.table$m2[t])) {
        out.term0 <- paste0(out.term0, ",m=", g.table$m[t])
      }
      if (is.na(g.table$m[t]) == F & is.na(g.table$m2[t]) == F) {
        out.term0 <- paste0(out.term0, ",m=c(", g.table$m[t], ",", g.table$m2[t], ")")
      }
      out.terms[t] <- paste0(out.term0, ")")
      if (g.table$type[t] == "linear") {
        out.terms[t] <- g.table$term[t]
      }
    }

    if (f == 1) {
      out.form <- list(stats::as.formula(paste0(yvar, " ~ ", paste(out.terms, collapse = " + "))))
    }
    if (f > 1) {
      out.form[[f]] <- stats::as.formula(paste0(" ~ ", paste(out.terms, collapse = " + ")))
    }
  }
  if (
    length(form.list) == 1) {
    return(out.form[[1]])
  } else {
    return(out.form)
  }
}



#' Fit the GAM SDM
#'
#' @description Fit a GAM that predicts abundance based on the data. One should  supply a formula directly
#' @details If verbose=T, printouts of the model selection process are produced.
# new version with expanded capabilities for handling additional 2D smooth variables
# this function has two methods of term selection that can be used individually or together.
# Select=T will allow the model to reduce smooth terms to less than one degree of freedom, those terms are then discarded.  However, it does not affect non-smooth terms. Select=T can cause models to fit very slowly.
# reduce=T applies to all terms and will eliminate terms that increase the gcv score in a stepwise-manner
#' @param data A data frame including observations and covariates
#' @param gam.formula GAM formula or appropriate character string, should include everything
#' @param reduce Logical; should the model selection based on UBRE/GCV be used
#' @param select Logical; should the built in process be allowed to remove terms with less than one degree of freedom
#' @param verbose Logical; should summary be printed after each iteration
#' @param family.gam a family that is accepted by the gam function
#' @param theta Numeric; an optional set theta value for nb models
#' @param link.fx character; the name of a valid link function
#'
#' @return Returns a fitted GAM model object
#' @export
#'
#' @examples
FitGAM <- function(data,
                   gam.formula = NULL,
                   reduce = F,
                   select = F,
                   verbose = F,
                   family.gam = "poisson",
                   theta = NA,
                   link.fx = NA) {

  # Assemble a formula based on the inputs
  gam.form <- stats::as.formula(gam.formula)
  species <- as.character(gam.form)[[2]]

  # need to detect an offset, in case one isn't supplied
  off.0 <- trimws(strsplit(as.character(gam.form)[[3]], split = "[+]")[[1]])
  off.1 <- unlist(lapply(strsplit(off.0, split = "[()]"), FUN = function(x) {
    return(ifelse(x[1] == "offset", x[2], NA))
  }))
  offset.gam <- off.1[is.na(off.1) == F]
  if (length(offset.gam) == 0) {
    offset.gam <- NULL
  }


  # autodetect things and construct a family, if link function isn't specified
  if (family.gam %in% c("gaussian", "quasigaussian") & is.na(link.fx)) {
    link.fx <- "identity"
  }
  if (family.gam %in% c("poisson", "quasipoisson", "nb") & is.na(link.fx)) {
    link.fx <- "log"
  }
  if (family.gam == "binomial") {
    if (is.na(link.fx)) {
      link.fx <- "logit"
    }
    data[, species] <- as.integer(data[, species] > 0)
  }
  if (is.na(theta) == F) {
    gam.fam <- eval(paste0(family.gam, "(theta=", theta, ",link=", link.fx, ")"))
  } else {
    gam.fam <- eval(paste0(family.gam, "(link=", link.fx, ")"))
  }

  # set up a list of terms for elimination
  xvar.gam <- trimws(unlist(strsplit(as.character(gam.form)[[3]], split = "[+]")))

  # pulls the offset out
  for (v in 1:length(xvar.gam)) {
    if (strsplit(xvar.gam[v], split = "[(]")[[1]][1] == "offset") {
      xvar.gam <- xvar.gam[-v]
      break
    }
  }
  if (is.null(offset.gam) == F) {
    offset.form <- paste0("offset(", offset.gam, ")")
  } else {
    offset.form <- NULL
  }

  # this bit is a complicated fix to problem with the formula characters
  for (i in 1:length(xvar.gam)) {
    phrase <- strsplit(xvar.gam, split = " ")[[i]]
    spot <- which(phrase == "bs") + 2
    if (length(spot) > 0) {
      phrase2 <- strsplit(phrase[spot], "")[[1]]
      lets <- which(phrase2 %in% letters)
      newphrase <- paste0("'", paste(phrase2[lets], collapse = ""), "',")
      phrase[spot] <- newphrase
      xvar.gam[i] <- paste(phrase, collapse = "")
    }
  }

  if (select) {
    out.gam <- mgcv::gam(gam.form, family = gam.fam, data = data, select = T)
    if (verbose) {
      print(summary(out.gam))
    }

    # test if any terms were eliminated, and refit without them, iterating as necessary
    worst.var <- which.min(summary(out.gam)$edf)
    worst.name <- names(summary(out.gam)$s.table[, 1])[worst.var]
    target.edf <- length(unlist(strsplit(worst.name, split = "[(,)]"))) - 1

    if (summary(out.gam)$edf[worst.var] < 1) {
      bad.vars <- worst.var
    }
    if (summary(out.gam)$edf[worst.var] >= 1) {
      bad.vars <- NA
    }

    while (is.na(bad.vars) == F) {
      # reconstruct the formula
      xvar.gam <- xvar.gam[-bad.vars]
      gam.form <- stats::as.formula(paste0(paste(species, "~", paste(c(xvar.gam, offset.form), collapse = "+"))))

      # re fit and check
      out.gam <- mgcv::gam(gam.form, family = gam.fam, data = data, select = T)
      worst.var <- which.min(summary(out.gam)$edf)
      if (summary(out.gam)$edf[worst.var] < 1) {
        bad.vars <- worst.var
      }
      if (summary(out.gam)$edf[worst.var] >= 1) {
        bad.vars <- NA
      }
      if (verbose) {
        print(summary(out.gam))
      }
    }
  } else {
    out.gam <- mgcv::gam(gam.form, family = gam.fam, data = data, select = F)
  }

  if (reduce) {
    # stepwise variable selection based on ubre
    while (reduce) {
      gcv_gam <- out.gam$gcv.ubre
      pvals <- summary(out.gam)$s.pv
      pvals <- c(pvals, summary(out.gam)$p.pv[-1])
      least_sig <- which.max(pvals)

      xvar.gam1 <- xvar.gam[-least_sig]
      gam.form1 <- stats::as.formula(paste0(paste(species, "~", paste(c(xvar.gam1, offset.form), collapse = "+"))))
      gam.form2.1 <- stats::as.formula(paste0(as.character(gam.form1)[1], as.character(gam.form1)[3]))

      out.gam1 <- mgcv::gam(gam.form1, family = gam.fam, data = data)

      gcv_gam1 <- out.gam1$gcv.ubre

      # Fixed bug causing xvars not to update
      if (gcv_gam > gcv_gam1) {
        gam.form <- gam.form1
        gcv_gam <- gcv_gam1
        xvar.gam <- xvar.gam1
        out.gam <- out.gam1
      }
      if (verbose == T) {
        print(summary(out.gam))
      }
      if (gcv_gam <= gcv_gam1) {
        reduce <- F
      }
    }
  }
  return(out.gam)
}


#' Fit the Hurdle GAM
#'
#' @description The ziplss models are different enough to need their own function, but otherwise behave the same
# due to some difficulties with the select option. It is only applied to the models individually, and
# reduce is applied to both.
#' @param data A data frame including observations and covariates
#' @param density.formula formula for the abundance component
#' @param prob.formula formula for the probability component
#' @param reduce Logical; should the model selection based on UBRE/GCV be used
#' @param select Logical; should the built in process be allowed to remove terms with less than one degree of freedom
#' @param verbose Logical; should a summary be printed after each iteration
#'
#' @return returns a fitted hurdle model object i.e. list of two models
#' @export
#'
#' @examples
FitHurdleGAM <- function(data,
                         density.formula = NULL,
                         prob.formula = NULL,
                         select = F,
                         reduce = F,
                         verbose = F) {
  dens.form0 <- stats::as.formula(density.formula)
  species <- as.character(dens.form0)[2]

  # term selection will be carried out separately for the prob model and density model
  if (select | reduce) {
    prob.gam <- FitGAM(
      data = data, gam.formula = stats::as.formula(paste0(species, as.character(prob.formula)[1], as.character(prob.formula)[2])),
      reduce = reduce, select = select, verbose = verbose, family.gam = "binomial", link.fx = "cloglog"
    )
    prob.form <- stats::as.formula(paste(as.character(stats::formula(prob.gam))[-2], collapse = ""))

    dens.gam <- FitGAM(
      data = data, gam.formula = dens.form0, family.gam = "poisson", reduce = reduce,
      select = select, verbose = verbose
    )
    dens.form <- stats::formula(dens.gam)
  } else {
    dens.form <- stats::as.formula(density.formula)
    prob.form <- stats::as.formula(prob.formula)
  }

  out.gam <- mgcv::gam(list(dens.form, prob.form), family = "ziplss", data = data, select = F)


  # additional term selections for the factors is carried out together
  if (reduce) {
    # stepwise variable selection based on ubre
    while (reduce) {
      gcv_gam <- out.gam$gcv.ubre

      gam.sum <- summary(out.gam)
      # need some extra checks to handle the fact that there are two formulas mixed up in the table
      dens.xvars <- AutodetectGAMTerms(out.gam)[[1]]
      prob.xvars <- AutodetectGAMTerms(out.gam)[[2]]

      smooth.names <- names(gam.sum$chi.sq)
      fac.names <- names(gam.sum$p.pv)
      reduce.table <- data.frame(names = c(smooth.names, fac.names), prob = c(gam.sum$s.pv, gam.sum$p.pv))
      reduce.table <- reduce.table[reduce.table$names %in% c("(Intercept)", "(Intercept).1") == F, ]
      drop.var <- reduce.table$names[which.max(reduce.table$prob)]
      drop.name <- strsplit(drop.var, split = "[(,)]")[[1]][2]
      drop.index <- which.max(reduce.table$prob)

      smooth.fac <- ifelse(drop.index <= length(smooth.names), "smooth", "fac")

      if (smooth.fac == "smooth") {
        prob.dens <- ifelse(strsplit(drop.var, split = "[(]")[[1]][1] == "s", "dens", "prob")
      }
      if (smooth.fac == "fac") {
        prob.dens <- ifelse(strsplit(drop.var, split = "[)]")[[1]][1] == "1", "dens", "prob")
      }
      if (prob.dens == "dens") {
        dens.xvars <- dens.xvars[dens.xvars$term != drop.name, ]
      }
      if (prob.dens == "prob") {
        prob.xvars <- prob.xvars[prob.xvars$term != drop.name, ]
      }

      new.form <- AssembleGAMFormula(yvar = species, gam.table = list(dens.xvars, prob.xvars), hgam = T)
      out.gam1 <- mgcv::gam(new.form, family = "ziplss", data = data)
      gcv_gam1 <- out.gam1$gcv.ubre

      # Fixed bug causing xvars not to update
      if (gcv_gam > gcv_gam1) {
        out.gam <- out.gam1
      }
      if (verbose == T) {
        print(summary(out.gam))
      }
      if (gcv_gam <= gcv_gam1) {
        reduce <- F
      }
    }
  }
  return(out.gam)
}


#' Auto-detect GAM terms
#'
#' @description This is a useful function that returns a table of the terms in a GAM, whether they are a one dimensional (linear) term, a two dimensional smooth, or a factor. Also detects the offset.
#' @param model A fitted GAM model object
#' @param hgam character describing how to handle hgams; use "d" for the density model, "p" for the prob model, or "b" or "all" for both
#'
#' @return returns a data frame with columns describing the name and type of all model components and relevant smoother parameters
#' @export
#'
#' @examples
AutodetectGAMTerms <- function(model, hgam = "all") {
  n.formulas <- 1
  if (model$family$family == "ziplss" & hgam %in% c("b", "both", "all")) {
    n.formulas <- 2
  }

  for (f in 1:n.formulas) {
    form.index <- 3

    # a lot of special handling for ziplss models
    if (model$family$family == "ziplss") {
      if (hgam %in% c("b", "both", "all")) {
        form1 <- stats::formula(model)[[f]]
      }
      if (hgam %in% c("d", "dens", "density")) {
        form1 <- stats::formula(model)[[1]]
      }
      if (hgam %in% c("p", "prob", "probability")) {
        form1 <- stats::formula(model)[[2]]
      }
    } else {
      form1 <- stats::formula(model)
    }

    terms <- trimws(strsplit(as.character(form1[[length(form1)]]), split = "[+]")[[2]])

    # loop through and figure out the information for the table
    type.dat <- data.frame(type = rep(NA, length(terms)), dims = 1, term = NA, term2 = NA, bs = NA, k = NA, m = NA, m2 = NA)
    for (t in 1:length(terms)) {
      x <- terms[t]
      x2 <- strsplit(x, split = "[(=)]")[[1]]
      # linear terms
      if (length(x2) == 1) {
        type.dat$type[t] <- "linear"
        type.dat$term[t] <- x2
      } else {
        # smoothed terms
        if (x2[1] %in% c("s", "te")) {
          dims <- length(strsplit(x2[2], split = ", ")[[1]]) - 1
          for (n in 1:dims) {
            type.dat[t, 2 + n] <- strsplit(x2[2], split = ", ")[[1]][n]
          }
          type.dat$type[t] <- "smooth"
          type.dat$dims[t] <- dims

          # find smoother basis
          formula.options <- strsplit(x, split = ",")[[1]]
          bs.spot <- which(unlist(lapply(strsplit(formula.options, "="), FUN = function(x) {
            return(trimws(x[1]))
          })) == "bs")
          if (length(bs.spot) > 0) {
            type.dat$bs[t] <- strsplit(formula.options[bs.spot], split = "\"")[[1]][2]
          }

          # find smoother k
          k.spot <- which(unlist(lapply(strsplit(formula.options, "="), FUN = function(x) {
            return(trimws(x[1]))
          })) == "k")
          if (length(k.spot) > 0) {
            type.dat$k[t] <- trimws(strsplit(strsplit(formula.options[k.spot], split = "=")[[1]][2], split = "[)]")[[1]])
          }

          # find penalty m, which can be complicated
          m.spot <- which(unlist(lapply(strsplit(formula.options, "="), FUN = function(x) {
            return(trimws(x[1]))
          })) == "m")
          if (length(m.spot) > 0) {
            if ("c" %in% strsplit(formula.options[m.spot], split = "")[[1]]) {
              m.spot <- c(m.spot, m.spot + 1)
            }
            for (n in 1:length(m.spot)) {
              m1 <- trimws(formula.options[m.spot][n])
              if (n == 1) {
                m2 <- trimws(strsplit(m1, split = "=")[[1]][2])
              } else {
                m2 <- m1
              }
              m3 <- trimws(strsplit(m2, split = "[()]")[[1]])
              type.dat[t, 6 + n] <- suppressWarnings(stats::na.omit(as.numeric(m3))[1])
            }
          }
        }
        # factor terms
        if (x2[1] == "as.factor") {
          type.dat$type[t] <- "factor"
          type.dat$term[t] <- x2[2]
        }
      }
    }
    # Make the table
    terms2 <- unlist(strsplit(x = names(model$model), split = "[()]"))
    off.term <- which(terms2 == "offset") + 1
    if (length(off.term) > 0) {
      type.dat <- rbind(type.dat, data.frame(type = "offset", dims = 1, term = terms2[off.term], term2 = NA, bs = NA, k = NA, m = NA, m2 = NA))
    }
    # if multiple formulas, need to make a list
    if (n.formulas > 1) {
      if (f == 1) {
        out.dat <- list(type.dat)
      } else {
        out.dat[[f]] <- type.dat
      }
    } else {
      out.dat <- type.dat
    }
  }
  return(out.dat)
}

#' GAM Stats
#'
#' @description Make a quick and dirty jackknife estimate of the deviance explained by each variable
# It can be rather slow for complicated models
#' @param model a GAM model object
#' @param data a data frame; usually the same one used to fit the GAM model
#'
#' @return a named vector of decimal values indicating the percent contribution
#' @export
#'
#' @examples
GAMStats <- function(model, # a gam model
                     data) { # the data used to fit the model

  # set up some basics
  if (strsplit(model$family$family, split = "[()]")[[1]][1] == "Negative Binomial") {
    gamfam <- "nb"
  } else {
    gamfam <- model$family
  }

  if (model$family$family == "ziplss") {
    species <- as.character(stats::formula(model)[[1]])[[2]]
    gamfam <- "ziplss"
  } else {
    species <- as.character(stats::formula(model))[2]
  }

  if (model$family$family == "binomial") {
    data[, species] <- as.integer(data[, species] > 0)
  }

  # Need some code to get the terms and format them
  terms <- AutodetectGAMTerms(model)

  if (gamfam[1] == "ziplss") {
    d.terms <- terms[[1]]
    p.terms <- terms[[2]]
    x.vars <- unique(c(d.terms$term[d.terms$type != "offset"], p.terms$term[p.terms$type != "offset"]))
  } else {
    x.vars <- unique(terms$term[terms$type != "offset"])
  }

  # Progress bar; this is usually quick, but can sometimes get bogged down by convergence issues
  pb <- utils::txtProgressBar(min = 0, max = length(x.vars), style = 3)

  # loop to estimate the deviance explained by each term]
  dev.vec <- rep(NA, times = length(x.vars))
  for (i in 1:length(x.vars)) {

    # build a formula and fit a new gam
    if (gamfam[1] == "ziplss") {
      d.terms1 <- d.terms[d.terms$term != x.vars[i], ]
      p.terms1 <- p.terms[p.terms$term != x.vars[i], ]

      new.gam.form <- AssembleGAMFormula(yvar = species, gam.table = list(d.terms1, p.terms1), hgam = T)
      try(new.gam <- mgcv::gam(new.gam.form, family = "ziplss", data = data, select = F))

    } else {
      terms1 <- terms[terms$term != x.vars[i], ]
      new.gam.form <- AssembleGAMFormula(yvar = species, gam.table = terms1, hgam = F)
      try(new.gam <- mgcv::gam(stats::as.formula(new.gam.form), family = gamfam, data = data))
    }

    # evaluate
    if (exists("new.gam")) {
      dev.vec[i] <- summary(new.gam)$dev.expl
      rm(new.gam)
      # update progress bar
      utils::setTxtProgressBar(pb, i)
    } else {
      utils::setTxtProgressBar(pb, length(gam.terms))
      print("Deviance could not be estimated for model; returning NAs")
      break
    }
  }
  close(pb)

  # Now need to label the deviances
  if (sum(is.na(dev.vec)) == 0) {
    dev.lost <- 1 - dev.vec / summary(model)$dev.expl
    dev.exp <- dev.lost * 1 / sum(dev.lost, na.rm = T) * 100


    if (gamfam[1] == "ziplss") {
      name.terms <- unique(rbind(terms[[1]], terms[[2]]))
    } else {
      name.terms <- terms[terms$type != "offset", c("term", "term2")]
    }

    dev.names <- vector(length = length(x.vars))
    for (i in 1:length(x.vars)) {
      x.var.row <- which(name.terms$term == x.vars[i])
      if (is.na(name.terms$term2[x.var.row]) == F) {
        dev.names[i] <- paste0(name.terms$term[x.var.row], "*", name.terms$term2[x.var.row])
      } else {
        dev.names[i] <- name.terms$term[x.var.row]
      }
    }
    names(dev.exp) <- dev.names

    # sometimes you end up with negative deviance, so correct that
    if (min(dev.exp) < 0) {
      dev.exp <- dev.exp - min(dev.exp)
      dev.exp <- dev.exp / sum(dev.exp) * 100
    }
  } else {
    dev.exp <- NA
  }
  return(dev.exp)
}


#' Make GAM abundance
#'
#' @description This function makes abundance rasters for the GAMs. When an offset is present in the model, generates a map for the mean value of the offset
#' @param model a mgcv GAM model object
#' @param r.stack  raster stack with all the necessary covariates
#' @param scale.factor numeric; scale factor that will multiply the model predictions
#' @param filename filename for the output
#' @param land raster; a mask to be applied to the output usually landmasses
#' @param mask raster; an additional mask to apply to the output, if needed
#'
#' @return raster; usually representing abundance
#' @export
#'
#' @examples
MakeGAMAbundance <- function(model,
                             r.stack,
                             scale.factor = 1,
                             filename = "",
                             land = NULL,
                             mask = NULL) {

  # correct a common mistake
  if (is.na(filename) | is.null(filename)) {
    filename <- ""
  }

  model.terms <- AutodetectGAMTerms(model)
  if (model$family$family == "ziplss") {
    dterms <- model.terms[[1]]
    pterms <- model.terms[[2]]
    model.terms <- unique(rbind(dterms, pterms))
  }

  # now, need to add a check for the offset, and make a dummy raster if there is one
  if ("offset" %in% model.terms$type) {
    off.name <- model.terms$term[which(model.terms$type == "offset")]
    off.val <- ifelse(is.list(model$offset), mean(model$offset[[1]]), mean(model$offset))
    off.raster <- raster::raster(
      ext = r.stack@extent, crs = r.stack@crs, nrow = r.stack@nrows, ncol = r.stack@ncols,
      vals = off.val
    )
    names(off.raster) <- off.name
    r.stack <- raster::stack(list(r.stack, off.raster))
  }

  # will also need to detect factors and format them into a list
  gam.factors <- model.terms$term[model.terms$type == "factor"]
  gam.factors2 <- list()
  if (length(gam.factors) > 0) {
    for (t in 1:length(gam.factors)) {
      range <- raster::subset(x = r.stack, subset = which(names(r.stack) == gam.factors[t]))
      gam.factors2[[t]] <- sort(unique(stats::na.omit(raster::getValues(range))))
    }
    names(gam.factors2) <- gam.factors
  } else {
    gam.factors2 <- NULL
  }

  # Detect the link function and make the predictions
  pred.type <- ifelse(model$family$link == "cloglog", "link", "response")[1]

  # predict raster function malfunctions if ziplss model has no factors, so need a special case
  if(model$family$family == "ziplss" & is.null(gam.factors2)){
    r.vals<-as.data.frame(raster::getValues(r.stack))
    if ("offset" %in% model.terms$type) {
      r.vals<-cbind(r.vals,data.frame(off.val))
      names(r.vals[ncol(r.vals)])<-off.name
    }
    pred.vals<-mgcv::predict.gam(model,newdata=r.vals,type="response")
    predict.raster<-raster::setValues(x = raster::raster(r.stack),values = pred.vals)
  }else{
    predict.raster <- raster::predict(r.stack, model,factors = gam.factors2, progress = "text",
                              overwrite = TRUE,type = pred.type, newdata.guaranteed = TRUE)
  }
  # Detects and apply the cloglog approximation if appropriate
  if (model$family$link[1] == "cloglog") {
    predict.raster <- exp(predict.raster)
  }

  predict.raster <- predict.raster * scale.factor

  if (filename != "") {
    raster::writeRaster(x = predict.raster, filename = filename, overwrite = TRUE)
  }

  # now apply masks if appropriate
  if (is.null(land) == F) {
    predict.raster <- raster::mask(predict.raster, land,
                                   inverse = TRUE, overwrite = TRUE,
                                   filename = filename
    )
  }
  if (is.null(mask) == F) {
    predict.raster <- raster::mask(predict.raster, mask,
                                   overwrite = TRUE,
                                   filename = filename
    )
  }
  return(predict.raster)
}

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
  if (is.null(cv.model.list) == F && is.na(cv.model.list)) {
    cv.model.list <- NULL
  }

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
