#' GAM Stats
#'
#' @description Make a quick and dirty jackknife estimate of the deviance explained by each variable
# It can be rather slow for complicated models
#' @param model a GAM model object
#' @param data a data frame; usually the same one used to fit the GAM model
#'
#' @return a named vector of decimal values indicating the percent contribution to deviance explained
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
      utils::setTxtProgressBar(pb, length(gam.terms)) # maybe change this to just 'terms'?
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

