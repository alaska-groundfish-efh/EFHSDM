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
#' gam.form <- formula("a_atf ~ s(lon,lat,bs = 'ds',m=c(1,.5), k=10) +
#' s(bdepth, bs='tp',m=1,k=4) + s(btemp, bs='tp',m=1,k=4) +
#' s(slope, bs='tp',m=1,k=4) + offset(logarea)")
#' data("region_data_goa")
#' region_data_goa$sponge <- as.integer(region_data_goa$sponge > 0)
#' region_data_goa$logarea <- log(region_data_goa$area)
#'
#' poisson.model <- FitGAM(gam.formula = gam.form, data = region_data_goa, family.gam = "poisson")

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
