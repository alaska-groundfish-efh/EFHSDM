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
