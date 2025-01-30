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
