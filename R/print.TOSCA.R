#' Print for class \code{'TOSCA'}.
#'
#' @param x An obj of class \code{'TOSCA'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @exportS3Method print TOSCA
#' @export print.TOSCA
#'
#' @examples
#' \dontrun{
#' data('exampleFit', package = 'TOSCA')
#'
#' print(exampleFit)
#' }
print.TOSCA <- function(x, ...) {
  stopifnot(inherits(x, "TOSCA"))

  cli::cli_rule(
    paste(
      crayon::bgMagenta(crayon::white(paste0("** TOSCA ** Sample: ", x$Input$Sample_name, " ")))
    )
  )

  cli::cli_h3(paste0("Clinical History (", TOSCA:::get_clinical_history_length(x), ")"))
  #cat('\n')

  cli::cli_rule(
    paste0(
      "{.field {TOSCA:::get_first_sample_name(x)}} ",
      "({.field {TOSCA:::get_first_sample_date(x)}})",
      " --> ",
      "{.field {.white {TOSCA:::get_therapies_ordered(x)}}}",
      " --> ",
      "{.field {TOSCA:::get_second_sample_name(x)}} ",
      "({.field {TOSCA:::get_second_sample_date(x)}})"
    )
  )

  cat('\n')

  cli::cli_h3(paste("Timed events: ",
                   "{.field {TOSCA:::get_number_of_cna(x)}} CNAs {.field {TOSCA:::get_possible_dormancy(x)}} dormancy",
                   "{.field {TOSCA:::get_driver(x)}} Drivers altering the mutation rate"))

  cat("\n")
  ## Print summary of the clinical history

  # "Patient is [age] at [Name of first sample]"
  # "Patient is treated with:
  #  - XXX cycles of DrugXXX (XX-XX-XX, XX-XX-XX,XX-XX-XX ...), which can [induce Dormancy | induce Mutagenesis]
  #  - XXX cycles of DrugXXX (XX-XX-XX, XX-XX-XX,XX-XX-XX ...)  which can [induce Dormancy | induce Mutagenesis associated with [mutation name]] "
  # Patient is [age] when sample [Name of second sample] is collected.
  #
  # The tumour at [Name of second sample] is clonal for [XXX CNAs | driver] associated to [alpha-beta | sbs mutations]
  cli::cli_h3("Inference results")
  if (!is.null(x$Fit)){
    TOSCA::get_fit_summary(x)
    # ---- Print Diagnostics ----
    cat("\n--- Sampling Diagnostics ---\n")

    # Convergence: Check all Rhat < 1.01
    max_rhat <- max(x$Fit$summary$rhat, na.rm = TRUE)
    converged <- all(x$Fit$summary$rhat <= 1.01, na.rm = TRUE)
    cat(sprintf("Convergence (Rhat < 1.01): %s (max Rhat = %.3f)\n",
                if (converged) "\u2705 Yes" else "\u274C No", max_rhat))

    # Divergent transitions
    n_chains = length(x$Fit$diagnostic_summary$num_divergent)
    n_iterations = nrow(x$Fit$posteriors$stan_posterior) / n_chains
    divergences <- x$Fit$diagnostic_summary$num_divergent
    total_div <- sum(divergences)
    total_samples <- n_iterations*n_chains
    div_pct <- 100 * total_div / total_samples
    cat(sprintf("Divergent transitions: %d / %d (%.2f%%)\n", total_div, total_samples, div_pct))

    # EBFMI check
    ebfmi <- x$Fit$diagnostic_summary$ebfmi
    low_ebfmi <- sum(ebfmi < 0.3)
    cat(sprintf("EBFMI < 0.3 in %d / %d chains\n", low_ebfmi, length(ebfmi)))
    #cat("--- End Diagnostics ---\n")

    cat("\n--- Posterior Predictive Checks ---\n")
    print.data.frame(TOSCA:::check_ppc(x))

    cat("\n--- Inference summary ---\n")
    print.data.frame(TOSCA::get_fit_summary(x))
  }else{
    print("Fit not available")
  }
  invisible(x)
}


