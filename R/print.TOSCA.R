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
      crayon::bgRed(crayon::yellow("* TOSCA *")),
      "{.field {x$Input$Therapies %>% pull(Name) %>% unique() %>% length}} therapies, ",
      "{.field {x$Input$Mutations %>% dplyr::filter(Type == 'alpha') %>% nrow()}} CNAs, ",
      "{.field {x$Input$Mutations %>% dplyr::filter(Type == 'driver') %>% nrow()}} drivers."
    )
  )

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


# # cli::cli_alert_info(paste0(" CNA segments: ", x$n_cna_clonal, " clonal, ", x$n_cna_sbclonal, " subclonal."))
#
# # cli::cli_alert_info(paste0("Mutation mapping (head): ", paste0(head(x$n_karyotype), ' (',
# # names(head(x$n_karyotype)), ')', collapse = '; ')))
#
# # cli::cli_alert_info(paste0("Mutation mapping (up to top 5): "))
#
# cli::cli_h3("Clonal CNAs")
# cat('\n')
# bar_print_console(x, top = 10)
# cat('\n')
#
# if(x$n_cna_subclonal > 0)
# {
#   cli::cli_h3("Subclonal CNAs (showing up to 10 segments)")
#   cat('\n')
#   bar_print_console_scl(x, top = 10)
#   cat('\n')
# }
#
# # Available analyses
# with_peaks = all(!is.null(x$peaks_analysis))
# with_CCF = all(!is.null(x$CCF_estimates))
# with_smoothing = all(!is.null(x$before_smoothing))
# with_arm_frag = all(!is.null(x$arm_fragmentation))
# with_wg_frag = all(!is.null(x$wg_fragmentation))
# with_drivers = x %>% has_driver_data()
# with_MAF = x %>% has_MAF_annotations()
#
# cli::cli_alert_info(paste0(
#   "Sample Purity: ",
#   paste0(x$purity  * 100, '% ~ Ploidy: ', x$ploidy, '.')
# ))
#
# if (with_drivers |
#     with_peaks |
#     with_CCF | with_smoothing | with_arm_frag | with_wg_frag | with_MAF)
#   cat('\n')
#
# if(with_drivers)
# {
#   nd = x$mutations %>% dplyr::filter(is_driver) %>% nrow()
#   cli::cli_alert_info("There are {.value {nd}} annotated driver(s) mapped to clonal CNAs.")
#
#   w_d = x$mutations %>%
#     dplyr::filter(is_driver) %>%
#     dplyr::select(chr, from, to, ref, alt, DP, NV, VAF, driver_label, is_driver) %>%
#     as.data.frame()
#
#   writeLines(paste0("      ",
#                     capture.output( w_d %>%
#                                       print(row.names = F)
#                     )))
# }
#
# if(with_MAF)
# {
#   cli::cli_alert("This sample seem to have MAF columns annotated")
# }
#
# ppass = function()
#   "{crayon::bgGreen(crayon::black(\" PASS \"))}"
# pfail = function()
#   "{crayon::bgRed(crayon::white(\" FAIL \"))}"
#
#
# if (with_peaks)
# {
#   if ("matches" %in% names(x$peaks_analysis))
#   {
#     prop = x$peaks_analysis$matches %>%
#       dplyr::group_by(QC) %>%
#       dplyr::summarise(prop = sum(weight), .groups = 'drop') %>%
#       dplyr::arrange(dplyr::desc(prop)) %>%
#       dplyr::filter(dplyr::row_number() == 1) %>%
#       dplyr::pull(prop)
#     prop = round(prop * 100, digits = 0)
#
#     pur_sc = round(x$peaks_analysis$score, digits = 4)
#     pur_ch = paste0(round(x$peaks_analysis$score * 100, digits = 0), '%')
#
#     if (x$peaks_analysis$QC == "PASS")
#       cli::cli_h1(
#         paste0(
#           ppass(),
#           " Peaks QC {crayon::bold(x$peaks_analysis$matching_strategy)}: {crayon::green(paste0(prop, '%'))}, {crayon::green(paste('\u03bb =', pur_sc))}. Purity correction: {.value {pur_ch}}."
#         )
#       )
#     else
#       cli::cli_h1(
#         paste0(
#           pfail(),
#           " Peaks QC {crayon::bold(x$peaks_analysis$matching_strategy)}: {crayon::red(paste0(prop, '%'))}, {crayon::red(paste('\u03bb =', pur_sc))}. Purity correction: {.value {pur_ch}}."
#         )
#       )
#
#     cat('\n')
#
#     xx = x$peaks_analysis$matches %>%
#       dplyr::mutate(QC = ifelse(QC == "PASS", ppass(), pfail()))
#
#     for (karyo in xx$karyotype %>% unique)
#     {
#       xxx = xx %>%
#         dplyr::filter(karyotype == karyo)
#
#       qc = paste(xxx$QC, sprintf("%-7s", round(xxx$offset, 3)), collapse = ' ')
#       # qc = paste0("[", xx$karyotype[1], "]",  qc)
#
#       n = sprintf("%-5s", x$n_karyotype[karyo])
#       p = x$n_karyotype[karyo] / sum(x$n_karyotype[xx$karyotype %>% unique])
#       p = format(p * 100, digits = 1)
#       p = sprintf("%3s", p)
#
#       cli::cli_alert_info(paste0(
#         crayon::blue(xxx$karyotype[1]),
#         " ~ n = {n} ({p}%) {clisymbols::symbol$arrow_right} ",
#         qc,
#         ""
#       ))
#
#     }
#   }
#
#   # General karyotypes
#   if('general' %in% (x$peaks_analysis %>% names))
#   {
#     gen =  x$peaks_analysis$general$summary %>%
#       dplyr::ungroup() %>%
#       dplyr::mutate(prop = round(n/sum(n) * 100, 0))
#
#     n_matched =  gen$matched %>% sum
#     n_mismatched =  gen$mismatched %>% sum
#
#     cli::cli_h1(
#       paste(
#         "General peak QC ({.field {sum(gen$n)}} mutations):", ppass(), n_matched, pfail(), n_mismatched, "- epsilon = {.value {x$peaks_analysis$general$params$epsilon}}."
#       )
#     )
#     cat('\n')
#
#     for(i in 1:nrow(gen))
#     {
#       qc = paste(
#         sprintf("%-7s", paste(ppass(), gen$matched[i])),
#         sprintf("%-7s", paste(pfail(), gen$mismatched[i]))
#       )
#
#       n = sprintf("%-5s", gen$n[i])
#       p = sprintf("%3s", gen$prop[i])
#
#       cli::cli_alert_info(
#         paste0(
#           crayon::blue(gen$karyotype[i]),
#           " ~ n = {n} ({p}%) {clisymbols::symbol$arrow_right} ",
#           qc,
#           ""
#         )
#       )
#     }
#   }
#
#   # Subclonal
#   if('subclonal' %in% names(x$peaks_analysis))
#   {
#     nsegs =  x$peaks_analysis$subclonal$expected_peaks$segment_id %>% unique %>% length()
#
#     general_table = x$peaks_analysis$subclonal$summary %>%
#       # filter(prop > 0) %>%
#       dplyr::group_by(segment_id) %>%
#       dplyr::summarise(BR = sum(model == 'branching' & prop > 0), LI = sum(model == 'linear'& prop > 0)) %>%
#       dplyr::arrange(dplyr::desc(BR), dplyr::desc(LI))
#
#     certainly_linear = general_table %>% dplyr::filter(BR == 0, LI > 0) %>% nrow()
#     certainly_branching = general_table %>% dplyr::filter(BR > 0, LI == 0) %>% nrow()
#     ambiguous = general_table %>% dplyr::filter(BR > 0, LI > 0) %>% nrow()
#
#     all_bad = nsegs - (certainly_linear + certainly_branching + ambiguous)
#
#     numcol = function(x)
#     {
#       if(x == 0) return(crayon::red("{0}"))
#       return(paste("{.field {", eval(parse(text=x)), "}}"))
#     }
#
#     cli::cli_h1(
#       paste(
#         "Subclonal peaks QC ({.field {nsegs}} segments, initial state {.field {x$most_prevalent_karyotype}}):",
#         crayon::underline("linear"), certainly_linear %>% numcol(),
#         crayon::underline("branching"), certainly_branching%>% numcol(),
#         crayon::underline("either"), ambiguous%>% numcol(),
#         crayon::underline("no support"), all_bad%>% numcol(),
#         "- epsilon = {.value {x$peaks_analysis$subclonal$params$epsilon}}."
#       )
#     )
#
#     puncertain = function()
#       "{crayon::bgYellow(crayon::white(\" UNKNOWN \"))}"
#
#     my_print = function(s, b){
#       # cat(s)
#
#       hb = x$peaks_analysis$subclonal$summary %>%
#         dplyr::filter(segment_id == s) %>%
#         dplyr::arrange(dplyr::desc(prop)) %>%
#         dplyr::mutate(prop = prop * 100) %>%
#         dplyr::mutate(prop = case_when(
#           prop == 0 ~ crayon::red('0'),
#           prop == 100 ~ crayon::green('100'),
#           TRUE ~ crayon::yellow(prop)
#         )) %>%
#         dplyr::mutate(label = paste0(model_id, " [", prop, "]")) %>%
#         dplyr::pull(label) %>%
#         paste(collapse = '; ')
#
#
#       n = x$peaks_analysis$subclonal$summary %>% filter(segment_id == s) %>% filter(row_number() == 1) %>% dplyr::pull(size)
#       cl = x$peaks_analysis$subclonal$summary %>% filter(segment_id == s) %>% filter(row_number() == 1) %>% dplyr::pull(clones)
#       cl = strsplit(cl, ' ')[[1]]
#       cl = paste0(cl[1], ' (', cl[2] %>% as.numeric *100, ') + ', cl[3], ' (', cl[4] %>% as.numeric * 100, ')')
#       cl = sprintf("%17s", cl)
#
#       cli::cli_alert_info(
#         paste0(
#           crayon::blue(sprintf(paste0("%", b, 's'), s)),
#           " ~ {n[1]} {crayon::yellow(cl)} : ",
#           hb,
#           ""
#         )
#       )
#     }
#
#     ml = general_table$segment_id %>% nchar() %>% max
#
#     if(certainly_linear > 0)
#     {
#       cli::cli_h3(paste(ppass(), "Linear models"))
#
#       general_table %>%
#         dplyr::filter(BR == 0, LI > 0) %>%
#         dplyr::pull(segment_id) %>%
#         lapply(my_print, b = ml)
#     }
#
#     if(certainly_branching > 0)
#     {
#       cli::cli_h3(paste(ppass(), "Branching models"))
#
#       general_table %>%
#         dplyr::filter(BR > 0, LI == 0) %>%
#         dplyr::pull(segment_id) %>%
#         lapply(my_print, b = ml)
#     }
#
#     if(ambiguous > 0)
#     {
#       cli::cli_h3(paste(puncertain(), "Either branching or linear models"))
#
#       general_table %>%
#         dplyr::filter(BR > 0, LI > 0) %>%
#         dplyr::pull(segment_id) %>%
#         lapply(my_print, b = ml)
#     }
#
#     if(all_bad > 0)
#     {
#       cli::cli_h3(paste(pfail(), "Bad models"))
#
#       general_table %>%
#         dplyr::filter(BR == 0, LI == 0) %>%
#         dplyr::pull(segment_id) %>%
#         lapply(my_print, b = ml)
#     }
#
#     # nlin = x$peaks_analysis$subclonal$expected_peaks %>%
#     #   filter(model == 'linear', matched) %>% nrow()
#     # nlin_f = x$peaks_analysis$subclonal$expected_peaks %>%
#     #   filter(model == 'linear', !matched) %>% nrow()
#     #
#     # nbr = x$peaks_analysis$subclonal$expected_peaks %>%
#     #   filter(model == 'branching', matched) %>% nrow()
#     # nbr_f = x$peaks_analysis$subclonal$expected_peaks %>%
#     #   filter(model == 'branching', !matched) %>% nrow()
#
#
#     # cli::cli_h3(
#     #   paste(
#     #     "Subclonal peak QC ({.field {nsegs}} segments):",
#     #     crayon::bold("linear"),
#     #     ppass(), nlin, pfail(), nlin_f,
#     #     "~",
#     #     crayon::bold("branching"),
#     #     ppass(), nbr, pfail(), nbr_f,
#     #     "- epsilon = {.value {x$peaks_analysis$subclonal$params$epsilon}}."
#     #   )
#     # )
#
#     #   S_table = x$peaks_analysis$subclonal$summary
#     #   S_table$size = strsplit(S_table$segment_id, split = '\\n') %>%
#     #     sapply(function(x) x[[2]])
#     #   S_table$clones = strsplit(S_table$segment_id, split = '\\n') %>%
#     #     sapply(function(x) x[[3]])
#     #   S_table$segment_id = strsplit(S_table$segment_id, split = '\\n') %>%
#     #     sapply(function(x) x[[1]]) %>%
#     #     strsplit(split = ' ') %>%
#     #     sapply(function(x) { paste(x[1], x[2], sep = '@') })
#     #
#     #   for (s in S_table$segment_id %>% unique())
#     #   {
#     #     hb = S_table %>%
#     #       filter(segment_id == s) %>%
#     #       pull(prop) %>%
#     #       length() == 1
#     #
#     #     bp = S_table %>%
#     #       filter(segment_id == s) %>%
#     #       mutate(prop = prop * 100) %>%
#     #       pull(prop) %>%
#     #       round(0) %>%
#     #       paste0('%')
#     #
#     #     bm = S_table %>%
#     #       filter(segment_id == s) %>%
#     #       pull(model)
#     #
#     #     qc = ifelse(hb,
#     #                 crayon::green(paste(bm, bp, collapse = ', ')),
#     #                 paste(bm, bp, collapse = ', '))
#     #
#     #     n = S_table %>% filter(segment_id == s) %>% pull(size)
#     #     n = sprintf("%20s", n)
#     #     cl = S_table %>% filter(segment_id == s) %>% pull(clones)
#     #     cl = sprintf("%17s", cl)
#     #
#     #     cli::cli_alert_info(
#     #       paste0(
#     #         crayon::blue(sprintf("%17s", s)),
#     #         " ~ {n[1]} {crayon::yellow(cl[1])} {clisymbols::symbol$arrow_right} ",
#     #         qc,
#     #         ""
#     #       )
#     #     )
#     #
#     # }
#   }
# }
#
# if (with_CCF)
# {
#   cat("\n")
#   cli::cli_alert_success(
#     "Cancer Cell Fraction (CCF) data available for karyotypes:{.value {names(x$CCF_estimates)}}."
#   )
#
#   lapply(x$CCF_estimates,
#          function(l) {
#            l$QC_table
#          }) %>%
#     Reduce(f = bind_rows) %>%
#     apply(1, function(z) {
#       if (z['QC'] == "PASS")
#         cli::cli_alert_success(
#           "{
#             crayon::bgGreen(
#               crayon::white(
#                 \" PASS \"))} CCF via {.value {crayon::bold(z['method'])}}."
#         )
#       else
#         cli::cli_alert_success(
#           "{crayon::bgRed(crayon::white(\" FAIL \"))} CCF via {.value {crayon::bold(z['method'])}}."
#         )
#     })
# }
#
# if (with_smoothing)
#   cli::cli_alert_success(
#     "These segments are smoothed; before smoothing there were {.value {x$before_smoothing$n_cna}} segments."
#   )
#
# if (with_arm_frag)
#   cli::cli_alert_success(
#     "Arm-level fragmentation analysis: {.value {sum(x$arm_fragmentation$table$significant)}} segments overfragmented."
#   )
#
# if (with_wg_frag)
# {
#   cond = x$wg_fragmentation$is_overfragmented
#   p = round(x$wg_fragmentation$pvalue, 5)
#   cli::cli_alert_success(
#     "Whole-genome fragmentation analysis: p = {.value {p}}: {.value {ifelse(cond, crayon::red('overfragmented'), crayon::green('not overfragmented'))}}."
#   )
# }
#
# ppass = function()
#   "{crayon::bgGreen(crayon::black(\" PASS \"))}"
# pfail = function()
#   "{crayon::bgRed(crayon::white(\" FAIL \"))}"
#
#
