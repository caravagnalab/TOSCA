# Plot clinical timeline + posterior times
#' Posterior distribution of the inferred times
#'
#' @param x TOSCA object
#'
#' @return Posterior distributio plot of the inferred times, mapped on the clinical history
#' @export
#'
#' @examples
#' library(TOSCA)
#' library(dplyr)
#' library(ggplot2)
#' data("exampleData_CNA")
#' mutations = exampleData_CNA$Mutations
#' parameters = exampleData_CNA$Parameters
#' samples = exampleData_CNA$Samples
#' therapies = exampleData_CNA$Therapies
#'
#' x = init(mutations=mutations, samples=samples, therapies=therapies, parameters=parameters)
#' fit = TOSCA::fit(x, model_name='CNA', n_iterations = 1000, n_chains = 4, warm_up = 500)
#' plot_timing(fit)
plot_timing = function(x)
{
  #clinical_timeline #= x$Input$Samples
  estimates = TOSCA:::get_inferred_parameters(x)

  endpoints = x$Input$Samples

  if (nrow(endpoints) > 2) endpoints = (endpoints %>% dplyr::arrange(as.Date(Date)))[2:3,]

  therapies = x$Input$Therapies

  # 1. time posterior plot
  timing_estimates = estimates %>%
    dplyr::select(starts_with('t_')) %>%
    apply(2, TOSCA:::convert_date_real, x=x) %>%
    dplyr::as_tibble()
  #times = timing_estimates$variable %>% unique()

  for (i in 1:ncol(timing_estimates)){
    timing_estimates[[i]] = as.Date(timing_estimates[[i]])
  }

  timing_estimates = timing_estimates %>% reshape2::melt() %>% dplyr::as_tibble()
  times = timing_estimates$variable %>% unique()
  therapy_names = therapies$Name %>% unique()
  var_colors = ggsci::pal_npg()(length(times)+length(therapy_names))
  times_colors = var_colors[1:length(times)]
  names(times_colors) = times
  clinical_colors = var_colors[length(times)+1:length(therapy_names)]
  names(clinical_colors) = therapy_names
  times_colors_df = data.frame("variable" = names(times_colors), "color"=times_colors)
  timing_estimates = dplyr::left_join(timing_estimates, times_colors_df, by = "variable")

  posterior_plot = ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = timing_estimates %>% dplyr::rename(Date = value),
      ggplot2::aes(Date, fill = color, ..density..),
      inherit.aes = FALSE,
      bins = 150
    ) +
    ggplot2::geom_point(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date), y = 0),
      inherit.aes = FALSE,
      size = 3
    ) +
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(legend.position = 'bottom')+
    ggplot2::scale_fill_identity()
  #scale_fill_manual(values = times_colors)

  hist_data <- ggplot2::ggplot_build(posterior_plot)$data[[1]]
  ymin <- 0
  ymax <- max(hist_data$y)

  # 5% above bottom
  ylab_pos <- ymin + 0.1 * (ymax - ymin)


  if ("t_dormancy_start" %in% timing_estimates$variable){
    MAP_dormancy_start = timing_estimates %>% dplyr::filter(variable == "t_dormancy_start") %>% dplyr::pull(value) %>% mean()
    MAP_dormancy_end = timing_estimates %>% dplyr::filter(variable == "t_dormancy_end") %>% dplyr::pull(value) %>% mean()
    timing_estimates = timing_estimates %>% dplyr::filter(!(variable %in% c("t_dormancy_start","t_dormancy_end")))

    posterior_plot = ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = timing_estimates %>% dplyr::rename(Date = value),
        ggplot2::aes(Date, fill = variable, ..density..),
        inherit.aes = FALSE,
        bins = 150
      ) +
      ggplot2::geom_point(
        data = endpoints,
        ggplot2::aes(x = as.Date(Start), y = 0),
        inherit.aes = FALSE,
        size = 3
      ) + my_ggplot_theme()+
      ggplot2::theme(legend.position = 'bottom')+
      ggplot2::scale_fill_manual(values = var_colors)+
      ggplot2::geom_rect(
        ggplot2::aes(xmin = MAP_dormancy_start,
                     xmax = MAP_dormancy_end),
        ymin = 0,
        ymax = Inf,
        fill = 'grey',
        colour = "white",
        size = 0.5,
        alpha = .5
      )
  }

  therapies = therapies %>% dplyr::mutate(Duration = as.Date(End)-as.Date(Start)) %>% dplyr::mutate(short=ifelse(Duration < 30, T, F))
  new_col = data.frame(Name = names(clinical_colors), colors = clinical_colors)
  therapies = dplyr::left_join(therapies,new_col, by="Name")

  posterior_plot = posterior_plot + ggplot2::geom_rect(
    data = therapies %>% dplyr::filter(short == F),
    ggplot2::aes(xmin = as.Date(Start),
                 xmax = as.Date(End),
                 fill=colors),
    ymin = 0,
    ymax = Inf,
    #fill = c,
    colour = "white",
    size = 0.5,
    alpha = .5
  ) +
    ggplot2::geom_segment(
      data = therapies %>% dplyr::filter(short == T),
      ggplot2::aes(x=as.Date(Start),
                   xend=as.Date(End),
                   y=0, yend=Inf), color = therapies %>% dplyr::filter(short == T) %>% dplyr::pull(colors), alpha=1)  +
    ggplot2::guides(color = "none")

  posterior_plot = posterior_plot +
    ggplot2::geom_segment(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date),
                   xend = as.Date(Date),
                   y=0,
                   yend = ylab_pos),
      size = .5, linetype="dashed"
    )+
    ggplot2::geom_label(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date), y = ylab_pos, label=Name),
      size = 3
    )

  dummy_guide <- function(
    labels = NULL,
    ...,
    title = NULL,
    key   = draw_key_point,
    guide_args = list(),
    min_value_time
  ) {
    aesthetics <- list(...)
    n <- max(lengths(aesthetics), 0)
    labels <- labels %||% seq_len(n)

    aesthetics$alpha <- aesthetics$alpha %||% rep(1, n)

    guide_args$override.aes <- guide_args$override.aes %||% aesthetics
    guide <- do.call(guide_legend, guide_args)

    ggplot2::update_geom_defaults("point", list(dummy = "x"))

    dummy_geom <- ggplot2::geom_point(
      data = data.frame(
        x = rep(min_value_time, n),
        y = rep(Inf, n),
        dummy = factor(labels, levels = labels)   # preserve order!
      ),
      ggplot2::aes(x, y, dummy = dummy),
      alpha = 0,
      key_glyph = key
    )

    fills <- aesthetics$fill

    dummy_scale <- ggplot2::discrete_scale(
      "dummy", "dummy_scale",
      palette = function(x) {
        stats::setNames(fills, labels)[x]
      },
      name = title,
      guide = guide
    )

    list(dummy_geom, dummy_scale)
  }


  posterior_plot +
    dummy_guide(
      labels = c(names(times_colors), names(clinical_colors)),
      fill   = c(times_colors, clinical_colors),
      colour = NA,
      title  = "Inferred times and treatments",
      key = draw_key_polygon,
      min_value_time = min(timing_estimates$value, na.rm = TRUE)
    )

}
