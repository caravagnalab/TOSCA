#' Plot clinical timeline + posterior times with density
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_timing_density = function(x)
{
  #clinical_timeline #= x$Input$Samples
  estimates = TOSCA:::get_inferred_parameters(x)

  endpoints = x$Input$Samples

  if (nrow(endpoints) > 2) endpoints = (endpoints %>% dplyr::arrange(as.Date(Date)))[2:3,]

  therapies = x$Input$Therapies

  # 1. time posterior plot
  timing_estimates = estimates %>%
    dplyr::select(starts_with('t_')) #%>%
    #apply(2, TOSCA:::convert_date_real, x=x) %>%
    #dplyr::as_tibble()
  #times = timing_estimates$variable %>% unique()

  if (x$Fit$model_info$dormancy) {
    dormancy_start = TOSCA:::convert_date_real(date = TOSCA:::get_start_therapy(x, class= "Chemotherapy inducing dormancy"), x=x)
    dormancy_end =  TOSCA:::convert_date_real(date = timing_estimates %>% pull(t_dormancy_end) %>% mean(), x=x)
    timing_estimates = timing_estimates %>% select(!c("t_dormancy_end", "t_mrca_tr",
                                                      colnames(timing_estimates)[grepl("t_cna_tr", colnames(timing_estimates))]))
  }

  timing_estimates = timing_estimates %>%
    apply(2, TOSCA:::convert_date_real, x=x) %>%
    dplyr::as_tibble()

  timing_estimates = TOSCA:::convert_timing_names(timing_estimates)

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


  # posterior_plot = ggplot2::ggplot() +
  #   ggplot2::geom_histogram(
  #     data = timing_estimates %>% dplyr::rename(Date = value),
  #     ggplot2::aes(Date, fill = color, y=..density..),
  #     inherit.aes = FALSE,
  #     bins = 150
  #   ) + #geom_density(data = timing_estimates %>% dplyr::rename(Date = value), aes(color = color)) +
  #   ggplot2::geom_point(
  #     data = endpoints,
  #     ggplot2::aes(x = as.Date(Date), y = 0),
  #     inherit.aes = FALSE,
  #     size = 3
  #   ) +
  #   TOSCA:::my_ggplot_theme()+
  #   ggplot2::theme(legend.position = 'bottom')+
  #   ggplot2::scale_fill_identity()
  #scale_fill_manual(values = times_colors)

  posterior_plot = ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = timing_estimates %>% dplyr::rename(Date = value),
      ggplot2::aes(Date, fill = color, y = ..density..),
      inherit.aes = FALSE,
      bins = 150,
      alpha = 0  # make histograms semi-transparent
      #color = "white"
    ) +
    ggplot2::geom_density(
      data = timing_estimates %>% dplyr::rename(Date = value),
      ggplot2::aes(Date, color = color, fill = color),
      inherit.aes = FALSE,
      size = 1, alpha = .8
    ) +
    ggplot2::geom_point(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date), y = 0),
      inherit.aes = FALSE,
      size = 3
    ) +
    TOSCA:::my_ggplot_theme() +
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()


  hist_data <- ggplot2::ggplot_build(posterior_plot)$data[[1]]
  ymin <- 0
  ymax <- max(hist_data$y)

  # 5% above bottom
  ylab_pos <- ymin + 0.1 * (ymax - ymin)


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
                   xend=as.Date(Start),
                   y=0, yend=Inf), color = therapies %>% dplyr::filter(short == T) %>% dplyr::pull(colors), alpha=1)  +
    ggplot2::guides(color = "none")

  # posterior_plot = posterior_plot +
  #   ggplot2::geom_segment(
  #     data = endpoints,
  #     ggplot2::aes(x = as.Date(Date),
  #                  xend = as.Date(Date),
  #                  y=0,
  #                  yend = ylab_pos),
  #     size = .5, linetype="dashed"
  #   )+
  #   ggplot2::geom_label(
  #     data = endpoints,
  #     ggplot2::aes(x = as.Date(Date), y = ylab_pos, label=Name),
  #     size = 3
  #   )

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

  if (nrow(therapies)>0){
    posterior_plot = posterior_plot +
      dummy_guide(
        labels = c(names(times_colors), names(clinical_colors)),
        fill   = c(times_colors, clinical_colors),
        colour = NA,
        title  = "Inferred times and treatments",
        key = draw_key_polygon,
        min_value_time = min(timing_estimates$value, na.rm = TRUE)
      )
  }else{
    posterior_plot = posterior_plot +
      dummy_guide(
        labels = names(times_colors),
        fill   = times_colors,
        colour = NA,
        title  = "Inferred times",
        key = draw_key_polygon,
        min_value_time = min(timing_estimates$value, na.rm = TRUE)
      )
  }

  if (x$Fit$model_info$dormancy) {

    dormancy_df = data.frame("event"= c("Dormancy"), "Start"=c(dormancy_start), "End"=c(dormancy_end))
    dormancy_df2 = data.frame("event"= c("Dormancy Start", "Dormancy End"), "date"=c(dormancy_start, dormancy_end))

    posterior_plot = posterior_plot +
      ggplot2::geom_rect(
        data = dormancy_df,
        ggplot2:::aes(xmin = as.Date(Start), xmax = as.Date(End)),
        ymin = 0,ymax = Inf, color = "white", fill="#c0c0c0ff", alpha = .5
      )+
      ggplot2::geom_segment(
        data = dormancy_df2,
        ggplot2::aes(x = as.Date(date),
                     xend = as.Date(date),
                     y=0,
                     yend = ylab_pos + 1.5*ylab_pos),
        size = .5, linetype="dashed", color = "#7d7d7dff"
      )+
      ggplot2::geom_label(
        data = dormancy_df2,
        ggplot2::aes(x = as.Date(date), y = ylab_pos + 1.5*ylab_pos, label=event),
        size = 3, fill="#c0c0c0ff"
      )
  }

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
  posterior_plot
}

#' Plot clinical timeline + posterior times with histogram
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_timing = function(x)
{
  #clinical_timeline #= x$Input$Samples
  estimates = TOSCA:::get_inferred_parameters(x)

  endpoints = x$Input$Samples

  if (nrow(endpoints) > 2) endpoints = (endpoints %>% dplyr::arrange(as.Date(Date)))[2:3,]

  therapies = x$Input$Therapies

  # 1. time posterior plot
  timing_estimates = estimates %>%
    dplyr::select(starts_with('t_')) #%>%
  #apply(2, TOSCA:::convert_date_real, x=x) %>%
  #dplyr::as_tibble()
  #times = timing_estimates$variable %>% unique()

  if (x$Fit$model_info$dormancy) {
    dormancy_start = TOSCA:::convert_date_real(date = TOSCA:::get_start_therapy(x, class= "Chemotherapy inducing dormancy"), x=x)
    dormancy_end =  TOSCA:::convert_date_real(date = timing_estimates %>% pull(t_dormancy_end) %>% mean(), x=x)
    timing_estimates = timing_estimates %>% select(!c("t_dormancy_end", "t_mrca_tr",
                                                      colnames(timing_estimates)[grepl("t_cna_tr", colnames(timing_estimates))]))
  }

  timing_estimates = timing_estimates %>%
    apply(2, TOSCA:::convert_date_real, x=x) %>%
    dplyr::as_tibble()

  timing_estimates = TOSCA:::convert_timing_names(timing_estimates)

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
      ggplot2::aes(Date, fill = color, y=..density..),
      inherit.aes = FALSE,
      bins = 150, color = "black",
    ) + #geom_density(data = timing_estimates %>% dplyr::rename(Date = value), aes(color = color)) +
    ggplot2::geom_point(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date), y = 0),
      inherit.aes = FALSE,
      size = 3
    ) +
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(legend.position = 'bottom',
                   axis.line.y = element_blank(),
                   panel.grid.major= element_blank(),
                   panel.grid.minor= element_blank(),
                   #axis.ticks.y= element_blank(),
                   #axis.text.y= element_blank(),
                   #axis.title.y= element_blank(),
                   panel.border = element_blank())+
    ggplot2::scale_fill_identity()
  #+
  #scale_fill_manual(values = times_colors)

  # posterior_plot = ggplot2::ggplot() +
  #   ggplot2::geom_histogram(
  #     data = timing_estimates %>% dplyr::rename(Date = value),
  #     ggplot2::aes(Date, fill = color, y = ..density..),
  #     inherit.aes = FALSE,
  #     bins = 150,
  #     alpha = 0  # make histograms semi-transparent
  #     #color = "white"
  #   ) +
  #   ggplot2::geom_density(
  #     data = timing_estimates %>% dplyr::rename(Date = value),
  #     ggplot2::aes(Date, color = color, fill = color),
  #     inherit.aes = FALSE,
  #     size = 1, alpha = .8
  #   ) +
  #   ggplot2::geom_point(
  #     data = endpoints,
  #     ggplot2::aes(x = as.Date(Date), y = 0),
  #     inherit.aes = FALSE,
  #     size = 3
  #   ) +
  #   TOSCA:::my_ggplot_theme() +
  #   ggplot2::theme(legend.position = 'bottom') +
  #   ggplot2::scale_fill_identity() +
  #   ggplot2::scale_color_identity()


  hist_data <- ggplot2::ggplot_build(posterior_plot)$data[[1]]
  ymin <- 0
  ymax <- max(hist_data$y)

  # 5% above bottom
  ylab_pos <- ymin + 0.1 * (ymax - ymin)


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
                   xend=as.Date(Start),
                   y=0, yend=Inf), color = therapies %>% dplyr::filter(short == T) %>% dplyr::pull(colors), alpha=1)  +
    ggplot2::guides(color = "none")

  rect_ind <- which(lapply(posterior_plot$layers, function(x) class(x$geom)[1]) == "GeomRect")
  ## manually change the order to put geom rect first
  posterior_plot$layers <- c(posterior_plot$layers[rect_ind], posterior_plot$layers[-rect_ind])


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

  if (nrow(therapies)>0){
    posterior_plot = posterior_plot +
      dummy_guide(
        labels = c(names(times_colors), names(clinical_colors)),
        fill   = c(times_colors, clinical_colors),
        colour = NA,
        title  = "Inferred times and treatments",
        key = draw_key_polygon,
        min_value_time = min(timing_estimates$value, na.rm = TRUE)
      )
  }else{
    posterior_plot = posterior_plot +
      dummy_guide(
        labels = names(times_colors),
        fill   = times_colors,
        colour = NA,
        title  = "Inferred times",
        key = draw_key_polygon,
        min_value_time = min(timing_estimates$value, na.rm = TRUE)
      )
  }

  if (x$Fit$model_info$dormancy) {

    dormancy_df = data.frame("event"= c("Dormancy"), "Start"=c(dormancy_start), "End"=c(dormancy_end))
    dormancy_df2 = data.frame("event"= c("Dormancy Start", "Dormancy End"), "date"=c(dormancy_start, dormancy_end))

    posterior_plot = posterior_plot +
      ggplot2::geom_rect(
        data = dormancy_df,
        ggplot2:::aes(xmin = as.Date(Start), xmax = as.Date(End)),
        ymin = 0,ymax = Inf, color = "white", fill="#c0c0c0ff", alpha = .5
      )
    # +
    #   ggplot2::geom_segment(
    #     data = dormancy_df2,
    #     ggplot2::aes(x = as.Date(date),
    #                  xend = as.Date(date),
    #                  y=0,
    #                  yend = ylab_pos + 1.5*ylab_pos),
    #     size = .5, linetype="dashed", color = "#7d7d7dff"
    #   )+
    #   ggplot2::geom_label(
    #     data = dormancy_df2,
    #     ggplot2::aes(x = as.Date(date), y = ylab_pos + 1.5*ylab_pos, label=event),
    #     size = 3, fill="#c0c0c0ff"
    #   )
  }
  posterior_plot = posterior_plot +
    ggplot2::geom_segment(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date),
                   xend = as.Date(Date),
                   y=0,
                   yend = ylab_pos),
      size = .5 #, linetype="dashed"
    )+
    ggplot2::geom_label(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date), y = ylab_pos, label=Name),
      size = 3,
      family = "Times New Roman"
    )
  posterior_plot + geom_hline(yintercept=0, color="black")
}


#' Plot clinical timeline + posterior times with MAP
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_timing_MAP = function(x)
{
  #clinical_timeline #= x$Input$Samples
  estimates = TOSCA:::get_inferred_parameters(x)

  endpoints = x$Input$Samples

  if (nrow(endpoints) > 2) endpoints = (endpoints %>% dplyr::arrange(as.Date(Date)))[2:3,]

  therapies = x$Input$Therapies

  # 1. time posterior plot
  timing_estimates = estimates %>%
    dplyr::select(starts_with('t_')) #%>%
  #apply(2, TOSCA:::convert_date_real, x=x) %>%
  #dplyr::as_tibble()
  #times = timing_estimates$variable %>% unique()

  if (x$Fit$model_info$dormancy) {
    dormancy_start = TOSCA:::convert_date_real(date = TOSCA:::get_start_therapy(x, class= "Chemotherapy inducing dormancy"), x=x)
    dormancy_end =  TOSCA:::convert_date_real(date = timing_estimates %>% pull(t_dormancy_end) %>% mean(), x=x)
    timing_estimates = timing_estimates %>% select(!c("t_mrca_tr",
                                                      colnames(timing_estimates)[grepl("t_cna_tr", colnames(timing_estimates))]))
  }

  timing_estimates = timing_estimates %>%
    apply(2, TOSCA:::convert_date_real, x=x) %>%
    dplyr::as_tibble()

  timing_estimates = TOSCA:::convert_timing_names(timing_estimates)

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

  # Compute medians per variable
  medians_df <- timing_estimates %>%
    group_by(variable) %>%
    summarise(median = median(value), .groups = "drop")

  # Reorder factor levels of 'variable' by descending median (convert Date -> numeric)
  timing_estimates <- timing_estimates %>%
    mutate(variable = factor(variable,
                             levels = medians_df$variable[order(-as.numeric(medians_df$median))]))

  # Also reorder medians_df to match
  medians_df <- medians_df %>%
    mutate(variable = factor(variable,
                             levels = levels(timing_estimates$variable)))

  # Timing Plot
  # timing_plot2 <- ggplot(timing_estimates) +
  #   geom_violin(aes(x = variable, y = value, color = variable, fill = variable),
  #               alpha = .5, scale = "width") +
  #   scale_color_manual(values = times_colors) +
  #   scale_fill_manual(values = times_colors) +
  #   geom_point(data = medians_df,
  #              aes(y = median, x = variable, color = variable),
  #              shape = 18, size = 2.5) +
  #   theme_bw() +
  #   labs(y = 'Date', x = '') +
  #   theme(
  #     axis.text.y = element_text(face = 'bold', size = 10),
  #     legend.position = "none"
  #   ) +
  #   coord_flip()
  all_dates <- c(
    as.Date(therapies$Start),
    as.Date(therapies$End),
    timing_estimates$value,
    as.Date(endpoints$Date)
  )
  xlims <- range(all_dates, na.rm = TRUE)
  xlims <- xlims + c(-10, 10)

  timing_plot2 <- ggplot(timing_estimates) +
    see::geom_violinhalf(aes(x = variable, y = value, fill = variable),
                color = "black",
                alpha = 1, scale = "width") +
    scale_color_manual(values = times_colors) +
    scale_fill_manual(values = times_colors) +
    # geom_point(data = medians_df,
    #            aes(y = median, x = variable, color = variable),
    #            shape = 18, size = 2.5) +
    theme_bw() +
    labs(y = 'Date', x = '') +
    theme(
      #axis.text.y = element_text(face = 'bold', size = 10),
      legend.position = "none"
    ) +
    # IMPORTANT: set the date limits with scale_y_date because 'value' is on y
    ylim(xlims[1],xlims[2])+
    #scale_y_date(limits = xlims, expand = c(0,0)) +
    coord_flip()

  timing_plot2 <- timing_plot2 +
    theme(
      # Clean theme: no grid, no border
      #panel.grid = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_blank(),
      # Hide x-axis labels (we keep only bottom one)
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      # Remove y axis line (keep only x at bottom plot)
      axis.line.x = element_blank(),
      axis.line.y = element_line(),
      legend.position = "none",
      text = element_text(family = "Times New Roman", size = 14, face = "plain"),
      plot.margin = unit(c(0,1,0,0), "cm") # top, right, bottom, left
    )

  # Clinical plot
  therapies = therapies %>% dplyr::mutate(Duration = as.Date(End)-as.Date(Start)) %>%
    dplyr::mutate(short=ifelse(Duration < 30, T, F)) %>% dplyr::mutate(Drug=Name)
  new_col = data.frame(Name = names(clinical_colors), colors = clinical_colors)
  therapies = dplyr::left_join(therapies,new_col, by="Name")

  dates= data.frame(x = seq.Date(from = xlims[1], to = xlims[2], length.out = 7))
  #data.frame(x = seq(xlims[1], xlims[2], by = 365*2))
  clinical_plot = ggplot() +
    #geom_line(yintercept = 0)+
    ggplot2::geom_rect(
    data = therapies, #%>% dplyr::filter(short == F),
    ggplot2::aes(xmin = as.Date(Start),
                 xmax = as.Date(End),
                 fill=Drug),
    ymin = 1,
    ymax = 1.1,
    #fill = c,
    colour = "white",
    size = 0.5,
    alpha = 1
  ) +
    scale_fill_manual(values = clinical_colors)+
    # ggplot2::geom_segment(
    #   data = therapies %>% dplyr::filter(short == T),
    #   ggplot2::aes(x=as.Date(Start),
    #                xend=as.Date(Start),
    #                y=0, yend=Inf, color = Name), alpha=1)  +
    # scale_color_manual(values = clinical_colors)+
    xlab("Date")+
    theme_minimal() +
    theme(
      # Remove everything related to y-axis
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),

      panel.grid.major.y = element_blank(), #
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(), #
      panel.grid.minor.x = element_blank(), #

      legend.position = "bottom",
      text = element_text(family = "Times New Roman"),
      plot.margin = unit(c(0,1,0,0), "cm"), # top, right, bottom, left
      legend.title = element_text(size = 10),  # legend title size
      legend.text  = element_text(size = 10)
    ) +
    # Force x-axis line at y=0
    scale_y_continuous(expand = c(0, 0))
    #geom_hline(yintercept = 0, color = "black")+ylim(0,2)+
  #clinical_plot


  if (x$Fit$model_info$dormancy) {

    dormancy_df = data.frame("event"= c("Dormancy"), "Start"=c(dormancy_start), "End"=c(dormancy_end))
    dormancy_df2 = data.frame("event"= c("Dormancy\nStart", "Dormancy\nEnd"), "date"=c(dormancy_start, dormancy_end))
    dormancy_df2$y = c(1.5, 1.7)

    clinical_plot = clinical_plot +
      ggplot2::geom_rect(
        data = dormancy_df,
        ggplot2:::aes(xmin = as.Date(Start), xmax = as.Date(End)),
        ymin = 1,ymax = 1.1, color = "white", fill="#c0c0c0ff", alpha = .5
      )
    # +
    #   ggplot2::geom_segment(
    #     data = dormancy_df2,
    #     ggplot2::aes(x = as.Date(date),
    #                  xend = as.Date(date),
    #                  y=1,
    #                  yend = 1.5),
    #     size = .5, linetype="dashed", color = "#7d7d7dff"
    #   )+
    #   ggplot2::geom_label(
    #     data = dormancy_df2,
    #     ggplot2::aes(x = as.Date(date), y = y, label=event),
    #     size = 3, fill="#c0c0c0ff"
    #   )
  }

  endpoints$shape = "1"
  medians_df$shape = "2"
  endpoints$color = c("S1", "S2")
  times_colors_df$variable = rownames(times_colors_df)
  medians_df = left_join(medians_df, times_colors_df)
  endpoints2 = rbind(endpoints, medians_df %>% dplyr::rename(Name=variable, Date = median) %>% dplyr::mutate(Date = as.character(Date)) ) %>%
    mutate(color = ifelse(color %in% c("S1", "S2"), color, Name))
  times_colors2 = c(times_colors, "S1"="black", "S2"="black") #"S1"="#a272bfff", "S2"="#683b81ff")
  #endpoints$color = c("S1", "S2")
  clinical_plot <- clinical_plot +
    #scale_x_date(limits = xlims, expand = c(0,0))+
    # theme(
    #   panel.grid = element_blank(),
    #   panel.border = element_blank(),
    #   axis.line.y = element_blank(),
    #   # collapse tick labels onto axis line
    #   axis.text.x = element_text(margin = margin(t = 0)),
    #   axis.ticks.length = unit(0, "pt")
    # )+
    geom_hline(yintercept = 1, color = "black")+ylim(.8,1.1)+
    #geom_vline(xintercept = xlims[1], color = "black")+
    geom_point(data = dates, aes(x=x, y=1))+
    geom_text(data = dates %>% mutate(d = as.character(x)) %>% rowwise() %>%
                mutate(d = paste0(strsplit(d, "-")[[1]][2], "-",strsplit(d, "-")[[1]][1])),
              aes(x=x, y=0.95, label = d),
              size=4,
              family = "Times New Roman")+
    ggplot2::geom_segment(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date),
                   xend = as.Date(Date),
                   y=1,
                   yend = .85),
      size = .5 #, linetype="dashed"
    )+
    ggplot2::geom_point(
      data = endpoints2,
      ggplot2::aes(x = as.Date(Date), y = 1, color = color, shape = shape, size=shape)
      #size = 3.7
    ) + scale_color_manual(values = times_colors2)+
    scale_shape_manual(values = c("1"=16, "2"=18))+
    scale_size_manual(values = c("1"=3.7, "2"=3))+
    ggplot2::geom_label(
      data = endpoints,
      ggplot2::aes(x = as.Date(Date), y = .85, label=Name),
      family = "Times New Roman",
      size = 4
    ) + xlim(xlims[1], xlims[2]) + guides(color = "none", shape = "none", size="none")

  #timing_estimates$value <- as.Date(timing_estimates$value)

  require(patchwork)
  combined <- timing_plot2 / clinical_plot + plot_layout(heights = c(1, .7))
  combined

}



