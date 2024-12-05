
#' @title Plot smooth terms for mgcv::gam objects
#' @param model fitted model object
#' @param s name of predictor to plot smooth term for
#' @param scale if "response", inverse link function will be applied
#' @param rug show data rug at bottom of plot?
#' @param ... additional arguments to geom_ribbon. Providing args (e.g. fill) will reset defaults for other args (alpha, color, linetype).
#'
#' @return a ggplot object
#' @export
plot_smoother <- function(model, s, scale = c("link", "response"), rug = TRUE, ...) {
  
  require("ggplot2")
  
  ## Obtain values from plot.gam function
  plot_out <- plot(model, select = 0, n1 = 200)
  
  ## Match smoother name to output from plot.gam
  s_names <- vapply(plot_out, function(x) ifelse(is.null(x$xlab), "", x$xlab), FUN.VALUE = "character")
  s_num <- match(s, s_names)
  if(is.na(s_num)) stop(paste(s, "is not a 1D smooth predictor in provided model"))
  
  ## Obtain and store plot values from plot.gam
  plot_data <- data.frame(
    x = seq(plot_out[[s_num]]$xlim[1], plot_out[[s_num]]$xlim[2], 
            length.out = nrow(plot_out[[s_num]]$fit)),
    fit = plot_out[[s_num]]$fit,
    se = plot_out[[s_num]]$se
  )
  
  ## Upper and lower bounds
  plot_data$lower <- plot_data$fit - plot_data$se
  plot_data$upper <- plot_data$fit + plot_data$se
  
  ## Apply transformation to response scale if necessary
  scale <- match.arg(scale)
  
  if(ncol(model$model) > 2 & scale == "response" & model$family$link != "identity") {
    warning("scale = 'response' should not be used for models with multiple predictors")
  }
  
  plot_data$fit   <- switch(scale, link = plot_data$fit, 
                            response = model$family$linkinv(plot_data$fit))
  plot_data$lower <- switch(scale, link = plot_data$lower, 
                            response = model$family$linkinv(plot_data$lower))
  plot_data$upper <- switch(scale, link = plot_data$upper, 
                            response = model$family$linkinv(plot_data$upper))
  
  ## Return ggplot
  ggplot(plot_data, aes(x, fit)) +
    {
      if (length(list(...)) == 0) {
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", alpha = 0.5, 
                    color = "black", linetype = "dashed")
      } else {
        geom_ribbon(aes(ymin = lower, ymax = upper), ...)
      }
    } +
    geom_line() +
    labs(x = plot_out[[s_num]]$xlab, y = plot_out[[s_num]]$ylab) + 
    {
      if(rug) {
        geom_rug(aes(x = model$model[[s]]), data = model$model, 
                 inherit.aes = FALSE)
      }
    }
  
}

#' @title Plot 2D smooth terms for mgcv::gam objects
#' @param model fitted model object
#' @param s label corresponding to smooth term (e.g. "s(x, y, k = ..., ...)") to plot
#' @param se Logical - should a seprate plot for smoother standard error be returned?
#' @param rug Logical - should points corresponding to observed data be plotted?
#' @param crop_ebs Logical - should smoother be restricted to EBS survey region?
#'
#' @return A ggplot2 object
#' @export
plot_smoother_2D <- function(model, s, se = TRUE, rug = TRUE, crop_ebs = TRUE) {
  
  require("ggplot2")
  
  ## Obtain values from plot.gam function
  plot_out <- plot(model, select = 0, n2 = 250, rug = rug)
  
  ## Match smoother name to output from plot.gam
  s_names <- attr(terms(formula(model)), "term.labels")
  s_num <- match(s, s_names)
  if(is.na(s_num)) stop(paste(s, "is not a smooth predictor in provided model"))
  if(length(plot_out[[s_num]][["y"]]) == 0) stop(paste(s, "is not a 2-dimensional smooth predictor in provided model"))
  
  ## Obtain and store plot values from plot.gam
  plot_data <- data.frame(
    x = rep(plot_out[[s_num]]$x, length(plot_out[[s_num]]$y)),
    y = rep(plot_out[[s_num]]$y, each = length(plot_out[[s_num]]$x)),
    fit = plot_out[[s_num]]$fit,
    se = plot_out[[s_num]]$se
  )
  
  ## Crop to EBS survey region
  if (crop_ebs) {
    
    require("sf")
    
    plot_data$id <- 1:nrow(plot_data)
    
    ebs_boundary <- get_ebs_shapefile(region = "EBS", type = "boundary")
    
    plot_data_sf <- st_as_sf(
      mutate(plot_data, x = x*1000, y = y*1000), 
      coords = c("x", "y"), crs = st_crs(ebs_boundary)
    )
    
    plot_data_sf <- st_filter(plot_data_sf, ebs_boundary)
    
    plot_data <- plot_data |> filter(id %in% plot_data_sf$id) |> dplyr::select(-id)
    
    ebs_boundary <- fortify(as_Spatial(ebs_boundary))
    
  }
  
  ## Rug points
  rug_data <- unique(plot_out[[s_num]]$raw)

  fit_plot <- plot_data |> 
    ggplot(aes(x, y, z = fit)) + 
    geom_contour_filled(color = "black", bins = 15) + 
    guides(fill = guide_colorsteps(barwidth = unit(0.5, "npc"), label.hjust = 1, label.vjust = 1)) + 
    theme(
      legend.position = "bottom", 
      legend.text = element_text(angle = 90), 
      legend.title = element_blank()
    ) + 
    scale_fill_viridis_d(option = "magma") + 
    labs(
      x = plot_out[[s_num]]$xlab, 
      y = plot_out[[s_num]]$ylab, 
      title = plot_out[[s_num]]$main
    ) + {
      if(rug) geom_point(aes(x, y), data = rug_data, inherit.aes = FALSE, shape = 21, fill = NA, color = "white")
    } + {
      if(crop_ebs) {
        geom_polygon(
          aes(x = long/1000, y = lat/1000, group = group),
          data = ebs_boundary, color = "black", fill = NA, 
          linewidth = 0.8, inherit.aes = FALSE
        )
      }
    }
  
  
  if (se) {
    
    se_plot <- plot_data |> 
      ggplot(aes(x, y, z = se)) + 
      geom_contour_filled(color = "black", bins = 15) + 
      guides(fill = guide_colorsteps(barwidth = unit(0.5, "npc"), label.hjust = 1, label.vjust = 1)) + 
      theme(
        legend.position = "bottom", 
        legend.text = element_text(angle = 90), 
        legend.title = element_blank()
      ) + 
      scale_fill_viridis_d(option = "magma") + 
      labs(
        x = plot_out[[s_num]]$xlab, 
        y = plot_out[[s_num]]$ylab, 
        title = paste0("SE(", plot_out[[s_num]]$main, ")")
      ) + {
        if(crop_ebs) {
          geom_polygon(
            aes(x = long/1000, y = lat/1000, group = group),
            data = ebs_boundary, color = "black", fill = NA, 
            linewidth = 0.8, inherit.aes = FALSE
          )
        }
      }
    
    fit_plot <- cowplot::plot_grid(fit_plot, se_plot)
    
  }
  
  fit_plot
}
