#' Plot the cumulative effects of lags
#'
#' Plots the cumulative effects over lags in a distributed lag model.
#'
#' @param model fitted `mgcv` model with a distributed lag term
#' @param facet_by which variable should be used to facet the plots
#' @param group_by which variable should be used to group the lines
#' @param term which term in the model is the distributed lag
#' @param group_labs label the groups? (Experimental)
#' @author David L Miller
#' @importFrom dplyr arrange group_by
#' @importFrom ggplot2 geom_line geom_hline scale_colour_viridis_c facet_wrap labs
#' @export
#' @examples
#' library(peach)
#' library(mgcv)
#' library(lubridate)
#'
#' mp <- process_nth_arrival(n=5)
#' mp <- mp[["M.persicae"]]
#'
#' # create the template
#' template <- mp[, c("Year", "Site")]
#' # want 1 Jan to 28 Feb
#' template$end <- ymd(paste0(mp$Year, "-02-28"))
#' template$start <- ymd(paste0(mp$Year, "-01-01"))
#' # note that there are NAs in these ranges, which process_temps
#' # will ignore leaving NA in those places
#' tt <- process_temps(template=template)
#'
#' # merge that into the observation data
#' mp <- merge(mp, tt)
#'
#' # fit a simple model
#' m <- gam(Arrive ~ te(temp, lag),
#'          control=list(keepData=TRUE),
#'          data=mp, method="REML",
#'          family=gaussian())
#'
#' plot_dlnm_cumulative(m)
plot_dlnm_cumulative <- function(model, facet_by="Site", group_by="Year",
                                 term=1, group_labs=FALSE){

  # exit if we didn't save the data
  if(is.null(model$data)){
    stop("Model was not run with `control=list(keepData=TRUE)`, no data to plot")
  }else{
    # determine the distributed lag term
    our_terms <- model$smooth[[term]]$term
    cov_name <- our_terms[1]
    lag_name <- our_terms[2]

    pred <- model$data
    # deal with NA values that were removed from the analysis
    if(!is.null(model$na.action)){
      pred <- pred[-model$na.action, ]
    }

    lags <- unique(pred[[lag_name]])[1,]

    pred2 <- pred[rep(1:nrow(pred), rep(length(lags), nrow(pred))), ]
    pred2[[lag_name]] <- as.numeric(1:length(lags))
    pred2$tmp2 <- NA

    for(ii in 1:nrow(pred2)){
      pred2$tmp2[ii] <- pred2[[cov_name]][ii, pred2[[lag_name]][ii]]
    }

    pred2[[cov_name]] <- pred2$tmp2


    pred2$p <- predict(model, newdata=pred2, type="iterms")[, term]

    xlim <- range(pred2[[lag_name]])
    ylim <- range(pred2$p)

    pred2 <- pred2 %>%
      group_by(.data[[facet_by]], .data[[group_by]]) %>%
      arrange(.data[["lag"]]) %>%
      mutate(cump = cumsum(.data[["p"]]))

    # make some group labels
    glabs <- pred2 %>%
      group_by(.data[[facet_by]], .data[[group_by]]) %>%
      summarize(mlag = max(.data[[lag_name]]),
                mp = .data[["cump"]][which.max(.data[[lag_name]])])

    res <- ggplot(pred2) +
      geom_line(aes(x=.data[[lag_name]], y=.data[["cump"]],
                    colour=.data[[group_by]], group=.data[[group_by]])) +
      geom_hline(yintercept=0, lty=2) +
      scale_colour_viridis_c() +
      facet_wrap(facet_by) +
      labs(x="Lag", y="Cumulative effect") +
      theme_minimal()

    if(group_labs){
      res <- res +
        geom_text(aes(x=.data[["mlag"]], y=.data[["mp"]], label=.data[[group_by]]),
                  data=glabs)

    }
    res
  }
}
