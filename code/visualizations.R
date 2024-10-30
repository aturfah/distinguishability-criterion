## Functions that implement the visualizations for the paper
## Date: 4/19/2024

library(latex2exp)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(tidyr)

source("code/PHM_algorithm.R")

#' Build dendrogram data from output of `PHM()`
#' 
#' @param phm_output Output of `PHM()`
.buildPHMDendrogramData <- function(phm_output, uniform_heights=F) {
  K <- length(phm_output)
  pmc <- phm_output[[K]]$pmc
  pmc_remains <- sapply(K:2, function(k) phm_output[[k]]$pmc)
  pmc_change <- sapply((K-1):1, function(k) phm_output[[k]]$pmc_change)
  height <- if (uniform_heights) {1:(K-1)} else {
    pmc / (pmc_remains)
  }
  merge_components <- t(sapply(K:2, function(k) phm_output[[k]]$merge_components))
  
  ## Figure out the components that merge together
  output <- data.frame()
  merge_tree <- lapply(1:K, function(k) k)
  height_tracker <- lapply(1:K, function(k) list(base=0, height=0))
  component_id_map <- 1:K
  base_height <- 0
  for (idx in 1:(K-1)) {
    ## Figure out which components merge
    mcs <- unname(merge_components[idx, ])
    hgt <- height[idx]
    
    ## Add this height to all components in height_tracker
    for (posn in 1:length(height_tracker)) height_tracker[[posn]]$height <- hgt # height_tracker[[posn]]$height + hgt
    
    ## For the components that were merged, add a vertical bar for them
    new_rows <- rbind(c(ID=component_id_map[mcs[1]], 
                        y=height_tracker[[mcs[1]]]$base,
                        yend=height_tracker[[mcs[1]]]$height,
                        pmc=pmc_remains[idx],
                        pmc_change=pmc_change[idx]),
                      c(ID=component_id_map[mcs[2]],
                        y=height_tracker[[mcs[2]]]$base,
                        yend=height_tracker[[mcs[2]]]$height,
                        pmc=pmc_remains[idx],
                        pmc_change=pmc_change[idx]))
    if (nrow(output) == 0)  {
      output <- data.frame(new_rows)
    } else {
      output <- rbind(new_rows, output)
    }
    
    ## Make note of the combined merge nodes so we have the merges stored somewhere
    merge_tree[[mcs[1]]] <- list(
      merge_tree[[mcs[1]]],
      merge_tree[[mcs[2]]]
    )
    merge_tree[[mcs[2]]] <- NULL
    
    ## Remove the component and re-map
    base_height <-  hgt
    height_tracker[[mcs[2]]] <- NULL
    height_tracker[[mcs[1]]]$base <- base_height
    
    
    component_id_map <- component_id_map[-mcs[2]]
  }
  
  ## Figure out where everything goes relative to base node
  order_x <- function(merge_res) {
    if (typeof(merge_res) == "list") {
      if (length(merge_res) == 1) {
        return(order_x(merge_res[[1]]))
      } else {
        return(c(order_x(merge_res[[1]]), order_x(merge_res[[2]])))
      }
    } else {
      return(merge_res)
    }
  }
  x_posns <- order_x(merge_tree)

  map_xposns <- function(vec) {
    sapply(vec, function(x) which(x_posns == x))
  }
  output <- output %>%
    mutate(x=ifelse(y==0, map_xposns(ID), NA))
  
  while(any(is.na(output))) {
    output <- output %>%
      left_join(output, by=c("y"="yend")) %>%
      dplyr::rename(x=x.x,
             ID=ID.x,
             pmc=pmc.x,
             pmc_change=pmc_change.x) %>%
      group_by(ID, y, yend, x, pmc, pmc_change) %>%
      summarize(x=ifelse(all(is.na(x)), mean(x.y), mean(x)),
                xend=x, .groups="keep") %>%
      ungroup() %>%
      arrange(-yend) %>%
      dplyr::select(-ends_with(".y"))
  }
  
  ## Add horizontal components
  horiz_comps <- output %>%
    group_by(yend, pmc, pmc_change) %>%
    summarize(
      xend=max(x),
      x=min(x),
      .groups="keep"
    ) %>%
    mutate(y=yend)
  output <- bind_rows(
    output, 
    horiz_comps
  )
  
  labels=output %>%
    filter(is.na(ID)) %>%
    mutate(xposn=(x+xend)/2,
           lab=ifelse(
             pmc_change < 1e-3,
             formatC(pmc_change, format = "e", digits = 2),
             round(pmc_change, 4)
           ))
  
  list(
    df=output,
    xlab=x_posns,
    labels=labels
  )
}

#' Generate the PHM dendrogram visualization
#' 
#' QQQ
#' 
#' @param phm_output Output from the `PHM()` function
#' @param colors Optional vector with colors for the mixture components
#' @param suppress_labels Suppress text boxes including Pmc reduction
#' @param suppress_axis Suppress axis labels with colors
plotPHMDendrogram <- function(phm_output, colors=NULL, suppress_labels=F, suppress_axis=F, uniform_heights=F) {
  pmc_dendro_data <- .buildPHMDendrogramData(phm_output, uniform_heights=uniform_heights)
  K <- length(pmc_dendro_data$xlab)
  if (is.null(colors) & !suppress_axis) colors <- brewer.pal(K, "Set1")

  if (suppress_axis) {
    AXIS_X <- element_text(size=6)
  } else {
    AXIS_X <- element_text(color=unname(colors[pmc_dendro_data$xlab]),
                           size=10)
  }
  
  offset <- 1
  if (uniform_heights) {
    offset <- 0
  }

  plt <- ggplot(pmc_dendro_data$df, aes(x=x, y=y+offset, xend=xend, yend=yend+offset)) +
    geom_segment() +
    # xlab("Mixture Component ID") +
    xlab("") +
    ylab("") +
    # ylab(TeX("$P_{MC} /$ Remaining $P_{MC}$")) +
    theme_bw() + theme(text=element_text(size=8),
                       axis.text.x=AXIS_X,
                       axis.text.y=element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid=element_blank(),
                       panel.spacing=unit(0, "lines"),
                       panel.border=element_blank()
    )
  
  if (!uniform_heights) {
    plt <- plt + scale_y_log10(expand=expansion(mult=c(0, 0.05)))
  } else {
    plt <- plt + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
  }

  if (!suppress_labels) {
    plt <- plt + geom_label(data=pmc_dendro_data$labels,
              aes(x=xposn, y=y+1, label=lab),
              size=2,
              label.size=0.15,
              label.padding = unit(0.15, "lines"),
              label.r = unit(0.1, "lines"))
  }
  if (!suppress_axis) {
    plt <- plt + scale_x_continuous(breaks=1:K,
                                    labels=rep("\U25A0", K))
  } else {
    plt <- plt + scale_x_continuous(breaks=1:K,
                                    labels=pmc_dendro_data$xlab)
  }
  plt
}

#' Plot $\Delta \P_{\rm{mc}}$ matrix
#' 
#' @param temp doot
#' 
#' @return pew
plotPmcMatrix <- function(phm_output, k=NULL, colors=NULL,
                          threshold=1e-3, threshold_replace="< 0.001", fmt_func=NULL,
                          include_pmc_title=T) {
  
  if (is.null(k)) k <- phm_output[[length(phm_output)]]$clusters
  if (is.null(colors)) colors <- brewer.pal(k, "Set1")
  if (is.null(fmt_func)) fmt_func<- function(x) round(x, 4)

  
  pmc_matrix <- phm_output[[k]]$pmc_matrix
  pmc <- sum(pmc_matrix)
  title_obj <- ggtitle(TeX(paste0("$P_{mc} =", fmt_func(pmc), "$")))

  colnames(pmc_matrix) <- paste0("V", 1:k)
  
  pmc_matrix_long <- (2 * pmc_matrix) %>%
    # Data wrangling
    as_tibble(.name_repair="unique") %>%
    rowid_to_column(var="X") %>%
    gather(key="Y", value="Z", -1) %>%
    mutate(Y=as.numeric(gsub("V","",Y)))  

  htmp <- pmc_matrix_long %>%
    mutate(Z.mod = ifelse(
        Z < threshold,
        threshold_replace,
        fmt_func(Z)
      ),
      Z.mod = ifelse(X == Y, "--", Z.mod)) %>%
    filter(X <= Y) %>%
    # Viz
    ggplot(aes(X, Y, fill=Z)) + 
    geom_tile(color="black") +
    geom_text(aes(label=Z.mod, fill=NULL), 
              size=2) +
    scale_fill_gradient(limits=c(0, pmc),
                        low="white", high="red",
                        aesthetics="fill") +
    ## Format the axis labels
    scale_x_continuous(breaks=1:k,
                       minor_breaks = NULL,
                       labels=rep("\U25A0", k)) +
    scale_y_continuous(breaks=1:k,
                       minor_breaks = NULL,
                       labels=rep("\U25A0", k)) +
    theme_bw() +
    coord_flip() +
    xlab("") + ylab("") +
    theme(text=element_text(size=8),
          legend.position="NULL",
          panel.grid=element_blank(),
          axis.text.x=element_text(color=colors),
          axis.text.y=element_text(color=colors))
  
  if (include_pmc_title == T) 
    htmp <- htmp + title_obj
  
  return(htmp)
}


#' Generate the distruvt plot from the hierarchical merge
#' 
#' QQQ
#' 
#' @param phm_output Output from the `PHM() function`
#' @param k Number of clusters for which to generate the distruct plot
#' @param labels Ground truth class labels for the observations (ordered factor vector)
#' @param colors Optinal vector with colors for the mixture components
#' @param include_title Whether to include the default title with Pmc and k values
plotPHMDistruct <- function(phm_output, k=length(phm_output),  
                         labels, colors=NULL, include_title=FALSE) {

  ## Validation
  if (is.null(colors)) colors <- brewer.pal(k, "Set1")
  post_mat <- phm_output[[k]]$posterior_matrix
  if (is.null(post_mat)) stop("Does phm_output have a posterior matrix computed?")
  
  
  labels_order=levels(labels)

  ## Format the title
  pmc_k <- phm_output[[k]]$pmc
  display_pmc <- ifelse(
    pmc_k < 1e-3,
    formatC(pmc_k, format = "e", digits = 2),
    round(pmc_k, 4)
  )
  
  ## Define break points for the distruct plot
  group_counts <- table(labels)
  label_positions <- sapply(1:length(group_counts), function(idx) {
    ifelse(idx != 1, sum(group_counts[1:(idx-1)]), 0) + group_counts[idx] / 2
  })

  ## Form data for plotting
  df <- data.frame(labels=labels, posterior=post_mat)
  
  group_columns_df <- df %>%
    rowwise() %>%
    mutate(max_column = which.max(c_across(where(is.numeric)))) %>%
    ungroup() %>%
    dplyr::select(-starts_with("posterior")) %>%
    group_by(labels, max_column) %>%
    summarize(count=n()) %>%
    group_by(labels) %>%
    filter(count == max(count))
  
  group_columns <- pull(group_columns_df, "max_column")
  names(group_columns) <- pull(group_columns_df, "labels")
  
  plt_data <- df %>%
    rowwise() %>%
    mutate(sort_val = c_across(where(is.numeric))[group_columns[labels]]) %>%
    ungroup() %>%
    arrange(labels, -sort_val) %>%
    mutate(row=1:n()) %>%
    pivot_longer(cols=starts_with("posterior")) %>%
    arrange(row, value)
    
  # plt_data <- df %>%
  #   arrange(labels, -across(paste0("posterior.", column_sort_order))) %>%
  #   mutate(row=1:n()) %>%
  #   pivot_longer(cols=starts_with("posterior")) 
  
  ## Positions for vertical lines
  vert_lines <- plt_data %>%
    group_by(labels) %>%
    summarize(minX=min(row),
              maxX=max(row)) %>%
    mutate(minX=minX-0.5, maxX=maxX+0.5,
           y=0, yend=1) %>%
    pivot_longer(cols=c(minX, maxX))
  
  ## Distruct plot
  plt <- plt_data %>%
    ggplot(aes(x=row, y=value, fill=name)) +
    geom_bar(position="stack", stat="identity") +
    xlab("") + ylab("") +
    geom_segment(data=vert_lines,
                 aes(x=value, xend=value, y=y, yend=yend, fill=NULL),
                 linetype="dashed", linewidth=0.3) +
    scale_fill_manual(values=colors) +
    scale_x_continuous(expand=c(0, 0), 
                       labels=labels_order,
                       breaks=label_positions
                       , minor_breaks = NULL) +
    scale_y_continuous(expand=c(0, 0),
                       breaks=NULL, minor_breaks = NULL, n.breaks=0) +
    theme_bw() + theme(title=element_text(size=6),
                       axis.text.x=element_text(hjust=0.5, size=6),
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank(),
                       panel.grid=element_blank(),
                       legend.position="none",
                       plot.margin = unit(c(0.01, 0.03, 0, 0.01), 
                                          "inches"))
  
  if (include_title == TRUE) {
    plt <- plt + ggtitle(TeX(paste0("K = ", k, "; ", "$P_{MC} =$", " ", display_pmc)))
  }
  
  return(plt)
}



# group_counts <- c(100, 300, 200)
# distbn_params_list <- list(
#   list(
#     prob=group_counts[1]/sum(group_counts),
#     mean=matrix(-3, nrow=1, ncol=1),
#     var=array(1, dim=c(1, 1, 1))
#   ),
#   list(
#     prob=group_counts[2]/sum(group_counts),
#     mean=matrix(0, nrow=1, ncol=1),
#     var=array(1, dim=c(1, 1, 1))
#   ),
#   list(
#     prob=group_counts[3]/sum(group_counts),
#     mean=matrix(5, nrow=1, ncol=1),
#     var=array(2, dim=c(1, 1, 1))
#   )
# )
# 
# data <- rbind(
#   .sampleMixture(list(distbn_params_list[[1]]), group_counts[1]),
#   .sampleMixture(list(distbn_params_list[[2]]), group_counts[2]),
#   .sampleMixture(list(distbn_params_list[[3]]), group_counts[3])
# )
# 
# labels <- factor(
#   c(rep("a", group_counts[1]), 
#     rep("b", group_counts[2]), 
#     rep("c", group_counts[3])),
#   levels=c("a", "c", "b"),
#   ordered=T
# )
# 
# # post_vals <- .computePosteriorProbMatrix(distbn_params_list, data)
# res_mclust <- Mclust(data)
# 
# phm_output <- PHM(res_mclust, data=data, mc_est=F)
# 
# plotPHMDendrogram(phm_output)
# plotPmcMatrix(phm_output, include_pmc_title = F)
# plotPHMDistruct(phm_output, labels=labels, colors=c("red", "cyan", "magenta"))
