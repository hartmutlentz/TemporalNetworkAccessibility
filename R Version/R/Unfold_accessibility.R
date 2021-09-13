library(Matrix)
library(pracma)
library(ggplot2)

###### Functions ########
#### prepare data
#' Preparation of network data
#'
#' \code{prepare_data} Prepares the network data in a way that a list of
#' sparse matrices can be created
#'
#' @param df data frame with three columns: from, to, time
#'
#' @return a modified data frame
#' @export
#'
#' @examples
prepare_dataframe <- function(df) {
  if (min(df[, 3]) == 0) {
    df[, 3] <- df[, 3] + 1
  }
  names(df) <- c("from", "to", "time")
  unq <- unique(c(df$from, df$to))
  df$from <- match(df$from, unq)
  df$to <- match(df$to, unq)
  return(df)
}

##### Create Network function
##### Input dataframe with integers
##### Output: List of sparse matrices

#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
create_temporal_network <- function(df) {
  ntimesteps <- max(df$time)
  nfarms <- max(c(max(df$from), max(df$to)))
  ntransports <- nrow(df)
  temp_matrix <- Matrix(data = 0, nrow = nfarms, ncol = nfarms, sparse = TRUE)
  temp_net <- list(temp_matrix)
  for (i in 1:ntimesteps) {
    temp_net[[i]] <- temp_matrix
  }
  for (i in 1:ntransports) {
    temp_net[[df$time[i]]][df$from[i], df$to[i]] <- 1
  }
  return(temp_net)
}


#### Calculate cumulative paths
#' Title
#'
#' @param sm
#'
#' @return A dataframe with the time, the cumulative path and the shortest paths
#' @export
#'
#' @examples
cumul_path <- function(sm) {
  p <- sm[[1]]
  # create identity matrix
  d <- bandSparse(dim(p)[1], dim(p)[2], 0, list(rep(1, dim(p)[1] + 1)))
  p <- p + d
  cumu <- nnzero(p)
  for (i in seq_len(length(sm))) {
    print(i)
    p <- p + p %*% sm[[i]]
    cumu[i + 1] <- nnzero(p)
  }
  h <- pracma::gradient(cumu)
  dfout <- data.frame(
    t = seq_len(length(h)),
    h = h,
    c = cumu
  )
  return(dfout)
}

#' Title
#'
#' @param A dataframe with the time, the cumulative path and the shortest paths
#'
#' @return A ggplot figure
#' @export
#'
#' @examples
plot_tn <- function(df) {
  f <- max(df$c) / max(df$h)
  plt <- ggplot(df) +
    geom_line(aes(x = t, y = c), size = 1.5) +
    geom_line(aes(x = t, y = h * f), color = "red") +
    geom_point(aes(x = t, y = h * f), color = "red") +
    ylab("cumulative #paths") +
    scale_y_continuous(sec.axis =
                         sec_axis(~ . / f, name = "# shortest paths")) +
    theme_bw() +
    theme(legend.position = "none")
  return(plt)
}
