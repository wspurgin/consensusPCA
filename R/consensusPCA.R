build_consensus <- function(data_set, consensus_size) {
  consensus_set <- list()
  sample_size <- nrow(data_set) / consensus_size
  for(i in 1:consensus_size) {
    d_indices <- sample(1:nrow(data_set), sample_size, replace = TRUE)
    d         <- data_set[d_indices, ]
    attribute <- paste("sample", i)
    consensus_set[[ attribute ]] <- prcomp(d)
  }
  consensus_set
}

reconcile <- function(consensus) {
  n     <- length(consensus)
  total <- consensus[[1]]
  total$x <- NULL
  class(total) <- "prcomp"
  for(i in 2:n) {
    curr           <- consensus[[i]]
    total$rotation <- (abs(total$rotation) + abs(curr$rotation))
    total$sdev     <- (total$sdev     + curr$sdev)
    total$center   <- (total$center   + curr$center)
  }

  total$rotation <- total$rotation / n
  total$sdev     <- total$sdev / n
  total$center   <- total$center / n

  # Orthonormalize the resultant, averaged rotation.
  r_names <- rownames(total$rotation)
  c_names <- colnames(total$rotation)
  total$rotation <- qr.Q(qr(total$rotation))
  rownames(total$rotation) <- r_names
  colnames(total$rotation) <- c_names
  total
}

cosine_difference <- function(a,b) {
  dot <- (a %*% b)
  norm_a <- norm(a, type = "2")
  norm_b <- norm(b, type = "2")
  inverse_theta <- dot / (norm_a * norm_b)
  if (!(inverse_theta <= 1 && inverse_theta >= -1)) {
    # Sometimes, especially when two vectors are the same, we get plagued by
    # roundoff putting our answer just out of the defined range of acos. So
    # we'll move it if necessary.
    if (inverse_theta > 1) {
      inverse_theta <- ifelse(abs(inverse_theta - 1) < 1e-6, 1, inverse_theta)
    } else {
      inverse_theta <- ifelse(abs(inverse_theta - -1) < 1e-6, -1, inverse_theta)
    }
  }
  acos( inverse_theta )
}

cosine_difference_degrees <- function(a,b) {
  cosine_difference(a, b) * (360 / 2 / pi)
}

# PCA Stability
# For each parital PCA, we'll produce a distance matrix by PC of
# cosine difference between the same PC on all other partial PCA
stabilization_analysis <- function(consensus_set) {
  consensus_size <- length(consensus_set)
  num_pc <- length(consensus_set[[1]]$center) #TODO error checking if 0-length consensus set
  pc_stability <- list()
  k <- consensus_size * consensus_size
  for( l in 1:num_pc) {
    name <- paste("PC", l, sep = "")
    pc_res <- list()
    pc_res$cosine_difference <- matrix(data = rep(0, k), nrow = consensus_size, ncol = consensus_size)
    for (i in 1:consensus_size) {
      partial_pca <- consensus_set[[i]]
      for (j in 1:consensus_size) {
        other <- consensus_set[[j]]
        pc_res$cosine_difference[j, i] <- cosine_difference_degrees( abs(partial_pca$rotation[,l]), abs(other$rotation[,l]) )
      }
    }
    pc_res$absolute_difference <- sum(pc_res$cosine_difference)
    pc_res$mean_difference <- mean(pc_res$cosine_difference)
    pc_res$stdev <- sd(pc_res$cosine_difference)
    pc_stability[[name]] <- pc_res
  }
  pc_stability
}
