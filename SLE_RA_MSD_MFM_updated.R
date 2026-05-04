############################################################
## SLE + RA joint clustering using MSD-MFM-ECM
## Updated from RA-only workflow
############################################################

## =========================================================
## 0. Load and preprocess data
## =========================================================

## Change this path if needed
raw_data <- read.csv("C:/Users/wujrt/Desktop/pyne_prject/pyne_prject/genes.csv")

## Original workflow: transpose, set gene names, scale genes
raw_data <- t(raw_data)
raw_data <- raw_data[, 1:9671]
colnames(raw_data) <- raw_data[1, ]
raw_data <- raw_data[-1, ]

X <- apply(raw_data[, 2:9671], 2, as.numeric)
X <- scale(X, center = TRUE, scale = TRUE)
Y <- raw_data[, 1]

## Split disease groups
SLE_X <- X[Y == "SLE", , drop = FALSE]
RA_X  <- X[Y == "RA",  , drop = FALSE]

## Combine SLE and RA for joint clustering
X_both <- rbind(SLE_X, RA_X)

disease_label <- c(
  rep("SLE", nrow(SLE_X)),
  rep("RA",  nrow(RA_X))
)

## Optional sample names
rownames(X_both) <- paste0(disease_label, "_", seq_len(nrow(X_both)))

rm(list = setdiff(ls(), c("X_both", "disease_label")))


############################################################
## 1. MSD-MFM-ECM model fitting function
############################################################

fit_msd_mfm_ecm <- function(
    X, gene_cluster,
    K = 3, q = 3,
    max_iter = 100, tol = 1e-5,
    nstart_kmeans = 20,
    verbose = TRUE,
    seed = 1,
    psi_floor = 1e-3,
    ridge = 1e-6
) {
  set.seed(seed)
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  ## ---- check + re-encode gene_cluster to 1..G contiguous ----
  gene_cluster <- as.integer(gene_cluster)
  if (length(gene_cluster) != p) stop("gene_cluster length must equal ncol(X).")
  if (any(is.na(gene_cluster))) stop("gene_cluster contains NA.")
  if (any(gene_cluster < 1)) stop("gene_cluster must be in 1..G with no 0.")
  gene_cluster <- as.integer(factor(gene_cluster))
  G <- max(gene_cluster)
  
  module_idx <- split(seq_len(p), gene_cluster)
  
  ## ---------- helpers ----------
  logsumexp <- function(v) {
    m <- max(v)
    m + log(sum(exp(v - m)))
  }
  
  module_means_from_gene_matrix <- function(X, gene_cluster, G) {
    Xt <- t(X)  ## p x n
    sums <- rowsum(Xt, group = gene_cluster, reorder = FALSE)  ## G x n
    counts <- tabulate(gene_cluster, nbins = G)
    t(sums / matrix(counts, nrow = G, ncol = ncol(sums)))
  }
  
  compute_B <- function(U, s, psi, gene_cluster, G, q) {
    wj <- (s^2) / psi
    sum_w_by_g <- as.numeric(tapply(wj, gene_cluster, sum))
    B <- matrix(0, q, q)
    for (g in 1:G) {
      ug <- matrix(U[g, ], ncol = 1)
      B <- B + sum_w_by_g[g] * (ug %*% t(ug))
    }
    B
  }
  
  ## Solve KKT for maximizing -1/2 u^T A u + b^T u subject to ||u|| = 1
  solve_unit_norm <- function(A, b, ridge_local = 1e-8) {
    A <- (A + t(A)) / 2
    A <- A + ridge_local * diag(nrow(A))
    
    eg <- eigen(A, symmetric = TRUE)
    Q <- eg$vectors
    d <- eg$values
    bt <- as.vector(t(Q) %*% b)
    
    f <- function(lambda) sum((bt / (d + lambda))^2) - 1
    lower <- max(1e-12, -min(d) + 1e-12)
    
    if (f(lower) < 0) {
      u <- Q %*% (bt / (d + lower))
      u <- as.vector(u / sqrt(sum(u^2) + 1e-12))
      return(u)
    }
    
    upper <- lower
    for (tt in 1:60) {
      upper <- upper * 2
      if (f(upper) < 0) break
    }
    lambda_star <- uniroot(f, lower = lower, upper = upper, tol = 1e-10)$root
    u <- Q %*% (bt / (d + lambda_star))
    as.vector(u / sqrt(sum(u^2) + 1e-12))
  }
  
  compute_loglik <- function(pi_k, delta, alpha, U, s, psi) {
    B <- compute_B(U, s, psi, gene_cluster, G, q)
    A <- diag(q) + B
    cholA <- chol(A)
    V <- chol2inv(cholA)
    
    logdetA <- 2 * sum(log(diag(cholA)))
    logdetSigma <- sum(log(psi)) + logdetA
    const <- -0.5 * (p * log(2 * pi) + logdetSigma)
    
    invpsi <- 1 / psi
    w_s_over_psi <- s * invpsi
    
    logp_nk <- matrix(0, n, K)
    for (k in 1:K) {
      mu_k <- delta + alpha[k, gene_cluster]
      E <- sweep(X, 2, mu_k, "-")
      
      rPsiInvr <- rowSums((E^2) * matrix(invpsi, n, p, byrow = TRUE))
      
      sum_rg <- matrix(0, n, G)
      for (g in 1:G) {
        idx <- module_idx[[g]]
        sum_rg[, g] <- E[, idx, drop = FALSE] %*% w_s_over_psi[idx]
      }
      
      T_iq <- sum_rg %*% U
      M <- t(V %*% t(T_iq))
      tVt <- rowSums(T_iq * M)
      quad <- rPsiInvr - tVt
      
      logp_nk[, k] <- log(pi_k[k] + 1e-12) + const - 0.5 * quad
    }
    sum(apply(logp_nk, 1, logsumexp))
  }
  
  ## ---------- Initialization ----------
  Xm0 <- module_means_from_gene_matrix(X, gene_cluster, G)
  Xm0s <- scale(Xm0)
  
  n_distinct <- nrow(unique(round(Xm0s, 10)))
  if (n_distinct < 2) stop("Too few distinct samples in module-mean space; cannot initialize.")
  K_eff <- min(K, n_distinct)
  if (K_eff < K && verbose) {
    cat(sprintf("Init: requested K=%d but only %d distinct points; using K=%d\n", K, n_distinct, K_eff))
  }
  K <- K_eff
  
  km_init <- kmeans(Xm0s, centers = K, nstart = nstart_kmeans)
  c_init <- km_init$cluster
  
  r <- matrix(0, n, K)
  r[cbind(seq_len(n), c_init)] <- 1
  pi_k <- colMeans(r)
  
  delta <- colMeans(X)
  
  alpha <- matrix(0, nrow = K, ncol = G)
  Xc0 <- sweep(X, 2, delta, "-")
  Xm_c0 <- module_means_from_gene_matrix(Xc0, gene_cluster, G)
  for (k in 1:K) {
    w <- r[, k]
    sw <- sum(w)
    if (sw > 1e-10) alpha[k, ] <- colSums(Xm_c0 * w) / sw
  }
  
  U <- matrix(0, nrow = G, ncol = q)
  for (g in 1:G) {
    v <- rnorm(q)
    U[g, ] <- v / sqrt(sum(v^2) + 1e-12)
  }
  
  s <- rep(1, p)
  
  alpha_bar <- r %*% alpha
  mu0 <- matrix(delta, n, p, byrow = TRUE) + alpha_bar[, gene_cluster, drop = FALSE]
  resid0 <- X - mu0
  psi <- pmax(colMeans(resid0^2), psi_floor)
  
  loglik_trace <- numeric(max_iter)
  converged <- FALSE
  best_ll <- -Inf
  best_par <- NULL
  
  ## ---------- ECM loop ----------
  for (iter in 1:max_iter) {
    
    ## E-step
    B <- compute_B(U, s, psi, gene_cluster, G, q)
    A <- diag(q) + B
    cholA <- chol(A)
    V <- chol2inv(cholA)
    
    logdetA <- 2 * sum(log(diag(cholA)))
    logdetSigma <- sum(log(psi)) + logdetA
    const <- -0.5 * (p * log(2 * pi) + logdetSigma)
    
    invpsi <- 1 / psi
    w_s_over_psi <- s * invpsi
    
    logp_nk <- matrix(0, n, K)
    M_ik <- array(0, dim = c(n, q, K))
    a_ikg <- array(0, dim = c(n, G, K))
    
    for (k in 1:K) {
      mu_k <- delta + alpha[k, gene_cluster]
      E <- sweep(X, 2, mu_k, "-")
      
      rPsiInvr <- rowSums((E^2) * matrix(invpsi, n, p, byrow = TRUE))
      
      sum_rg <- matrix(0, n, G)
      for (g in 1:G) {
        idx <- module_idx[[g]]
        sum_rg[, g] <- E[, idx, drop = FALSE] %*% w_s_over_psi[idx]
      }
      
      T_iq <- sum_rg %*% U
      M <- t(V %*% t(T_iq))
      M_ik[, , k] <- M
      a_ikg[, , k] <- M %*% t(U)
      
      tVt <- rowSums(T_iq * M)
      quad <- rPsiInvr - tVt
      
      logp_nk[, k] <- log(pi_k[k] + 1e-12) + const - 0.5 * quad
    }
    
    for (i in 1:n) {
      lse <- logsumexp(logp_nk[i, ])
      r[i, ] <- exp(logp_nk[i, ] - lse)
    }
    
    ## CM-step 1: pi
    pi_k <- colMeans(r)
    
    ## CM-step 2: alpha
    mean_s_g <- as.numeric(tapply(s, gene_cluster, mean))
    Xc <- sweep(X, 2, delta, "-")
    Xm <- module_means_from_gene_matrix(Xc, gene_cluster, G)
    
    for (k in 1:K) {
      w <- r[, k]
      sw <- sum(w)
      if (sw < 1e-12) next
      adj <- Xm - sweep(a_ikg[, , k], 2, mean_s_g, "*")
      alpha[k, ] <- colSums(adj * w) / sw
    }
    
    ## CM-step 3: delta
    alpha_bar <- r %*% alpha
    Abar_ng <- matrix(0, n, G)
    for (k in 1:K) Abar_ng <- Abar_ng + r[, k] * a_ikg[, , k]
    
    AlphaGene <- alpha_bar[, gene_cluster, drop = FALSE]
    FactorGene <- Abar_ng[, gene_cluster, drop = FALSE] * matrix(s, n, p, byrow = TRUE)
    delta <- colMeans(X - AlphaGene - FactorGene)
    
    ## CM-step 4: s
    uV_u <- numeric(G)
    for (g in 1:G) uV_u[g] <- as.numeric(t(U[g, ]) %*% V %*% U[g, ])
    
    for (g in 1:G) {
      idx <- module_idx[[g]]
      if (length(idx) == 0) next
      
      a_mat <- matrix(0, n, K)
      for (k in 1:K) a_mat[, k] <- a_ikg[, g, k]
      
      b_mat <- uV_u[g] + a_mat^2
      denom_g <- sum(r * b_mat)
      if (denom_g < 1e-12) next
      
      for (j in idx) {
        num_j <- 0
        for (k in 1:K) {
          mu_kj <- delta[j] + alpha[k, gene_cluster[j]]
          e <- X[, j] - mu_kj
          num_j <- num_j + sum(r[, k] * (e * a_mat[, k]))
        }
        s[j] <- num_j / denom_g
      }
    }
    
    ## CM-step 5: U
    invpsi <- 1 / psi
    U_new <- U
    
    for (g in 1:G) {
      idx <- module_idx[[g]]
      if (length(idx) == 0) next
      
      w_g <- sum((s[idx]^2) * invpsi[idx])
      if (w_g < 1e-12) next
      
      A_g <- matrix(0, q, q)
      b_g <- rep(0, q)
      
      for (k in 1:K) {
        Mk <- M_ik[, , k]
        
        ebar <- rep(0, n)
        for (j in idx) {
          mu_kj <- delta[j] + alpha[k, gene_cluster[j]]
          ebar <- ebar + (s[j] * invpsi[j]) * (X[, j] - mu_kj)
        }
        
        W <- r[, k] * ebar
        b_g <- b_g + as.vector(t(Mk) %*% W)
        
        sw <- sum(r[, k])
        A_g <- A_g + w_g * (sw * V + t(Mk) %*% (Mk * r[, k]))
      }
      
      U_new[g, ] <- solve_unit_norm(A_g + ridge * diag(q), b_g, ridge_local = ridge)
    }
    U <- U_new
    
    ## CM-step 6: psi
    psi_new <- numeric(p)
    
    uV_u <- numeric(G)
    for (g in 1:G) uV_u[g] <- as.numeric(t(U[g, ]) %*% V %*% U[g, ])
    
    for (g in 1:G) {
      idx <- module_idx[[g]]
      if (length(idx) == 0) next
      
      ug <- U[g, ]
      a_mat <- matrix(0, n, K)
      for (k in 1:K) a_mat[, k] <- as.vector(M_ik[, , k] %*% ug)
      
      b_mat <- uV_u[g] + a_mat^2
      
      for (j in idx) {
        sj <- s[j]
        acc <- 0
        for (k in 1:K) {
          mu_kj <- delta[j] + alpha[k, gene_cluster[j]]
          e <- X[, j] - mu_kj
          acc <- acc + sum(r[, k] * (e^2 - 2 * sj * e * a_mat[, k] + (sj^2) * b_mat[, k]))
        }
        psi_new[j] <- acc / n
      }
    }
    psi <- pmax(psi_new, psi_floor)
    
    ## Observed log-likelihood after full update
    ll_new <- compute_loglik(pi_k, delta, alpha, U, s, psi)
    loglik_trace[iter] <- ll_new
    if (verbose) cat(sprintf("Iter %3d  logLik = %.3f\n", iter, ll_new))
    
    if (ll_new > best_ll) {
      best_ll <- ll_new
      best_par <- list(pi = pi_k, delta = delta, alpha = alpha, U = U, s = s, psi = psi, r = r)
    }
    
    if (iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) <
        tol * (1 + abs(loglik_trace[iter - 1]))) {
      converged <- TRUE
      loglik_trace <- loglik_trace[1:iter]
      break
    }
  }
  
  if (!is.null(best_par)) {
    pi_k <- best_par$pi
    delta <- best_par$delta
    alpha <- best_par$alpha
    U <- best_par$U
    s <- best_par$s
    psi <- best_par$psi
    r <- best_par$r
  }
  
  cl <- max.col(r)
  
  list(
    cluster = cl,
    responsibilities = r,
    pi = pi_k,
    delta = delta,
    alpha = alpha,
    U = U,
    s = s,
    psi = psi,
    loglik = loglik_trace,
    best_loglik = best_ll,
    converged = converged,
    settings = list(K = K, q = q, G = G, psi_floor = psi_floor, ridge = ridge)
  )
}


############################################################
## 2. Gene module construction by MAD + hierarchical clustering
############################################################

fit_gene_cluster_mad_hc <- function(
    X, p_keep = 3000, r0 = 0.3, min_size = 30,
    method_cor = "pearson", plot_mad = TRUE
) {
  X <- as.matrix(X)
  
  mad_g <- apply(X, 2, mad, na.rm = TRUE)
  if (plot_mad) {
    plot(mad_g, pch = 16, cex = 0.4,
         main = "MAD per gene", ylab = "MAD", xlab = "Gene index")
  }
  
  keep <- order(mad_g, decreasing = TRUE)[1:min(p_keep, ncol(X))]
  X_sub <- X[, keep, drop = FALSE]
  
  X_gene <- scale(t(X_sub))
  C <- cor(t(X_gene), use = "pairwise.complete.obs", method = method_cor)
  D <- sqrt(pmax(1 - C, 0))
  diag(D) <- 0
  hc <- hclust(as.dist(D), method = "average")
  
  gene_cluster <- cutree(hc, h = 1 - r0)
  tab <- table(gene_cluster)
  gene_cluster[gene_cluster %in% names(tab[tab < min_size])] <- 0
  
  ## Remove unassigned genes
  keep_assigned <- gene_cluster != 0
  X_assigned <- X_sub[, keep_assigned, drop = FALSE]
  gc_assigned <- gene_cluster[keep_assigned]
  
  ## Re-encode module IDs as 1..G
  old_ids <- sort(unique(gc_assigned))
  new_ids <- seq_along(old_ids)
  map <- setNames(new_ids, old_ids)
  gc_recode <- as.integer(map[as.character(gc_assigned)])
  
  out <- list(
    X_assigned = X_assigned,
    gene_cluster = gc_recode,
    keep = keep[keep_assigned],
    hc = hc,
    gene_cluster_old = gc_assigned,
    map_old2new = map,
    settings = list(p_keep = p_keep, r0 = r0, min_size = min_size, method_cor = method_cor)
  )
  return(out)
}


############################################################
## 3. BIC
############################################################

compute_BIC_msd_mfm <- function(fit, X, gene_cluster, K, q) {
  n <- nrow(X)
  p <- ncol(X)
  G <- length(unique(gene_cluster))
  
  loglik <- fit$best_loglik
  d <- (K - 1) + K * G + G * (q - 1) + 3 * p
  BIC <- -2 * loglik + d * log(n)
  
  return(list(BIC = BIC))
}


############################################################
## 4. Construct modules on SLE + RA jointly
############################################################

fit0 <- fit_gene_cluster_mad_hc(
  X = X_both,
  p_keep = 3000,
  r0 = 0.3,
  min_size = 20,
  plot_mad = TRUE
)

X_bothS <- fit0$X_assigned
gene_cluster <- fit0$gene_cluster

table(gene_cluster)


############################################################
## 5. BIC grid search
############################################################

Ks <- 2:12
qs <- 2:10

results <- data.frame(
  K = integer(),
  q = integer(),
  loglik = numeric(),
  BIC = numeric()
)

for (K in Ks) {
  for (q in qs) {
    cat(sprintf("Fitting K = %d, q = %d\n", K, q))
    
    fit <- fit_msd_mfm_ecm(
      X = X_bothS,
      gene_cluster = gene_cluster,
      K = K,
      q = q,
      max_iter = 200,
      tol = 1e-5,
      verbose = FALSE
    )
    
    bic <- compute_BIC_msd_mfm(
      fit = fit,
      X = X_bothS,
      gene_cluster = gene_cluster,
      K = K,
      q = q
    )
    
    results <- rbind(
      results,
      data.frame(
        K = K,
        q = q,
        loglik = fit$best_loglik,
        BIC = bic$BIC
      )
    )
  }
}

results_sorted <- results[order(results$BIC), ]
results_sorted$delta_BIC <- results_sorted$BIC - min(results_sorted$BIC)
print(results_sorted)


############################################################
## 6. Plot delta BIC
############################################################

library(dplyr)
library(ggplot2)

df_plot <- results_sorted %>%
  arrange(q, K)

ggplot(df_plot, aes(x = K, y = delta_BIC, color = factor(q), group = q)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  labs(
    x = "Number of clusters (K)",
    y = expression(Delta * BIC),
    color = "q",
    title = expression(Delta * BIC ~ "for different combinations of K and q")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(face = "bold", size = 11)
  )


############################################################
## 7. Fit final model using best BIC
############################################################

best_model <- results_sorted[1, ]
print(best_model)

fit_final <- fit_msd_mfm_ecm(
  X = X_bothS,
  gene_cluster = gene_cluster,
  K = best_model$K,
  q = best_model$q,
  max_iter = 200,
  tol = 1e-5,
  verbose = FALSE
)

print(table(fit_final$cluster))
print(fit_final$alpha)

## Disease distribution across discovered clusters
cluster_result <- data.frame(
  disease = disease_label,
  cluster = fit_final$cluster
)

print(table(cluster_result$disease, cluster_result$cluster))
print(prop.table(table(cluster_result$disease, cluster_result$cluster), margin = 1))
print(prop.table(table(cluster_result$disease, cluster_result$cluster), margin = 2))


############################################################
## 8. Plot alpha heatmap and cluster sizes for SLE + RA
############################################################

plot_alpha_and_cluster <- function(fit_final,
                                   disease_label = NULL,
                                   digits = 2,
                                   main_heat = "(a) Module effects across SLE + RA subtypes",
                                   main_bar  = "(b) Sample size of SLE + RA subtype",
                                   xlab_heat = "Module",
                                   ylab_heat = "Subtype",
                                   xlab_bar  = "Subtype",
                                   ylab_bar  = "Number of samples",
                                   width_ratio = c(2.6, 1),
                                   cex_cell = 0.60,
                                   cex_axis = 0.85,
                                   cex_main = 1.10) {
  
  alpha <- as.matrix(fit_final$alpha)
  if (!is.matrix(alpha)) stop("fit_final$alpha must be a matrix.")
  
  row_hc <- hclust(dist(alpha))
  col_hc <- hclust(dist(t(alpha)))
  A <- alpha[row_hc$order, col_hc$order, drop = FALSE]
  
  if (is.null(rownames(A)) || any(rownames(A) == "")) {
    rownames(A) <- paste0("Subtype ", row_hc$order)
  }
  if (is.null(colnames(A)) || any(colnames(A) == "")) {
    colnames(A) <- paste0("M", col_hc$order)
  }
  
  zlim <- range(A, finite = TRUE)
  cols <- colorRampPalette(c("blue", "white", "red"))(101)
  
  cluster_tab <- table(fit_final$cluster)
  K <- nrow(alpha)
  cluster_tab <- cluster_tab[as.character(seq_len(K))]
  cluster_tab[is.na(cluster_tab)] <- 0
  
  layout(matrix(c(1, 2), 1, 2), widths = width_ratio)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  ## Heatmap
  par(mar = c(7, 6.5, 3.5, 4.0))
  image(
    x = seq_len(ncol(A)),
    y = seq_len(nrow(A)),
    z = t(A[nrow(A):1, , drop = FALSE]),
    col = cols,
    zlim = zlim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  title(main = main_heat, cex.main = cex_main, font.main = 2)
  axis(1, at = seq_len(ncol(A)), labels = colnames(A), las = 2, cex.axis = cex_axis)
  axis(2, at = seq_len(nrow(A)), labels = rev(rownames(A)), las = 2, cex.axis = cex_axis)
  mtext(xlab_heat, side = 1, line = 4.5, font = 2, cex = 1.05)
  mtext(ylab_heat, side = 2, line = 4.5, font = 2, cex = 1.05)
  box()
  
  vals <- round(A, digits)
  norm <- (A - zlim[1]) / diff(zlim)
  txt_col <- ifelse(norm < 0.25 | norm > 0.75, "white", "black")
  
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(ncol(A))) {
      y_pos <- nrow(A) - i + 1
      text(x = j, y = y_pos,
           labels = format(vals[i, j], nsmall = digits),
           cex = cex_cell, font = 2, col = txt_col[i, j])
    }
  }
  
  usr <- par("usr")
  x0 <- usr[2] + 0.08
  x1 <- usr[2] + 0.28
  y0 <- usr[3]
  y1 <- usr[4]
  nleg <- length(cols)
  ys <- seq(y0, y1, length.out = nleg + 1)
  
  par(xpd = NA)
  for (kk in seq_len(nleg)) rect(x0, ys[kk], x1, ys[kk + 1], col = cols[kk], border = NA)
  rect(x0, y0, x1, y1, border = "black")
  ticks <- pretty(zlim, n = 5)
  tick_pos <- y0 + (ticks - zlim[1]) / diff(zlim) * (y1 - y0)
  axis(4, at = tick_pos, labels = ticks, las = 2, cex.axis = 0.75)
  mtext("alpha", side = 4, line = 2.2, font = 2, cex = 1.0)
  par(xpd = FALSE)
  
  ## Barplot
  par(mar = c(7, 4.8, 4.5, 1.5))
  bp <- barplot(
    cluster_tab,
    col = "grey80",
    border = "black",
    ylim = c(0, max(cluster_tab) * 1.18),
    xlab = "",
    ylab = "",
    main = ""
  )
  mtext(main_bar, side = 3, line = 1.5, font = 2, cex = cex_main)
  axis(1, at = bp, labels = names(cluster_tab), las = 1, cex.axis = cex_axis)
  mtext(xlab_bar, side = 1, line = 4.5, font = 2, cex = 1.05)
  mtext(ylab_bar, side = 2, line = 3.5, font = 2, cex = 1.05)
  text(bp, cluster_tab, labels = as.integer(cluster_tab), pos = 3, cex = 0.95, font = 2)
  
  invisible(list(alpha_ordered = A, cluster_tab = cluster_tab, zlim = zlim))
}

plot_alpha_and_cluster(fit_final, disease_label = disease_label, digits = 2)


############################################################
## 9. Optional: stacked disease composition per cluster
############################################################

cluster_disease_tab <- table(cluster_result$cluster, cluster_result$disease)
barplot(
  t(cluster_disease_tab),
  beside = FALSE,
  legend.text = TRUE,
  main = "Disease composition of each discovered cluster",
  xlab = "Cluster",
  ylab = "Number of samples"
)


############################################################
## 10. Gene module enrichment analysis
############################################################

analyze_gene_modules <- function(
    X,
    gene_cluster,
    species = c("human", "mouse"),
    id_type = c("SYMBOL", "ENSEMBL"),
    do_GO = TRUE,
    do_KEGG = TRUE,
    ont = "BP",
    min_size = 10,
    q_cutoff = 0.05,
    top_terms = 10,
    compute_score = TRUE,
    plot_each_module = FALSE
) {
  species <- match.arg(species)
  id_type <- match.arg(id_type)
  
  X <- as.matrix(X)
  if (is.null(colnames(X))) stop("X must have colnames = gene IDs, e.g. SYMBOL or ENSEMBL.")
  if (length(gene_cluster) != ncol(X)) stop("length(gene_cluster) must equal ncol(X).")
  gene_cluster <- as.integer(gene_cluster)
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Install Bioconductor package: clusterProfiler")
  if (!requireNamespace("enrichplot", quietly = TRUE)) stop("Install Bioconductor package: enrichplot")
  
  if (species == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Install Bioconductor package: org.Hs.eg.db")
    OrgDb <- get("org.Hs.eg.db", asNamespace("org.Hs.eg.db"))
    kegg_org <- "hsa"
  } else {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Install Bioconductor package: org.Mm.eg.db")
    OrgDb <- get("org.Mm.eg.db", asNamespace("org.Mm.eg.db"))
    kegg_org <- "mmu"
  }
  
  gene_id <- colnames(X)
  modules <- split(gene_id, gene_cluster)
  module_sizes <- sapply(modules, length)
  
  sym2ent <- clusterProfiler::bitr(
    gene_id,
    fromType = id_type,
    toType = "ENTREZID",
    OrgDb = OrgDb
  )
  
  universe_ent <- unique(sym2ent$ENTREZID)
  
  get_entrez <- function(genes_m) {
    ent <- sym2ent$ENTREZID[match(genes_m, sym2ent[[id_type]])]
    unique(stats::na.omit(ent))
  }
  
  go_list <- list()
  kegg_list <- list()
  summary_tbl <- data.frame(
    module = names(modules),
    size = as.integer(module_sizes),
    n_mapped_entrez = NA_integer_,
    top_GO = NA_character_,
    top_KEGG = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (m in names(modules)) {
    genes_m <- modules[[m]]
    ent_m <- get_entrez(genes_m)
    summary_tbl$n_mapped_entrez[summary_tbl$module == m] <- length(ent_m)
    
    if (length(ent_m) < min_size) next
    
    if (do_GO) {
      ego <- clusterProfiler::enrichGO(
        gene = ent_m,
        universe = universe_ent,
        OrgDb = OrgDb,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        qvalueCutoff = q_cutoff,
        readable = TRUE
      )
      go_list[[m]] <- ego
      if (!is.null(ego) && nrow(ego@result) > 0) {
        top_desc <- ego@result$Description[1:min(top_terms, nrow(ego@result))]
        summary_tbl$top_GO[summary_tbl$module == m] <- paste(top_desc, collapse = "; ")
      }
      if (plot_each_module && !is.null(ego) && nrow(ego@result) > 0) {
        print(enrichplot::dotplot(ego, showCategory = min(top_terms, nrow(ego@result))) +
                ggplot2::ggtitle(paste0("Module ", m, " GO-", ont)))
      }
    }
    
    if (do_KEGG) {
      ek <- try(
        clusterProfiler::enrichKEGG(
          gene = ent_m,
          organism = kegg_org,
          universe = universe_ent,
          pAdjustMethod = "BH",
          qvalueCutoff = q_cutoff
        ),
        silent = TRUE
      )
      if (!inherits(ek, "try-error")) {
        kegg_list[[m]] <- ek
        if (!is.null(ek) && nrow(ek@result) > 0) {
          top_desc <- ek@result$Description[1:min(top_terms, nrow(ek@result))]
          summary_tbl$top_KEGG[summary_tbl$module == m] <- paste(top_desc, collapse = "; ")
        }
        if (plot_each_module && !is.null(ek) && nrow(ek@result) > 0) {
          print(enrichplot::dotplot(ek, showCategory = min(top_terms, nrow(ek@result))) +
                  ggplot2::ggtitle(paste0("Module ", m, " KEGG")))
        }
      }
    }
  }
  
  module_score <- NULL
  if (compute_score) {
    module_score <- sapply(modules, function(g) rowMeans(X[, g, drop = FALSE]))
    module_score <- as.matrix(module_score)
    colnames(module_score) <- paste0("M", names(modules))
  }
  
  out <- list(
    modules = modules,
    module_sizes = module_sizes,
    mapping = sym2ent,
    GO = go_list,
    KEGG = kegg_list,
    summary = summary_tbl,
    module_score = module_score
  )
  return(out)
}


############################################################
## 11. Run enrichment on SLE + RA modules
############################################################
library(clusterProfiler)
res <- analyze_gene_modules(
  X = X_bothS,
  gene_cluster = gene_cluster,
  species = "human",
  id_type = "SYMBOL",
  do_GO = TRUE,
  do_KEGG = TRUE,
  ont = "BP",
  min_size = 10,
  q_cutoff = 0.05,
  top_terms = 10,
  compute_score = TRUE,
  plot_each_module = FALSE
)

print(res$module_sizes)
print(res$summary)

## Example: top GO terms for module 3, if available
if ("3" %in% names(res$GO) && nrow(res$GO[["3"]]@result) > 0) {
  print(res$GO[["3"]]@result[1:min(10, nrow(res$GO[["3"]]@result)),
                              c("Description", "p.adjust", "Count")])
}

## Module score heatmap
print(dim(res$module_score))
heatmap(scale(res$module_score), Colv = NA)


############################################################
## 12. Save key outputs
############################################################

write.csv(results_sorted, "SLE_RA_BIC_results.csv", row.names = FALSE)
write.csv(cluster_result, "SLE_RA_cluster_assignments.csv", row.names = FALSE)
write.csv(res$summary, "SLE_RA_module_enrichment_summary.csv", row.names = FALSE)

saveRDS(
  list(
    fit0 = fit0,
    fit_final = fit_final,
    results_sorted = results_sorted,
    cluster_result = cluster_result,
    enrichment = res
  ),
  file = "SLE_RA_MSD_MFM_results.rds"
)
