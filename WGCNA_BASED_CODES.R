############################################################
## 12. Comparison model: WGCNA modules on previous 979 genes
############################################################

fit_gene_cluster_wgcna <- function(
    X,
    power = NULL,
    minModuleSize = 20,
    mergeCutHeight = 0.25,
    networkType = "signed",
    corType = "pearson",
    plot_power = TRUE,
    seed = 1
) {
  set.seed(seed)
  
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Please install WGCNA: install.packages('WGCNA')")
  }
  
  WGCNA::allowWGCNAThreads()
  
  X <- as.matrix(X)
  datExpr <- X
  
  good <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
  if (!good$allOK) {
    datExpr <- datExpr[good$goodSamples, good$goodGenes, drop = FALSE]
  }
  
  if (is.null(power)) {
    powers <- c(1:10, seq(12, 30, by = 2))
    
    sft <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = networkType,
      corFnc = "cor",
      corOptions = list(
        use = "pairwise.complete.obs",
        method = "pearson"
      ),
      verbose = 5
    )
    
    if (plot_power) {
      plot(
        sft$fitIndices[, 1],
        -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        xlab = "Soft Threshold Power",
        ylab = "Scale Free Topology Model Fit, signed R^2",
        type = "n",
        main = "WGCNA soft-threshold selection"
      )
      text(
        sft$fitIndices[, 1],
        -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        labels = powers,
        col = "red"
      )
      abline(h = 0.8, col = "blue", lty = 2)
    }
    
    power <- sft$powerEstimate
    
    if (is.na(power)) {
      warning("WGCNA did not find a powerEstimate. Using power = 6.")
      power <- 6
    }
  }
  
  cat("Using WGCNA soft-threshold power =", power, "\n")
  
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    networkType = networkType,
    TOMType = networkType,
    corType = corType,
    minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    numericLabels = TRUE,
    pamRespectsDendro = TRUE,
    saveTOMs = FALSE,
    minKMEtoStay = 0,
    reassignThreshold = 0,
    verbose = 3
  )
  
  gene_cluster <- net$colors
  
  keep_assigned <- gene_cluster != 0
  
  if (sum(keep_assigned) == 0) {
    stop("WGCNA assigned all genes to grey module. Try smaller minModuleSize or fixed power.")
  }
  
  X_assigned <- datExpr[, keep_assigned, drop = FALSE]
  gc_assigned <- gene_cluster[keep_assigned]
  
  gc_recode <- as.integer(factor(gc_assigned))
  
  list(
    X_assigned = X_assigned,
    gene_cluster = gc_recode,
    keep = which(keep_assigned),
    wgcna_net = net,
    power = power,
    original_wgcna_colors = gc_assigned,
    settings = list(
      input_genes = ncol(X),
      retained_genes = ncol(X_assigned),
      power = power,
      minModuleSize = minModuleSize,
      mergeCutHeight = mergeCutHeight,
      networkType = networkType,
      corType = corType,
      minKMEtoStay = 0,
      reassignThreshold = 0
    )
  )
}


############################################################
## 13. Run WGCNA on previous 979 genes
############################################################

fit0_wgcna <- fit_gene_cluster_wgcna(
  X = X_bothS,
  power = NULL,
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  networkType = "signed",
  corType = "pearson",
  plot_power = TRUE,
  seed = 1
)

X_bothS_wgcna <- fit0_wgcna$X_assigned
gene_cluster_wgcna <- fit0_wgcna$gene_cluster

cat("Input genes from MAD-HC preprocessing:", ncol(X_bothS), "\n")
cat("Genes retained after WGCNA grey removal:", ncol(X_bothS_wgcna), "\n")
print(table(gene_cluster_wgcna))


############################################################
## 14. BIC grid search for WGCNA comparison model
############################################################

results_wgcna <- data.frame(
  K = integer(),
  q = integer(),
  loglik = numeric(),
  BIC = numeric()
)
Ks1=c(1:5)
for (K in Ks1) {
  for (q in qs) {
    cat(sprintf("[WGCNA on 979 genes] Fitting K = %d, q = %d\n", K, q))
    
    fit_w <- fit_msd_mfm_ecm(
      X = X_bothS_wgcna,
      gene_cluster = gene_cluster_wgcna,
      K = K,
      q = q,
      max_iter = 200,
      tol = 1e-5,
      verbose = FALSE
    )
    
    bic_w <- compute_BIC_msd_mfm(
      fit = fit_w,
      X = X_bothS_wgcna,
      gene_cluster = gene_cluster_wgcna,
      K = K,
      q = q
    )
    
    results_wgcna <- rbind(
      results_wgcna,
      data.frame(
        K = K,
        q = q,
        loglik = fit_w$best_loglik,
        BIC = bic_w$BIC
      )
    )
  }
}

results_wgcna_sorted <- results_wgcna[order(results_wgcna$BIC), ]
results_wgcna_sorted$delta_BIC <- results_wgcna_sorted$BIC - min(results_wgcna_sorted$BIC)

print(results_wgcna_sorted)


############################################################
## 15. Fit final WGCNA comparison model
############################################################

best_model_wgcna <- results_wgcna_sorted[1, ]
print(best_model_wgcna)

fit_final_wgcna <- fit_msd_mfm_ecm(
  X = X_bothS_wgcna,
  gene_cluster = gene_cluster_wgcna,
  K = best_model_wgcna$K,
  q = best_model_wgcna$q,
  max_iter = 200,
  tol = 1e-5,
  verbose = FALSE
)

cluster_result_wgcna <- data.frame(
  disease = disease_label,
  cluster = fit_final_wgcna$cluster
)

print(table(fit_final_wgcna$cluster))
print(table(cluster_result_wgcna$disease, cluster_result_wgcna$cluster))
print(prop.table(table(cluster_result_wgcna$disease, cluster_result_wgcna$cluster), margin = 1))
print(prop.table(table(cluster_result_wgcna$disease, cluster_result_wgcna$cluster), margin = 2))


############################################################
## 16. Compare MAD-HC modules vs WGCNA modules
############################################################

comparison_summary <- data.frame(
  model = c("MAD_HC_modules", "WGCNA_on_MADHC_genes"),
  best_K = c(best_model$K, best_model_wgcna$K),
  best_q = c(best_model$q, best_model_wgcna$q),
  best_loglik = c(fit_final$best_loglik, fit_final_wgcna$best_loglik),
  best_BIC = c(best_model$BIC, best_model_wgcna$BIC),
  n_modules = c(length(unique(gene_cluster)), length(unique(gene_cluster_wgcna))),
  n_genes = c(ncol(X_bothS), ncol(X_bothS_wgcna))
)

print(comparison_summary)

plot_alpha_and_cluster(
  fit_final_wgcna,
  disease_label = disease_label,
  digits = 2,
  main_heat = "(a) WGCNA module effects across clusters",
  main_bar = "(b) WGCNA cluster sample size",
  cex_bar_text = 0.60
)


############################################################
## 17. Save comparison results
############################################################

saveRDS(
  list(
    fit0_mad_hc = fit0,
    fit_final_mad_hc = fit_final,
    results_mad_hc = results_sorted,
    cluster_result_mad_hc = cluster_result,
    
    fit0_wgcna = fit0_wgcna,
    fit_final_wgcna = fit_final_wgcna,
    results_wgcna = results_wgcna_sorted,
    cluster_result_wgcna = cluster_result_wgcna,
    
    comparison_summary = comparison_summary
  ),
  file = "SLE_RA_MSD_MFM_comparison_MADHC_vs_WGCNA_979genes.rds"
)