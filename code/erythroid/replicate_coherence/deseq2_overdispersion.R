rm(list=ls())

## ----------------------------------------------------------
## Fit dispersion–mean trend: alpha_tr(mu_bar) = a1/mu_bar + a0
## ----------------------------------------------------------

fit_dispersion_trend <- function(mean_norm,
                                 alpha_hat,
                                 max_iter   = 10,
                                 min_disp   = 1e-8,
                                 link       = c("log", "identity")) {
  link <- match.arg(link)
  
  stopifnot(length(mean_norm) == length(alpha_hat))
  ok <- is.finite(mean_norm) & is.finite(alpha_hat) &
    mean_norm > 0 & alpha_hat > 0
  if (sum(ok) < 20L) {
    stop("Not enough finite positive genes to fit trend.")
  }
  
  d <- data.frame(
    disp     = alpha_hat[ok],
    inv_mean = 1 / mean_norm[ok]
  )
  
  use <- rep(TRUE, nrow(d))
  fit <- NULL
  
  for (iter in seq_len(max_iter)) {
    # fit only on currently used points
    fit <- glm(
      disp ~ inv_mean,
      data   = d[use, , drop = FALSE],
      family = Gamma(link = link)
    )
    
    fitted_sub <- pmax(fitted(fit), min_disp)
    
    # ratio must be computed on the same subset
    ratio_sub <- d$disp[use] / fitted_sub
    
    use_new_sub <- ratio_sub > 1e-4 & ratio_sub < 15
    
    # map back into global 'use'
    use_idx <- which(use)
    use[use_idx] <- use_new_sub
    
    # stop if trimming stabilized or we have too few points left
    if (all(use_new_sub) || sum(use) < 20L) {
      break
    }
  }
  
  # build trend on *all* genes (including not-ok as NA)
  if (link == "identity") {
    a0 <- unname(coef(fit)[1])
    a1 <- unname(coef(fit)[2])
    trend_ok <- a0 + a1 / mean_norm[ok]
  } else {  # link == "log"
    b0 <- unname(coef(fit)[1])
    b1 <- unname(coef(fit)[2])
    # always positive: alpha_tr(mu) = exp(b0 + b1 / mu)
    trend_ok <- exp(b0 + b1 / mean_norm[ok])
  }
  
  trend <- rep(NA_real_, length(mean_norm))
  trend[ok] <- trend_ok
  
  # enforce positivity
  if (any(trend_ok <= 0, na.rm = TRUE)) {
    min_pos <- min(trend_ok[trend_ok > 0], na.rm = TRUE)
    trend[trend <= 0] <- min_pos
  }
  
  list(
    trend = trend,
    ok    = ok,
    fit   = fit
  )
}


## ----------------------------------------------------------
## Shrink dispersions toward the trend (DESeq2-like)
## ----------------------------------------------------------

shrink_dispersions <- function(mean_norm,
                               alpha_hat,
                               df,
                               min_prior_var = 0.25,
                               outlier_sd    = 2,
                               link          = c("log", "identity")) {
  link <- match.arg(link)
  
  stopifnot(length(mean_norm) == length(alpha_hat))
  if (df <= 3)
    warning("df <= 3: prior-variance estimate may be less reliable.")
  
  # 1) Fit trend
  tr <- fit_dispersion_trend(mean_norm, alpha_hat, link = link)
  alpha_tr <- tr$trend
  
  # 2) Log-residuals around trend
  log_resid <- log(alpha_hat) - log(alpha_tr)
  slr <- mad(log_resid[is.finite(log_resid)],
             center   = 0,
             constant = 1)
  s2_lr <- slr^2
  
  # 3) Approx sampling variance of log-dispersion
  var_sample <- trigamma(df / 2)
  
  # 4) Prior variance on log(alpha)
  prior_var <- max(s2_lr - var_sample, min_prior_var)
  prior_sd  <- sqrt(prior_var)
  
  # 5) Flag dispersion outliers: far above trend
  is_outlier <- log(alpha_hat) > log(alpha_tr) + outlier_sd * slr
  if(any(is.na(is_outlier))){
    is_outlier[is.na(is_outlier)] <- FALSE
  }
  
  # 6) Normal–Normal shrinkage in log-space
  log_like  <- log(alpha_hat)
  log_prior <- log(alpha_tr)
  
  post_prec <- 1 / var_sample + 1 / prior_var
  post_mean <- (log_like / var_sample + log_prior / prior_var) / post_prec
  
  alpha_shrunk <- exp(post_mean)
  alpha_shrunk[is_outlier] <- alpha_hat[is_outlier]
  
  list(
    disp_raw     = alpha_hat,
    mean_norm    = mean_norm,
    disp_trend   = alpha_tr,
    disp_shrunk  = alpha_shrunk,
    prior_var    = prior_var,
    prior_sd     = prior_sd,
    var_sample   = var_sample,
    slr          = slr,
    is_outlier   = is_outlier
  )
}

#################

csv_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/erythroid_replicate_coherence/"
spliced_disp <- read.csv(paste0(csv_folder, "ery_overdisp_S.csv"), 
                         row.names = 1)
unspliced_disp <- read.csv(paste0(csv_folder, "ery_overdisp_U.csv"), 
                           row.names = 1)
spliced_depth <- read.csv(paste0(csv_folder, "spliced_depth.csv"), 
                          row.names = 1)
unspliced_depth <- read.csv(paste0(csv_folder, "unspliced_depth.csv"), 
                            row.names = 1)

# NOTE: 
# We are inversing since: 
# 1) The overdispersion's we're getting are inversed: https://github.com/linnykos/veloUncertainty/blob/main/code/erythroid/replicate_coherence/v4_ery_overdisp-2.R
# 2) Those values are 1/(result of glmGamPoi), which are in the parameterization of "larger = more variance" already.
## See https://rdrr.io/github/const-ae/glmGamPoi/man/glm_gp.html
df <- data.frame(
  gene = rownames(spliced_depth),
  spliced_disp = 1/spliced_disp[,1],
  spliced_depth = spliced_depth[,1],
  unspliced_disp = 1/unspliced_disp[,1],
  unspliced_depth = unspliced_depth[,1]
)

keep_disp <- intersect(
  which(!is.infinite(df$spliced_disp)),
  which(!is.infinite(df$unspliced_disp))
)
df <- df[keep_disp,]

keep_depth <- intersect(
  which(df$spliced_depth > 0.05),
  which(df$unspliced_depth > 0.05)
)
df <- df[keep_depth,]

dim(df) # 2030 genes

# plot(x = df$spliced_depth, y = df$spliced_disp)
# plot(x = df$spliced_depth, y = log(df$spliced_disp))
# plot(x = df$unspliced_depth, y = log(df$unspliced_disp))

#######################################

df_disp_spliced <- shrink_dispersions(
  mean_norm = df$spliced_depth,
  alpha_hat = df$spliced_disp,
  min_prior_var = 0.1,
  df        = 5, # purposely low to induce shrinkage
  link      = "log"
)

par(mfrow = c(1,2))
plot(x = df$spliced_disp, 
     y = df_disp_spliced$disp_shrunk, 
     asp = TRUE)

plot(x = log10(df$spliced_depth), 
     y = df$spliced_disp,
     ylab = "Dispersion",
     pch = 16,
     col = rgb(0.9, 0.9, 0.9))

ordering <- order(log10(df$spliced_depth))
lines(x = log10(df$spliced_depth[ordering]),
       y = df_disp_spliced$disp_trend[ordering],
       col = "blue")

idx <- which(abs(df$spliced_disp - df_disp_spliced$disp_shrunk) > 0.5)
points(x = log10(df$spliced_depth[idx]),
       y = df$spliced_disp[idx],
       col = "black",
       pch = 16)
points(x = log10(df$spliced_depth[idx]),
       y = df_disp_spliced$disp_shrunk[idx],
       col = "firebrick",
       pch = 16)

######

df_disp_unspliced <- shrink_dispersions(
  mean_norm = df$unspliced_depth,
  alpha_hat = df$unspliced_disp,
  min_prior_var = 0.1,
  df        = 5, # purposely low to induce shrinkage
  link      = "log"
)

par(mfrow = c(1,2))
plot(x = df$unspliced_disp, 
     y = df_disp_unspliced$disp_shrunk, 
     asp = TRUE)

plot(x = log10(df$unspliced_depth), 
     y = df$unspliced_disp,
     ylab = "Dispersion",
     pch = 16,
     col = rgb(0.9, 0.9, 0.9))

ordering <- order(log10(df$unspliced_depth))
lines(x = log10(df$unspliced_depth[ordering]),
      y = df_disp_unspliced$disp_trend[ordering],
      col = "blue")

idx <- which(abs(df$unspliced_disp - df_disp_unspliced$disp_shrunk) > 0.5)
points(x = log10(df$unspliced_depth[idx]),
       y = df$unspliced_disp[idx],
       col = "black",
       pch = 16)
points(x = log10(df$unspliced_depth[idx]),
       y = df_disp_unspliced$disp_shrunk[idx],
       col = "firebrick",
       pch = 16)

####################################

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/code/erythroid/replicate_coherence/csv/"

df$spliced_disp_shrunk <- df_disp_spliced$disp_shrunk
df$spliced_disp_trend <- df_disp_spliced$disp_trend
df$unspliced_disp_shrunk <- df_disp_unspliced$disp_shrunk
df$unspliced_disp_trend <- df_disp_unspliced$disp_trend

df <- df[,c("gene", 
            "spliced_depth", "spliced_disp", "spliced_disp_shrunk", "spliced_disp_trend",
            "unspliced_depth", "unspliced_disp", "unspliced_disp_shrunk", "unspliced_disp_trend")]

# flip it back prior to saving
disp_cols <- grep("disp", colnames(df))
for(j in disp_cols){
  df[,j] <- 1/df[,j]
}

write.csv(df, 
          file = paste0(out_folder, "erthyroid_shrunk_disp.csv"))
