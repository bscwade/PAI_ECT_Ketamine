require(ggplot2)
require(gridExtra)
require(reshape2)

basedir <- '~/Google Drive/My Drive/MGH/Studies/Mike_PAI/'

cohens_curve <- function(results, pai_threshold){
  
  pai_indices <- which(results$pai>=pai_threshold)
  
  # report effect size of difference between end-QIDS between those assigned to optimal versus non-optimal treatment
  t_test <- t.test(results$QIDS_end[pai_indices] ~ results$factual_better[pai_indices],
                   alternative='greater')
  md <- abs(t_test$estimate[1] - t_test$estimate[2])
  sdp <- sd_pooled(results$QIDS_end[results$factual_better==TRUE][pai_indices], results$QIDS_end[results$factual_better==FALSE][pai_indices])
  cohens_d <- md/sdp
  
  # get sample size for optimial and non-optimal group subsets
  tab_ono <- table(results$factual_better[pai_indices])
  n_nonoptimal <- tab_ono[names(tab_ono)==FALSE]
  n_optimal <- tab_ono[names(tab_ono)==TRUE]
  
  df <- data.frame(
    cd=cohens_d,
    md=md,
    n_nonoptimal=n_nonoptimal,
    n_optimal=n_optimal,
    threshold_value=pai_threshold,
    p_value=t_test$p.value
  )
  
  return(df)
  
}

dirs <- dir(paste0(basedir, 'results/mlsc_pai_results/varsel_07/'))

# load PAI results from MLSC
load_results <- function(basedir, resdir){
  results_list <- list()
  for(i in 1:470){
    print(i)
    fid <- grep(paste0('results_', i,'.Rdata'), dir(paste0(basedir, 'results/mlsc_pai_results/varsel_07/', resdir, '/'), full.names = TRUE), value = TRUE)
    print(fid)
    load(fid)
    results_list[[i]] <- results
    rm(results)
  }
  results_concat <- do.call(rbind, results_list)
  return(results_concat)
}


results_list <- lapply(dirs, function(d){
  
  res_tmp <- load_results(basedir = basedir, resdir = d)
  return(res_tmp)
  
})

cohens_data_list <- lapply(1:length(results_list), function(f){
  
  # cohens curve
  c_curve_data <- lapply(seq(0, 3, 0.1), function(t){
    tryCatch(cohens_curve(results = results_list[[f]], pai_threshold = t), error=function(e) NA)
  })
  
  cd_df <- do.call(rbind, c_curve_data)
  cd_df$p_adj <- p.adjust(cd_df$p_value, method='fdr')
  
  cd_df$match <- rep(dirs[[f]], nrow(cd_df))
  
  return(cd_df)
  
})

cdf <- do.call(rbind, cohens_data_list)


# generate break positions
breaks <- c(seq(0, 3, by=0.5))
# and labels
labels <- as.character(breaks)

# add column for adjusted significance
cdf$significant <- as.factor(ifelse(cdf$p_adj<=0.05, "Significant", "Non-Significant"))

p1 <- ggplot(cdf, aes(x=threshold_value, y=md, group=match, color=match)) + 
  geom_line() +
  geom_point(size=3, aes(shape=significant)) + 
  scale_shape_manual(values=c(1, 19)) + 
  scale_x_continuous(limits = c(0, 3), breaks = breaks, labels = labels, name = "PAI Threshold") + 
  ylab("Mean Difference: QIDS") + 
  ggtitle("Mean QIDS Difference by PAI Threshold") + 
  theme_bw() 

p2 <- ggplot(cdf, aes(x=threshold_value, y=cd, group=match, color=match)) + 
  geom_line() +
  geom_point(size=3, aes(shape=significant)) + 
  scale_shape_manual(values=c(1, 19)) + 
  scale_x_continuous(limits = c(0, 3), breaks = breaks, labels = labels, name = "PAI Threshold") + 
  ylab("Cohen's D") + 
  ggtitle("Effect Size by PAI Threshold") + 
  theme_bw() 

grid.arrange(p1, p2, ncol=2)

## Find minimal detectable effect
cdf_ss <- cdf[cdf$match=='match_all',]

cd_detectable <- sapply(1:nrow(cdf_ss), function(r){
  
  n <- min(c(cdf_ss$n_nonoptimal[r], cdf_ss$n_optimal[r]))
  
  power <- tryCatch(
    pwr.t.test(n = n, 
               sig.level = 0.05, 
               power = 0.8, 
               type='two.sample', 
               alternative = 'greater'),
    error=function(e) NA)
  
  es <- tryCatch(power$d, error=function(e) NA)
  
  return(es)
  
})

plot(cdf_ss$threshold_value, cd_detectable, type='o', col='red', pch=19)
text(cdf_ss$threshold_value, cd_detectable+.1, labels=signif(cd_detectable, 2))








