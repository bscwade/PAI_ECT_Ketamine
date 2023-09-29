# Draws from results computed stored on MLSC using nested cross validation
# use varsel 07 @ 1000 trees

require(ggplot2)
require(reshape2)
require(ranger)
require(treeshap)
require(fastDummies)
require(effectsize)
require(gridExtra)
require(wakefield)
require(pacman)
require(MatchIt)
require(knitr)
require(tableone)
require(caret)
require(nestedcv)
require(MLmetrics)

basedir <- '~/Google Drive/My Drive/MGH/Studies/Mike_PAI/'
n_trees <- 1000
treatment_a <- 'ect'
treatment_b <- 'ket'

# load results from MLSC
load(paste0(basedir, 'results/mlsc_pai_results/main_results_qids_psychosis_inpatient_age.Rdata'))
rm(basedir)
basedir <- '~/Google Drive/My Drive/MGH/Studies/Mike_PAI/'


## Summary of cohort
# summary of data before matching
summary_vars <- names(data[, !names(data) %in% c('group', grep('during|after', names(data), value=T))])
cat_vars <- names(which(sapply(summary_vars, function(x) length(unique(data[[x]]))) <= 2))
tab <- CreateTableOne(data=data,
                      vars=summary_vars,
                      factorVars = cat_vars,
                      strata = 'treatment')

print("Summary of Unmatched Data")
print(kableone(tab))

# save as csv
tab_save <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab_save, file = '~/Desktop/demographics_table_not_matched.csv')
rm(tab)

# summary of data after matching
tab <- CreateTableOne(data=df_match,
                      vars=summary_vars,
                      factorVars = cat_vars,
                      strata = 'treatment')

print("Summary of Matched Data")
print(kableone(tab))

# save as csv
tab_save <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab_save, file = '~/Desktop/demographics_table_matched.csv')

# load PAI results from MLSC
load_results <- function(basedir){
  results_list <- list()
  for(i in 1:470){
    print(i)
    fid <- grep(paste0('results_', i,'.Rdata'), dir(paste0(basedir, 'results/mlsc_pai_results/varsel_07/match_qids_inpatient_age/'), full.names = TRUE), value = TRUE)
    print(fid)
    load(fid)
    results_list[[i]] <- results
    rm(results)
  }
  results_concat <- do.call(rbind, results_list)
  return(results_concat)
}

report <- function(results, pai_threshold){
  
  print(sprintf("Results Using PAI Threshold: %s", pai_threshold))
  
  pai_indices <- which(results$pai>=pai_threshold)
  
  # Plot end qids scores by those who received their predicted best treatment versus those who received their sub-optimal treatment
  #boxplot(results$QIDS_end[pai_indices] ~ results$factual_better[pai_indices], names=c('Non-optimal', 'Optimal'), ylab = 'QIDS End Treatment', xlab='') 
  plot_data <- results[pai_indices, ]
  p <- ggplot(plot_data, aes(factual_better, QIDS_end, fill=factual_better)) + 
    geom_boxplot() + 
    theme_bw()
  plot(p)
  
  
  print("T-test of End QIDS by Optimal and Non-optimal Treatments")
  print(t.test(results$QIDS_end[pai_indices] ~ results$factual_better[pai_indices], alternative = 'greater'))
  
  # report effect size of difference between end-QIDS between those assigned to optimal versus non-optimal treatment
  t_test <- t.test(results$QIDS_end[pai_indices] ~ results$factual_better[pai_indices], alternative = 'greater')
  md <- abs(t_test$estimate[1] - t_test$estimate[2])
  sdp <- sd_pooled(results$QIDS_end[results$factual_better==TRUE][pai_indices], results$QIDS_end[results$factual_better==FALSE][pai_indices])
  cohens_d <- md/sdp
  print(sprintf('Effect size (Cohen D): %s', cohens_d))
  
  print("Check whether predicted optimal treatment is simply confounded with treatment received")
  kable(table(results$factual_treatment, results$factual_better))
  print(chisq.test(table(results$factual_treatment, results$factual_better)))
  
}

cohens_curve <- function(results, pai_threshold){
  
  pai_indices <- which(results$pai>=pai_threshold)
  
  # report effect size of difference between end-QIDS between those assigned to optimal versus non-optimal treatment
  t_test <- t.test(results$QIDS_end[pai_indices] ~ results$factual_better[pai_indices], alternative = 'greater')
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

shaps_importance <- function(fit_base, data_clf){
  
  model_unified <- ranger.unify(fit_base$final_fit$finalModel, data_clf)
  treeshap_res <- treeshap(model_unified, data_clf, verbose = FALSE)
  plot_feature_importance(treeshap_res, max_vars = 10)
  
}

shaps_waterfall_plots <- function(fit_base, data_clf){
  
  model_unified <- ranger.unify(fit_base$final_fit$finalModel, data_clf)
  treeshap_res <- treeshap(model_unified, data_clf, verbose = FALSE)
  
  # shap plots
  set.seed(0)
  indices <- sample(x=470, size=10, replace=FALSE)
  plots <- lapply(indices, function(p){
    plt <- plot_contribution(treeshap_res, obs = p, title = sprintf('SHAP Breakdown for patient %s', p))
    plot(plt)
  })
  
}

shaps_pdp <- function(fit_base, data_clf, x){
  
  model_unified <- ranger.unify(fit_base$final_fit$finalModel, data_clf)
  treeshap_res <- treeshap(model_unified, data_clf, verbose = FALSE)
  
  for(i in x){
    plot(plot_feature_dependence(treeshap_res, i))
  }
  
}


shaps_interaction_plots <- function(fit_base, data_clf, x, treatment_a){
  
  model_unified <- ranger.unify(fit_base$final_fit$finalModel, data_clf)
  inter <- treeshap(model_unified, data_clf, interactions = T, verbose = FALSE)
  
  for(i in 1:length(x)){
    
    if(length(unique(data_clf[, x[i]])) > 2){
      plot(plot_interaction(inter, x[i], paste0('treatment_', treatment_a)))
    }else{
      plt <- plot_interaction(inter, x[i], paste0('treatment_', treatment_a))
      colnames(plt$data) <- c(x[i], paste0('treatment_', treatment_a),'interaction')
      plt$data[[x[i]]] <- as.factor(plt$data[[x[i]]])
      plt$data[[paste0('treatment_', treatment_a)]] <- as.factor(plt$data[[paste0('treatment_', treatment_a)]])
      p <- ggplot(plt$data, aes(x=get(x[i]), y=interaction, fill=get(paste0('treatment_', treatment_a)))) + 
        geom_boxplot() + 
        guides(fill=guide_legend(title=paste0('treatment_', treatment_a))) +
        xlab(sprintf('%s', x[i])) +
        theme_bw()
      plot(p)
      
    }
  }
  
}


## permutation results for "global" model
r2_observed <- R2_Score(y_pred = ncv$final_fit$finalModel$predictions, y_true = data_clf$QIDS_end)

fid_permutations <- dir(paste0(basedir, 'results/mlsc_pai_results/varsel_07/match_all/permutations/'), full.names = TRUE)
r2_permutation_distribution <- sapply(fid_permutations, function(f){
  load(f)
  return(R2_permutation)
})
hist(r2_permutation_distribution)
p <- (sum(r2_permutation_distribution >= r2_observed) + 1) / (length(r2_permutation_distribution) + 1)
##

results <- load_results(basedir = basedir)

## performance of LOOCV model for PAI
R2_Score(y_pred=results$pred_factual, y_true = results$QIDS_end)


hist(results$pai)

report(results = results, pai_threshold = 0)

# cohens curve
c_curve_data <- lapply(seq(0, max(results$pai), 0.1), function(t){
  tryCatch(cohens_curve(results = results, pai_threshold = t), error=function(e) NA)
})

cd_df <- do.call(rbind, c_curve_data)
cd_df$p_adj <- p.adjust(cd_df$p_value, method='fdr')

# generate break positions
breaks <- c(seq(0, 3, by=0.5))
# and labels
labels <- as.character(breaks)

ggplot(cd_df, aes(x=threshold_value, y=cd)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 3), breaks = breaks, labels = labels, name = "PAI Threshold") + 
  ylab("Cohen's D") + 
  ggtitle("Effect Size by PAI Threshold") + 
  theme_bw() 
  

#plot(cd_df$threshold_value, cd_df$cd, type="b" , bty="l", col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=1)


# Proportion of patients for whom their treatment was predicted optimal by treatment received
prop.table(table(results$factual_treatment, results$factual_better), margin = 1)

# proportion of patients who received their optimal treatment
table(results$factual_better)

mean(results[results$factual_better==TRUE, 'QIDS_end'])
mean(results[results$factual_better==FALSE, 'QIDS_end'])

print(kable(cd_df))

# ## SHAPS
# shaps_importance(fit_base = ncv, data_clf = data_clf)
shaps_waterfall_plots(fit_base = ncv, data_clf = data_clf)
# shaps_pdp(fit_base = ncv, data_clf = data_clf, x =  c('QIDS_0', 'BASIS_relationships_score',
#                                                              'BASIS_self_harm_score', 'BASIS_emotional_lability_score',
#                                                              'BASIS_psychosis_score', 'BASIS_substance_abuse_score',
#                                                              'age', 'n_diags', 'PD', 'treatment_ect'))

# shaps_interaction_plots(fit_base = ncv, data_clf = data_clf,
#                         x = c('QIDS_0', 'BASIS_relationships_score',
#                               'BASIS_self_harm_score', 'BASIS_emotional_lability_score',
#                               'BASIS_psychosis_score', 'BASIS_substance_abuse_score',
#                               'age', 'bipolar','insomnia', 'Female', 'PD', 'n_diags', 'MOCA_Score'),
#                         treatment_a = treatment_a)



