#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])

require(ranger)
require(fastDummies)
require(wakefield)
require(pacman)
require(MatchIt)
require(knitr)
require(tableone)
require(caret)
require(nestedcv)
require(MLmetrics)

basedir <- '/autofs/cluster/neuromod/personal_folders/benjamin_wade/projects/pai/'
#basedir <- '~/Google Drive/My Drive/MGH/Studies/Mike_PAI/'
n_trees <- 1000
treatment_a <- 'ect'
treatment_b <- 'ket'

format_data <- function(treatment_a, treatment_b, basedir){
  
  setwd(basedir)
  
  # get data
  treatment_a_data <- read.csv(grep(treatment_a, dir('data/', full.names = T), value = T, ignore.case = T))
  treatment_b_data <- read.csv(grep(treatment_b, dir('data/', full.names = T), value = T, ignore.case = T))
  
  # add treatment identifier
  treatment_a_data$treatment <- rep(treatment_a, nrow(treatment_a_data))
  treatment_b_data$treatment <- rep(treatment_b, nrow(treatment_b_data))
  
  # Calculate end-treatment QIDS
  treatment_a_data$QIDS_end <- sapply(1:nrow(treatment_a_data), function(s){
    Q0 <- treatment_a_data$QIDS_0[s]
    PC <- treatment_a_data$QIDS_min_pct_chg_acute[s]/100
    return( Q0 + (PC * Q0) )
  })
  
  treatment_b_data$QIDS_end <- sapply(1:nrow(treatment_b_data), function(s){
    Q0 <- treatment_b_data$QIDS_0[s]
    PC <- treatment_b_data$QIDS_min_pct_chg_acute[s]/100
    return( Q0 + (PC * Q0) )
  })
  
  # merge and format data including boolean group column for matching
  vars_keep <- intersect(names(treatment_a_data), names(treatment_b_data))
  vars_keep <- vars_keep[!vars_keep %in% c('QIDS_remit', 'QIDS_response', 'QIDS_min_pct_chg_acute', 'MDD', 'female_below_50')]
  data <- rbind(treatment_a_data[, vars_keep], treatment_b_data[, vars_keep])
  data <- subset(data, select = -n_before_meds) # drop incorrect n_before_meds
  
  # read in other/corrected demographic measures
  data_other <- read.csv('data/Three_arm_CMI_data.csv')
  data_other <- data_other[data_other$Study.ID %in% data$Study.ID, ]
  data_other <- data_other[, c('Study.ID', 'n_before_meds', grep('Hispanic|Race', names(data_other), value = T))]
  data <- merge(data, data_other, by = 'Study.ID')
  
  # drop nzv variables
  data <- data[, !names(data) %in% c('Study.ID', 'MDD_bipolar', 'MDD_PTSD', 'MDD_alc_dep', 'Not_Hispanic', names(data[, sapply(1:ncol(data), function(x) length(unique(data[, x])))==1]))]
  data$group <- as.logical(data$treatment == treatment_b)
  
  return(data)
  
}

match_data <- function(data){
  
  # match samples on baseline QIDS
  set.seed(1234)
  #match.it <- matchit(group ~ QIDS_0, data = data, method="nearest", ratio=1) 
  #match.it <- matchit(group ~ QIDS_0 + BASIS_psychosis_score + age, data = data, method="nearest", ratio=1) # ns
  #match.it <- matchit(group ~ QIDS_0 + Patient_Status + age, data = data, method="nearest", ratio=1) # significant 
  #match.it <- matchit(group ~ QIDS_0 + BASIS_psychosis_score + Patient_Status, data = data, method="nearest", ratio=1) #ns
  match.it <- matchit(group ~ QIDS_0 + BASIS_psychosis_score + Patient_Status + age, data = data, method="nearest", ratio=1) # significant
  a <- summary(match.it)
  df_match <- match.data(match.it)[1:ncol(data)]
  
  return(df_match)
  
}

make_data_clf <- function(df_match, treatment_b){
  data_clf <- dummy_cols(df_match, select_columns = c('treatment'))
  data_clf <- data_clf[, !names(data_clf) %in% c('treatment',  paste0('treatment_', treatment_b), 'group', 'Gender', grep('during|after', names(data_clf), value=T))]
  return(data_clf)
}

ranger_nzv_filter <- function (y, x, nfilter = NULL, type = c("index", "names", "full"), 
                               num.trees = 1000, mtry = ncol(x) * 0.2, ...) 
{
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' must be installed to use this filter", 
         call. = FALSE)
  }
  type <- match.arg(type)
  
  # nzv filter
  x <- x[, !names(x) %in% caret::nzv(x, names=T)]
  
  nfilter <- round(ncol(x) * 0.7)
  
  fit <- ranger::ranger(x = x, y = y, importance = "permutation", 
                        num.trees = num.trees, mtry = mtry, write.forest = FALSE, 
                        verbose = FALSE, num.threads = 1, ...)
  vi <- fit$variable.importance
  names(vi) <- if (type == "index") 
    1:ncol(x)
  else colnames(x)
  if (type == "full") 
    return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) 
    vi <- vi[1:min(nfilter, length(vi))]
  out <- names(vi)
  if (type == "index") 
    out <- as.integer(out)
  out
}

data <- format_data(treatment_a = treatment_a, treatment_b = treatment_b, basedir = basedir)
df_match <- match_data(data)
data_clf <- make_data_clf(df_match = df_match, treatment_b = treatment_b)

# get x, y
y <- data_clf$QIDS_end
x <- subset(data_clf, select = -QIDS_end)

# get permuted outcome
set.seed(index)
y_permuted <- sample(y, size = length(y), replace = FALSE)

# nested CV using caret
# tg <- expand.grid(mtry = c(2, 3, 5), 
#                   splitrule = c('variance', 'extratrees'),
#                   min.node.size=c(2, 5, 10))

ctrl <- trainControl(method="repeatedcv",
                     number=10,
                     repeats=1,
                     returnResamp="all",
                     classProbs=FALSE,
                     savePredictions=TRUE)
set.seed(index)
ncv <- nestcv.train(y = y_permuted, x = x,
                    method = "ranger",
                    savePredictions = "final",
                    filterFUN = ranger_nzv_filter,
                    filter_options = list(type='names'),
                    #tuneGrid = tg, 
                    num.trees=n_trees,
                    importance='permutation',
                    trControl = ctrl,
                    cv.cores = 1)

R2_permutation <- MLmetrics::R2_Score(y_pred = ncv$final_fit$finalModel$predictions, y_true = y)

save(R2_permutation, file = paste0(basedir, sprintf('results/permutations/qids_psychosis_inpatient_age/permutation_%s.Rdata', index)))



