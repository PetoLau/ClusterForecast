# Forecasting methods ----

simCtree <- function(Y, freq = 48, h = 2) {
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  matrix_train <- data.table(Load = Y,
                             fourier(data_ts, K = c(2, 2)))
  
  fuur_test <- as.data.frame(fourier(data_ts, K = c(2, 2), h = h))
  
  matrix_test <- data.table(fuur_test)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_train,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.925,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, matrix_test)
  
  return(as.vector(pred_tree))
}

simSnaive <- function(Y, freq = 48, h = 2) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  pred <- snaive(data_ts, h = h)$mean

  return(as.vector(pred))
}

simEs <- function(Y, freq = 48, h = 2) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  pred <- es(data_ts, model = c("AAA", "ANA", "AAdA"), ic = "AICc", h = h,
             cfType = "MAE", intervals = "none", silent = "all")$forecast
  
  return(as.vector(pred))
}

MLRpred <- function(Y, model_matrix, h = 2, freq = 48) {
  
  model_matrix$Load <- Y
  lm_m <- lm(Load ~ 0 + Daily:Weekly, data = model_matrix)
  
  return(as.vector(predict(lm_m, model_matrix[head(tail(1:length(Y), 7*freq), h), .(Daily, Weekly)])))
}

simLinReg <- function(Y, FUN = MLRpred, freq = 48, h = 2) {
  
  N <- length(Y)
  daily <- rep(1:freq, N/freq)
  week <- N / (freq*7)
  weekly <- rep(rep(c(1:7), week), each = freq)
  
  matrix_train <- data.table(Daily = as.factor(daily),
                             Weekly = as.factor(weekly))
  
  pred <- FUN(Y, matrix_train, freq = freq, h = h)
  
  return(pred)
}

# All together ----
predAlldisagg <- function(Y, freq = 48, h = 2) {
  
  pred_naive <- simSnaive(Y, freq = freq, h = h)
  pred_es <- simEs(Y, freq = freq, h = h)
  pred_mlr <- simLinReg(Y, FUN = MLRpred, freq = freq, h = h)
  pred_ctree <- simCtree(Y, freq = freq, h = h)

  return(list(NAIVE = pred_naive,
              MLR = pred_mlr,
              ES = pred_es,
              CTREE = pred_ctree
         ))
}

ForecastClusterDisagg <- function(dataset, parallel = FALSE, freq = 48, h = 2, ncores = 8) {
  
  if(parallel == FALSE) {
    pred_clusters <- lapply(1:nrow(dataset), function(i) predAlldisagg(dataset[i,], freq = freq, h = h))
  } else {
    cl <- makeForkCluster(ncores, outfile = "")
    registerDoParallel(cl)
    pred_clusters <- parLapply(cl, 1:nrow(dataset), function(i) predAlldisagg(dataset[i,], freq = freq, h = h))
    if(!is.null(cl)){
      parallel::stopCluster(cl)
      cl <- c()
    }
  }
  
  return(pred_clusters)
}
