# main file ----
library(data.table)
library(cluster)
library(clusterCrit)
library(TSrepr) # devtools::install_github("PetoLau/TSrepr")
library(smooth)
library(forecast)
library(party)
# Parallel computations can be used for Linux based machines - load these packages then please:
# library(parallel)
# library(doParallel)

# functions for forecasting and clustering
source("optimClusters.R")
source("forecasting.R")

seas <- 48 # frequency
win <- 21 # size of window in days
h <- 2 # length of forecast
week <- 7
dataset_load <- firstSetTrain # your dataset - in every row of a data.table is a time series (consumer)
win_hour <- win*24
n_days <- ncol(dataset_load)/seas - win # how many days of forecasts can we produce
n_hours <- ncol(dataset_load)/h - win*seas/h # how many forecasts can we produce

# main cycle ----
for(i in 0:(n_hours-1)) {
  
  # subsetting dataset
  datas <- data.matrix(dataset_load[, ((i*h)+1):((i+win_hour)*h), with = F])
  datas_test <- data.matrix(dataset_load[, (((i+win_hour)*h)+1):(((i+win_hour)*h)+h), with = F]) # possible to use for direct accuracy computing
  
  # Normalisation of time series
  data_norm_list <- lapply(1:nrow(datas), function(i) norm_z_list(datas[i,]))
  res_mean <- unlist(sapply(data_norm_list, `[`, c('mean')))
  res_sd <- unlist(sapply(data_norm_list, `[`, c('sd')))
  res_data <- matrix(unlist(sapply(data_norm_list, `[`, c('norm_values'))), nrow = nrow(datas), byrow = T)
  
  # Clustering is used only every new day
  if (i %% 24 == 0) {
    
    # Compute time series representations
    ts_repr <- repr_matrix(res_data, func = repr_lm, args = list(method = "l1", freq = c(seas, seas*7)))

    # Possible to use other time series representations:
    # ts_repr <- repr_matrix(res_data, func = repr_gam, args = list(freq = c(seas, seas*7)))
    # ts_repr <- repr_matrix(datas, func = repr_feaclip, windowing = T, win_size = seas)
    # ts_repr <- repr_matrix(res_data, func = repr_seas_profile, args = list(func = median, freq = seas))
  
    # Clustering - K-means or K-medoids is possible to use
    clus_res <- clusterOptimKmedoids(ts_repr, 20, 23, parallel = F)
    # clus_res <- clusterOptimKmeans(ts_repr, 40, 43, parallel = F)
    clusters <- clus_res$clustering
  }
  
  # Create time series for forecasting - mean (K-means) or median (K-medoids)
  # centroids <- t(sapply(sort(unique(clusters)), function(i) colMeans(check.matrix(res_data[clusters %in% i,]))))
  centroids <- t(sapply(sort(unique(clusters)), function(i) apply(check.matrix(res_data[clusters %in% i,]), 2, median)))

  # Forecasting
  centr_forecasts <- ForecastClusterDisagg(centroids, parallel = F, freq = seas, h = h) # set parallel = T, for parallel computation
  
  # Denormalising forecasts
  clust_ind <- lapply(sort(unique(clusters)), function(i) which(clusters == i))
  
  diss_forecasts <- lapply(1:length(centr_forecasts[[1]]), function(l)
    t(sapply(1:nrow(datas), function(i)
      denorm_z(centr_forecasts[[which(sapply(1:length(clust_ind), function(j)
        any(clust_ind[[j]] == i))== TRUE)]][[l]], res_mean[i], res_sd[i]))))
  
  # Benchmark forecasting - train model for every consumer separatelly
  forecasts_original <- t(sapply(1:nrow(datas), function(x) simSnaive(datas[x,], freq = seas, h = h)))
  forecasts_all_ctree <- t(sapply(1:nrow(datas), function(x) simCtree(datas[x,], freq = seas, h = h)))
  
  # save forecasts to files
  write.table(data.frame(i+1, 1:nrow(forecasts_original), forecasts_original), "for_orig_naive_disagg.csv", row.names = F, col.names = F, quote = F, append = T)
  write.table(data.frame(i+1, 1:nrow(forecasts_all_ctree), forecasts_all_ctree), "for_orig_ctree_disagg.csv", row.names = F, col.names = F, quote = F, append = T)
  
  write.table(data.frame(i+1, 1:nrow(diss_forecasts[[1]]), diss_forecasts[[1]]), "for_cent_naive_disagg.csv", row.names = F, col.names = F, quote = F, append = T)
  write.table(data.frame(i+1, 1:nrow(diss_forecasts[[2]]), diss_forecasts[[2]]), "for_cent_mlr_disagg.csv", row.names = F, col.names = F, quote = F, append = T)
  write.table(data.frame(i+1, 1:nrow(diss_forecasts[[3]]), diss_forecasts[[3]]), "for_cent_es_disagg.csv", row.names = F, col.names = F, quote = F, append = T)
  write.table(data.frame(i+1, 1:nrow(diss_forecasts[[4]]), diss_forecasts[[4]]), "for_cent_ctree_disagg.csv", row.names = F, col.names = F, quote = F, append = T)
  
  print(i)
}
