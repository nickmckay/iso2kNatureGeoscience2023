iso2kPCA <- function(data, recordName, data_summary){
  
  data <- data
  recordName <- recordName
  data_summary <- data_summary
  dat_din <- dineof(data, delta.rms = 1e-04)$Xa
  pca_out <- rda(dat_din, scale=TRUE)
  inertia <- pca_out$CA$eig/pca_out$tot.chi
  pca_sites <- scores(pca_out, choices=1:4)$sites
  pca_spp <- data.frame(scores(pca_out, choices=1:4)$species)
  pca_sts1 <- pca_sites[,1]  
  pca_sts2 <- pca_sites[,2]
  pca_sts3 <- pca_sites[,3]
  pca_sps1 <- pca_spp[,1]
  pca_sps2 <- pca_spp[,2]
  pca_sps3 <- pca_spp[,3]
  
  pca_spp <- cbind(pca_spp, recordName)
  PCA_results <<- merge(data_summary, pca_spp, by="recordName")
  PCA_timeseries <<- pca_sites
  data_DIN <<- dat_din
  inertia <<- inertia
  pca_out <<- pca_out
}