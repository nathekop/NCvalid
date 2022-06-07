NCI2 <- function(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100,NCstart=TRUE) {
  nw = NCvalid(x,kmax,kmin,method,corr,nstart)
  nw = nw$NCI2
}
