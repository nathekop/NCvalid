NCvalid <- function(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100) {
  dm=dim(x)
  d = as.vector(dist(x))
  crr = rep(0,kmax-kmin+2)
  if (method == 'hclust_complete') {
    hh = hclust(dist(x),method = 'complete')
  } else if (method == 'hclust_average') {
    hh = hclust(dist(x),method = 'average')
  }  else if (method == 'hclust_single') {
    hh = hclust(dist(x),method = 'single')
  }
  if (kmin == 2){
    lb = 2
  } else {
    lb = kmin-1
  }
  for (k in lb:(kmax+1))
  {
    if (method == 'kmeans'){
      km.out=kmeans(x,k,nstart =nstart)
      cluss = km.out$cluster
    }
    else {
      cluss = cutree(hh,k)
    }

    xnew = matrix(rep(0,dm[1]*dm[2]),dm[1],dm[2])
    centroid = matrix(1:(dm[2]*k),k,dm[2])
    for (j in 1:k)
    {
      if (nrow(x[cluss==j,])>1){
        centroid[j,] = colMeans(x[cluss==j,])
      }
      else {
        centroid[j,] = x[cluss==j,]
      }
      xnew[cluss==j,] = rep(centroid[j,],each = sum(cluss==j))
    }
    d2 = as.vector(dist(xnew))
    crr[k-kmin+2]= cor(d,d2,method=corr)
  }
  K = length(crr)
  NWI = ((crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)]))/((crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)]))
  NWI2 = (crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)])-(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)])
  NWI = data.frame(cbind("k"= kmin:kmax, "NCI1" = NWI))
  NWI2 = data.frame(cbind("k"=kmin:kmax,"NCI2"=NWI2))
  crr = data.frame(cbind("k" = (kmin-1):(kmax+1), "NC" = crr))
  my_list <- list("NC" = crr, "NCI1" = NWI, "NCI2" = NWI2)
  return(my_list)
}
