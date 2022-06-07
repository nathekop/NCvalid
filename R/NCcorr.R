NCcorr <- function(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100,NCstart=TRUE) {
  dm=dim(x)
  crr = rep(0,kmax-kmin+1)
  dis = dist(x)
  if (NCstart) {
    dtom = sqrt(rowSums((x-colMeans(x))^2))
    sd(dtom)/(max(dtom)-min(dtom))
  }
  if (method == 'hclust_complete') {
    hh = hclust(dist(x),method = 'complete')
  }
  else if (method == 'hclust_average') {
    hh = hclust(dist(x),method = 'average')
  }
  else if (method == 'hclust_single') {
    hh = hclust(dist(x),method = 'single')
  }
  for (k in kmin:kmax)
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
      if (length(x[cluss==j,])>dm[2]){
        centroid[j,] = colMeans(x[cluss==j,])
      }
      else {
        centroid[j,] = x[cluss==j,]
      }
      xnew[cluss==j,] = rep(centroid[j,],each = sum(cluss==j))
    }
    d = as.vector(dist(x))
    d2 = as.vector(dist(xnew))
    crr[k-kmin+1]= cor(d,d2,method=corr)
  }
  crr = cbind("k"=(kmin:kmax),"NC" = crr)
  return(crr)
}
