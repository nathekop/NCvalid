NCvalid <- function(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100 ,NCstart = TRUE) {
  dm=dim(x)
  d = as.vector(dist(x))
  crr = rep(0,kmax-kmin+2)
  if (NCstart) {
    dtom = sqrt(rowSums((x-colMeans(x))^2))
    crr[1] = sd(dtom)/(max(dtom)-min(dtom))
  }
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
      if (is.null(nrow(x[cluss==j,]))){
        centroid[j,] = as.numeric(x[cluss==j,])
      } else if (nrow(x[cluss==j,])==1){
        centroid[j,] = as.numeric(x[cluss==j,])
      } else {
        centroid[j,] = colMeans(x[cluss==j,])
      }
      xnew[cluss==j,] = rep(centroid[j,],each = sum(cluss==j))
    }
    d2 = as.vector(dist(xnew))
    crr[k-kmin+2]= cor(d,d2,method=corr)
  }
  K = length(crr)
  NWI = ((crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)]))/pmax(0,(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)]))
  NWI2 = (crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)])-(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)])
  NWI3 = NWI
  if (max(NWI)<Inf){
    if (min(NWI) == -Inf){
      NWI3[NWI==-Inf] = min(NWI[is.finite(NWI)])
    }
  }
  if (max(NWI)==Inf){
    NWI3[NWI==Inf] = max(NWI[is.finite(NWI)])+NWI2[NWI==Inf]
    NWI3[WPI<Inf] = NWI[NWI<Inf] + NWI2[NWI<Inf] 
    if (min(NWI) == -Inf){
      NWI3[NWI==-Inf] = min(NWI[is.finite(NWI)])+NWI2[NWI==-Inf]
    }
  }
  
  NWI = data.frame(cbind("k"= kmin:kmax, "NCI1" = NWI))
  NWI2 = data.frame(cbind("k"=kmin:kmax,"NCI2"=NWI2))
  NWI3 = data.frame(cbind("k"=kmin:kmax,"NCI"=NWI3))
  crr = data.frame(cbind("k" = (kmin-1):(kmax+1), "NC" = crr))
  my_list <- list("NC" = crr, "NCI" = NWI3, "NCI1" = NWI, "NCI2" = NWI2)
  return(my_list)
}
