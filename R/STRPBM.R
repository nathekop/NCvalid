STRPBM <- function(x,kmax,kmin=2,method='kmeans',nstart=100) {
  dm=dim(x)
  STR = rep(0,kmax-kmin+1)
  PBM = rep(0,kmax-kmin+1)
  EK = rep(0,kmax-kmin+3)
  DK = rep(0,kmax-kmin+3)
  md = rep(0,kmax-kmin+3)
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
    
    
    centroid = matrix(1:(dm[2]*k),k,dm[2])
    for (j in 1:k)
    {
      if (is.null(nrow(x[cluss==j,]))){
        centroid[j,] = as.numeric(x[cluss==j,])
        dtovj = 0
      } else if (nrow(x[cluss==j,])==1){
        centroid[j,] = as.numeric(x[cluss==j,])
        dtovj = 0
      } else {
        centroid[j,] = colMeans(x[cluss==j,])
        dtovj = apply(x[cluss==j,],1,function(x)sqrt(sum((x-centroid[j,])^2)))
      }
     
      EK[k-kmin+2] = EK[k] + sum(dtovj)
    }
    ddd = dist(centroid)
    md[k-kmin+2] = max(ddd)
    DK[k-kmin+2] = max(ddd)/min(ddd)
  }
  E0 = sum(sqrt(rowSums((x-colMeans(x))^2)))
  if (kmin == 2){
    EK[1] = E0
  }
  EKK = E0/EK
  STR = (EKK[2:(length(EKK)-1)]-EKK[1:(length(EKK)-2)])*(DK[3:(length(DK))]-DK[2:(length(DK)-1)])
  STR = data.frame(cbind("k"= kmin:kmax, "STR" = STR))
  
  PBM = EKK[2:(length(EKK)-1)]*md[2:(length(EKK)-1)]/(kmin:kmax)
  PBM = data.frame(cbind("k"= kmin:kmax, "PBM" = PBM))
  
  my_list = list("STR"=STR,"PBM" = PBM)
  return(my_list)
}
