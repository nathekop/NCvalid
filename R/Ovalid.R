Ovalid <- function(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100,indexlist='all',NCstart = TRUE) {
  dm=dim(x)
  dis = dist(x)
  PB = rep(0,kmax-kmin+1)
  CSL = rep(0,kmax-kmin+1)
  CH = rep(0,kmax-kmin+1)
  GD33 = rep(0,kmax-kmin+1)
  GD43 = rep(0,kmax-kmin+1)
  GD53 = rep(0,kmax-kmin+1)
  SC = rep(0,kmax-kmin+1)
  DB = rep(0,kmax-kmin+1)
  DBs = rep(0,kmax-kmin+1)
  SF = rep(0,kmax-kmin+1)
  DI = rep(0,kmax-kmin+1)
  STR = rep(0,kmax-kmin+1)
  PBM = rep(0,kmax-kmin+1)
  crr = rep(0,kmax-kmin+2)
  NCI1 = rep(0,kmax-kmin+1)
  NCI2 = rep(0,kmax-kmin+1)
  NCI3 = rep(0,kmax-kmin+1)
  if (method == 'hclust_complete') {
    hh = hclust(dist(x),method = 'complete')
  } else if (method == 'hclust_average') {
    hh = hclust(dist(x),method = 'average')
  } else if (method == 'hclust_single') {
    hh = hclust(dist(x),method = 'single')
  }
  for (k in kmin:kmax){
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

    clusdist = NULL
    for (j in 1:k)
    {
      dd = as.matrix(dist(x[cluss==j,]))
      CSL[k-kmin+1] = CSL[k-kmin+1] + sum(apply(dd,2,max))/sum(cluss==j)
      if (sum(cluss==j) > 1){
        clusdist = rbind(clusdist,xnew[cluss==j,][1,])
      } else {
        clusdist = rbind(clusdist,xnew[cluss==j,])
      }

    }
    dd2 = as.matrix(dist(clusdist))
    dd2 = matrix(dd2[dd2 > 0],k-1,k)
    CSL[k-kmin+1] = CSL[k-kmin+1]/sum(apply(dd2,2,min))

    if (sum(indexlist == "all") ==1 | "CH"%in% indexlist | "GDI33" %in% indexlist | "GDI43" %in% indexlist | "GDI53" %in% indexlist){
      intall = intCriteria(as.matrix(x),cluss,c("Calinski_Harabasz","GDI33","GDI43","GDI53"))
      CH[k-kmin+1] = intall$calinski_harabasz
      GD33[k-kmin+1] = intall$gdi33
      GD43[k-kmin+1] = intall$gdi43
      GD53[k-kmin+1] = intall$gdi53
    }

    if (sum(indexlist == "all") ==1 | "SC"%in% indexlist){
      ss = silhouette(cluss, dis)
      SC[k-kmin+1]  = mean(ss[,3])
    }

    if (sum(indexlist == "all") ==1 | "DB"%in% indexlist){
      dd = index.DB(x, cluss)
      DB[k-kmin+1] = dd$DB
    }

    if (sum(indexlist == "all") ==1 | "DI"%in% indexlist){
      DI[k-kmin+1] = dunn(distance = dis, cluss, method = "euclidean")
    }
    
    if (sum(indexlist == "all") ==1 | "SF"%in% indexlist | "DBs"%in% indexlist ){
      ddb = rep(0,k)
      wcdd = rep(0,k)
      sizecluss = rep(0,k)
      for (jj in 1:k){
        sizecluss[jj] = sum(cluss==jj)
      }
      for (j in 1:k){
        ddb[j] = max(dd$S[j]+dd$S[-j])/min(dist(rbind(centroid[j,],centroid[-j,]))[1:k-1])
        wcdd[j] = sum(dist(rbind(centroid[j,],x[cluss==j,]))[1:sizecluss[j]])/sizecluss[j]
      }
      DBs[k-kmin+1] = sum(ddb)/k
      bcd = sum(dist(rbind(colMeans(x),centroid))[1:k]*sizecluss)/(dim(x)[1]*k)
      wcd = sum(wcdd)
      SF[k-kmin+1] = 1-1/exp(bcd+wcd)
    }

    d = as.vector(dis)
    d3 = as.vector(dist(xnew))
    d3[d3>0] = 1
    PB[k-kmin+1] = cor(d,d3,method=corr)
  }
  CSL = data.frame(cbind("k"=kmin:kmax,"CSL"=CSL))
  GD33 = data.frame(cbind("k"=kmin:kmax,"GD33"=GD33))
  GD43 = data.frame(cbind("k"=kmin:kmax,"GD43"=GD43))
  GD53 = data.frame(cbind("k"=kmin:kmax,"GD53"=GD53))
  SC = data.frame(cbind("k"=kmin:kmax,"SC"=SC))
  DB = data.frame(cbind("k"=kmin:kmax,"DB"=DB))
  DBs = data.frame(cbind("k"=kmin:kmax,"DBs"=DBs))
  SF = data.frame(cbind("k"=kmin:kmax,"SF"=SF))
  DI = data.frame(cbind("k"=kmin:kmax,"DI"=DI))
  CH = data.frame(cbind("k"=kmin:kmax,"CH"=CH))
  PB = data.frame(cbind("k"=kmin:kmax,"CH"=PB))

if (sum(indexlist == 'all')==1 | 'NC' %in% indexlist | 'NCI1' %in% indexlist | 'NCI2' %in% indexlist  | 'NCI' %in% indexlist){
  nw = NCvalid(x,kmax,kmin,method,corr,nstart,NCstart)
  crr = nw$NC
  NCI1 = nw$NCI1
  NCI2 = nw$NCI2
  NCI = nw$NCI
}
if (sum(indexlist == "all") ==1 | "STR"%in% indexlist | "PBM"%in% indexlist){
  sss = STRPBM(x,kmax,kmin,method,nstart)
  STR = sss$STR
  PBM = sss$PBM
}

  my_list <- list("NC"=crr, "NCI" = NCI, "NCI1" = NCI1, "NCI2" = NCI2, "CH" = CH, "CSL"=CSL, "DB"=DB, "DBs"=DBs, "DI"=DI, "GD33" = GD33, "GD43" = GD43, "GD53" = GD53, "PB"=PB, "PBM"=PBM, "SF"=SF,  "SC"=SC, "STR"=STR)
  if (sum(indexlist == "all")==1){
    return(my_list)
  } else {
    return(my_list[indexlist])
  }
}









