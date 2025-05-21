packages <- c("data.table", "R.matlab","DSL","clusterability","MASS","spatstat.geom")
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
  library(pkg, character.only = TRUE)
})


library(data.table)
library(R.matlab)
library(DSL)
library(clusterability)
library(MASS)
library(spatstat.geom)


Allfeatures <- function(data){
  data = data.frame(data)
  distm = dist(data)
  d = dim(data)
  First_PC_SD = PCAVarComp1(data)
  SNR = SignalToNoiseRatio(data)
  DBSCAN_Clusters =SpreadMeasure3(data)
  OPO_Res_GDeg_PO_Median = Gdeg(data)
  Network_Density = net_dens(data)
  
  avg_nnd = avg_nnd(data)
  
  TempMF1 = MF1(distm)
  
  Dist_Skew = TempMF1[1]
  Dist_Kurtosis = TempMF1[2]
  dist_NRE = TempMF1[3]
  
  SilvermanStatistic = clusterability::clusterabilitytest(data,'silverman',reduction='none')$critbw
  
  
  return(c(First_PC_SD,SNR,DBSCAN_Clusters,OPO_Res_GDeg_PO_Median,Network_Density,avg_nnd,Dist_Skew,Dist_Kurtosis,dist_NRE,SilvermanStatistic))
  
}

avg_nnd <- function(ds) {
  #ds <- ds[, 1:(dim(ds)[2]-1)]
  dst <- enn(ds, 0.15*nrow(ds))
  graph <- igraph::graph.adjacency(dst, mode = "undirected", weighted = TRUE)
  hs <- igraph::knn(graph)
  return(round(mean(hs$knn), 4))
}

#dataWl = produceData(2,2,20,1)

Gdeg <- function(dat){
  dat <- as.data.frame(dat)
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  relations <- FNN::knn.index(dat,pref.k)
  vert <- 1:dim(dat)[1]
  #tdatt[,c(2,1)]
  g <- igraph::graph_from_data_frame(relations, directed=TRUE, vertices=vert)
  deg.vals <- igraph::degree(g)
  pot.outliers <- order(deg.vals,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  gdegpo = RatiosBetweenTwoGroups(deg.vals, pot.outliers)
  return(gdegpo)
}

RatiosBetweenTwoGroups <-  function(dx, grp){
  ## dx is some 1D quantity like, degree, number of components, residuals, density etc
  ## grp is the group that it belongs to like outliers, pot.outliers
  
  F2 <- ifelse(median(dx[-grp],na.rm = TRUE)==0, 0, median(dx[grp],na.rm = TRUE)/median(dx[-grp],na.rm = TRUE))
  
  out <- c(F2)
  return(out)
}


PCAVarComp1 <- function(dat){
  pca.dat <- prcomp(dat, center = TRUE, scale=TRUE )
  firstpc.sd <- pca.dat$sdev[1]/sum(pca.dat$sdev)
  return(firstpc.sd)
}

SignalToNoiseRatio <- function(dat){
  mean.all <- apply(dat,2, mean)
  sd.all <- apply(dat,2,sd)
  return(mean(mean.all/sd.all))
}

SpreadMeasure3 <- function(dat){
  dat <- data.frame(dat)
  nn <- max(floor(dim(dat)[1]/200),10)
  knn.out <- FNN::knn.dist(dat,nn)
  epsilon <- quantile(knn.out[,2],probs =0.9)
  dbscan.ex <- dbscan::dbscan(dat, eps=epsilon)
  output <- length(unique(dbscan.ex$cluster))-1
  return(output)
}

net_dens <- function(ds) {
  #ds <- ds[, 1:(dim(ds)[2]-1)]
  dst <- enn(ds, 0.15*nrow(ds))
  graph <- igraph::graph.adjacency(dst, mode = "undirected", weighted = TRUE)
  density <- igraph::graph.density(graph)
  return(round(density, 4))
}

FindNonBinaryVariables <- function(dat){
  un.vals <- apply(dat,2,function(x) length(unique(x)))
  cols <- which(un.vals >2)
  return(cols)
}

enn <- function(ds, e) {
  #ds <- ds[, 1:(dim(ds)[2]-1)]
  ds = ds
  dst <- as.matrix(dist(ds))
  for(i in 1:nrow(ds)) {
    a <- names(sort(dst[i,])[1:(e+1)])
    b <- rownames(ds)
    dst[i, setdiff(rownames(ds), intersect(a, b))] <- 0
  }
  return(dst)
}

ResidualsOfProxisAndOutliers <- function(dat){
  dat <- as.data.frame(dat)
  
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  z.n <- FNN::knn.dist(dat,pref.k)
  z1 <- z.n[,pref.k]
  pot.outliers <- order(z1,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  
  cols <- FindNonBinaryVariables(dat)
  if(length(cols)==0){  ## all columns are binary
    cols = 1:dim(dat)[2]
  }
  
  if(length(cols)==1){  ## only one numeric column in the dataset
    cols = 1:dim(dat)[2]
  }
  
  if(dim(dat)[1]/length(cols)< 10){  ## less than 10 observations per attribute
    num.cols.to.be.chosen <- min(floor(dim(dat)[1]/10), 50)
    vars <- apply(dat[,cols],2,var)
    cols <- order(vars, decreasing=TRUE)[1:num.cols.to.be.chosen]
  }
  
  if(length(cols)>1){
    set.seed(101)
    cols <- sample(cols,length(cols))
  }
  
  out <- data.frame(Res_KNOut_SD=numeric())
  
  
  
  if(length(cols)>1){
    kk <- min(length(cols),50)
    output.mat <- matrix(0, nrow=kk, ncol=35)
    for(ii in 1:kk){
      y <- dat[cols[ii]]
      x <- dat[cols[-ii]]
      model <- lm(unlist(y)~., data=x)
      pot.outliers.residuals <- order(abs(model$residuals),decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
      output.mat[ii,17] <- length(intersect(pot.outliers,pot.outliers.residuals))/length(pot.outliers)
    }
    #out <- apply(output.mat,2,mean)
    out <- mean(output.mat)
  }else{  ## only one column anyway
    residuals <- dat[,cols]-mean(dat[,cols])
    pot.outliers.residuals <- order(abs(residuals),decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
    out3 <- length(intersect(pot.outliers,pot.outliers.residuals))/length(pot.outliers)
    out <- out3
  }
  return(out)
}


MF1 <- function(distm){
  m = (distm-min(distm))/(max(distm)-min(distm))
  mf4 = moments::skewness(m)
  mf5 = moments::kurtosis(m)
  mf6 = sum(m <= 0.1)
  mf7 = sum(m > 0.1 & m <= 0.2)
  mf8 = sum(m > 0.2 & m <= 0.3)
  mf9 = sum(m > 0.3 & m <= 0.4)
  mf10 = sum(m > 0.4 & m <= 0.5)
  mf11 = sum(m > 0.5 & m <= 0.6)
  mf12 = sum(m > 0.6 & m <= 0.7)
  mf13 = sum(m > 0.7 & m <= 0.8)
  mf14 = sum(m > 0.8 & m <= 0.9)
  mf15 = sum(m > 0.9 & m <= 1)
  NRE = NA
  NRE = try(philentropy::KL(rbind(c(mf6,mf7,mf8,mf9,mf10,mf11,mf12,mf13,mf14,mf15)/length(m),rep(1/10,10))))
  return(c(mf4,mf5,NRE))
  
}



#{'Num_Observations'}    {'First_PC_SD'}    {'SNR'}    {'Spread_0'}    {'OPO_Res_KNOut_SD'}    {'OPO_GDeg_PO_Med…'}    {'net_dens'}    {'dist_skew'}    {'dist_kurtosis'}  {'silvN'}


normaliseBOXCOXZ <- function(feats,prelim){
  xmin = prelim[,1]
  feats = feats - xmin + 1
  lam =prelim[,2]
  bx = c()
  mu = prelim[,3]
  sig = prelim[,4]
  zc = c()
  for(i in 1:length(feats)){
    if(feats[i]<0){feats[i]=0}
    if(lam[i] != 0 ){
      bx[i] = (feats[i]^lam[i]-1)/lam[i]
    }
    else{
      bx[i] = log(feats[i])
    }
    zc[i] = (bx[i] - mu[i])/sig[i]
  }
  return(zc)
}

