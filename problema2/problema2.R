library(mlbench)
library (cclust)
library(Rmixmod)
library(ggplot2)
rm(list = ls())
par(mfrow=c(1,1))
printf <- function(...) invisible(print(sprintf(...)))

N <- 2000

data.1 <- mlbench.cassini(N,relsize = c(1,1,0.25))
dataM <- matrix(c(data.1$x[,1],data.1$x[,2]),byrow=F,ncol=2)

plot(data.1,main="Original Data with its true classification")
plot(data.1$x[,1],data.1$x[,2],main="Original Data without classification",xlab="",ylab="")


#apartat 2
#Trobar el millor k-means amb k=3 basant-nos amb l'índex CH

print("Apartat 2")

kmeans.3 <- cclust (data.1$x,3,iter.max=100,method="kmeans",dist="euclidean")
ch.3 <- clustIndex(kmeans.3,data.1$x, index="calinski")
for (i in 1:100) {
  r <- cclust (data.1$x,3,iter.max=100,method="kmeans",dist="euclidean")
  s <- clustIndex(kmeans.3,data.1$x, index="calinski")
  if(s > ch.3){
    ch.3 <- s
    kmeans.3 <- r
  }
}


plot(data.1$x[,1],data.1$x[,2],col=(kmeans.3$cluster+1),main="K-means with k=3",xlab="",ylab="")
points(kmeans.3$centers,col=seq(1:kmeans.3$ncenters)+1,cex=2,pch=19)

printf("El CH del millor 3-means és: %f",ch.3)


#Apartat 3
#Trobar quina és la millor K per a k-means

print("Apartat 3")
do.kmeans <- function (whatK)
{
  r <- cclust (data.1$x,whatK,iter.max=100,method="kmeans",dist="euclidean")
  (clustIndex(r,data.1$x, index="calinski"))
}


res <- vector("numeric",10)
for (K in 2:10)
  res[K] <- max (r <- replicate (20, do.kmeans(K)))

Kmax <- which.max(res)
printf("La millor K per a k-means és: %d",Kmax)
plot(res, type="l",main="K-means performance",xlab="k",ylab="CH")


#Trobem el millor k-means amb k=kmax basant-nos en CH
kmeans.max <- cclust (data.1$x,Kmax,iter.max=100,method="kmeans",dist="euclidean")
ch.max <- clustIndex(kmeans.max,data.1$x, index="calinski")
for (i in 1:100) {
  r <- cclust (data.1$x,Kmax,iter.max=100,method="kmeans",dist="euclidean")
  s <- clustIndex(kmeans.max,data.1$x, index="calinski")
  if(s > ch.max){
    ch.max <- s
    kmeans.max <- r
  }
}

#El representem
plot(data.1$x[,1],data.1$x[,2],col=(kmeans.max$cluster+1),xlab="",ylab="",main=paste("K-means with k=",Kmax))
points(kmeans.max$centers,col=seq(1:kmeans.max$ncenters)+1,cex=2,pch=19)

#apartat4
#Provem EM amb diferents k's per veure quin és el millor, amb diferents criteris

print("Apartat 4")

#Intentem trobar la millor K segons l'index CH

CH <- function(x,centr,part){
  # C-H = (SSB/(K-1)) / (SSW/(N-K))
  SSB <- 0
  SSW <- 0
  N <- length(x[,1])
  K <- length(centr[,1])
  d <- ncol(x)
  m = rep(0,d)
  for (i in d){
    m[i] = mean(x[,i])
  }
  for (k in 1:K) {
    for (i in which(part == k)) {
      SSW = SSW + norm(matrix(x[i,] - centr[k,]),'F')^2
    }
  }
  for (k in 1:K) {
    Nk <- length(which(part == k))
    SSB = SSB + Nk *norm(matrix(centr[k,] - m),'F')^2 
  }
  return ( (SSB/(K-1)) / (SSW/(N-K)) )
  
}


do.EM <- function(k){
  fammodel <- mixmodGaussianModel (family="general", equal.proportions=FALSE)
  z <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = k)
  return(CH(x=dataM,centr=z@bestResult@parameters@mean,part=z@bestResult@partition))
}

#Provem per totes les k's
resEM <- vector("numeric",5)
for (K in 2:5)
  resEM[K] <- max (r <- replicate (10, do.EM(K)))
KmaxEM <- which.max(resEM)
plot(resEM, type="l",main="EM Performance",xlab="k",ylab="CH",xlim=c(2,5))
printf("La millor K segons CH és: %d",KmaxEM)


#Ens quedem amb la millor execució
fammodel <- mixmodGaussianModel (family="general", equal.proportions=FALSE)
EM.max <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = KmaxEM)
chEM <- CH(dataM,EM.max@bestResult@parameters@mean,EM.max@bestResult@partition)
for (i in 1:10) {
  r <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = KmaxEM)
  s <- CH(dataM,EM.max@bestResult@parameters@mean,EM.max@bestResult@partition)
  if(s > chEM){
    chEM <- s
    EM.max <- r
  }
}
plot(EM.max)

#Ara provem amb l'index BIC, com més baix millor
do.EM.BIC <- function(k){
  fammodel <- mixmodGaussianModel (family="general", equal.proportions=FALSE)
  z <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = k)
  return(z@bestResult@criterionValue)
}

#Probem per una selecció de k's
resEM.BIC <- vector("numeric",5)
for (K in 2:5)
  resEM.BIC[K] <- min (r <- replicate (10, do.EM.BIC(K)))
resEM.BIC[1] = Inf
Kmin <- which.min(resEM.BIC)
plot(resEM.BIC, type="l",main="EM Performance",xlab="k",ylab="BIC",xlim=c(2,5))
printf("La millor k segons BIC és: %d", Kmin)

#Ens quedem amb la millor execució
EM.min <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = Kmin)
BIC <- CH(dataM,EM.min@bestResult@parameters@mean,EM.min@bestResult@partition)
for (i in 1:20) {
  r <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = Kmin)
  s <- CH(dataM,EM.min@bestResult@parameters@mean,EM.min@bestResult@partition)
  if(s < BIC){
    BIC <- s
    EM.min <- r
  }
}
plot(EM.min)


#Finalment provem de trobar la k que minimitza el nombre d'elements 
#que no "encaixen" amb la partició original (La veritat)

clusterDiff <- function(originalPartition,partition){
  n1 <- length(which(originalPartition == 1))
  n2 <- length(which(originalPartition == 2))
  n3 <- length(which(originalPartition == 3))
  max1 <- max(table(partition[which(originalPartition == 1)]))
  max2 <- max(table(partition[which(originalPartition == 2)]))
  max3 <- max(table(partition[which(originalPartition == 3)]))
  
  return((n1-max1) + (n2 - max2) + (n3 - max3))
}


do.EM.Diff <- function(k){
  fammodel <- mixmodGaussianModel (family="diagonal", equal.proportions=FALSE)
  z <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = k)
  return(clusterDiff(data.1$classes,z@bestResult@partition))
}

#Provem amb totes les k's
resEM.Diff <- vector("numeric",5)
for (K in 2:5)
  resEM.Diff[K] <- min (r <- replicate (10, do.EM.Diff(K)))
resEM.Diff[1] = Inf
KminDiff <- which.min(resEM.Diff)
plot(resEM.Diff, type="l",main="EM Performance",xlab="k",ylab="Cluster Diff",xlim=c(2,5))
printf("La k que fa que el clustering s'assembli més al clustering original és: %d",KminDiff)


#Ens quedem amb la millor execució
fammodel <- mixmodGaussianModel (family="diagonal", equal.proportions=FALSE) #canviem a diagonal, ja que ara assumim que sabem la veritat sobre les dades.
EM.Diff <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = KminDiff)
Diff <- clusterDiff(data.1$classes,EM.Diff@bestResult@partition)
for (i in 1:100) {
  r <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = KminDiff)
  s <- clusterDiff(data.1$classes,EM.Diff@bestResult@partition)
  if(s < Diff){
    Diff <- s
    EM.Diff <- r
  }
}
plot(EM.Diff)
printf("El nombre d'elements que són diferents entre el millor clustering obtingut i l'original és de: %d",Diff)
