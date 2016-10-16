library(mlbench)
library (cclust)
library(Rmixmod)
library(ggplot2)
library(clusterCrit)
rm(list = ls())
par(mfrow=c(1,1))

N <- 2000

data.1 <- mlbench.cassini(N,relsize = c(1,1,0.25))
dataM <- matrix(c(data.1$x[,1],data.1$x[,2]),byrow=F,ncol=2)

plot(data.1,main="Original Data with its true classification")
plot(data.1$x[,1],data.1$x[,2],main="Original Data without classification",xlab="",ylab="")


#apartat 2
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


#Apartat 3
do.kmeans <- function (whatK)
{
  r <- cclust (data.1$x,whatK,iter.max=100,method="kmeans",dist="euclidean")
  (clustIndex(r,data.1$x, index="calinski"))
}


res <- vector("numeric",10)
for (K in 2:10)
  res[K] <- max (r <- replicate (20, do.kmeans(K)))

(Kmax <- which.max(res))
plot(res, type="l",main="K-means performance",xlab="k",ylab="CH")


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

print(Kmax)
plot(data.1$x[,1],data.1$x[,2],col=(kmeans.max$cluster+1),xlab="",ylab="",main=paste("K-means with k=",Kmax))
points(kmeans.max$centers,col=seq(1:kmeans.max$ncenters)+1,cex=2,pch=19)

#apartat4

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
  fammodel <- mixmodGaussianModel (family="diagonal", equal.proportions=FALSE)
  z <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = k,criterion="BIC")
  return(CH(x=dataM,centr=z@bestResult@parameters@mean,part=z@bestResult@partition))
}

resEM <- vector("numeric",10)
for (K in 2:10)
  resEM[K] <- max (r <- replicate (10, do.EM(K)))
(KmaxEM <- which.max(resEM))
print(KmaxEM)
plot(resEM, type="l",main="EM Performance",xlab="k",ylab="CH",xlim=c(2,10))


fammodel <- mixmodGaussianModel (family="diagonal", equal.proportions=FALSE)
EM.max <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = KmaxEM)
chEM <- CH(dataM,EM.max@bestResult@parameters@mean,EM.max@bestResult@partition)
for (i in 1:100) {
  r <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = KmaxEM)
  s <- CH(dataM,EM.max@bestResult@parameters@mean,EM.max@bestResult@partition)
  if(s > chEM){
    chEM <- s
    EM.max <- r
  }
}
plot(EM.max)
print(chEM)
print(ch.max)

EMCrit <- intCriteria(dataM,as.vector(EM.max@bestResult@partition,"integer"),"all")
KmeansCrit <- intCriteria(dataM,as.vector(kmeans.max$cluster,"integer"),"all")

