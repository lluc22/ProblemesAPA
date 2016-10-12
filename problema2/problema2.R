library(mlbench)
library (cclust)
library(Rmixmod)
library(ggplot2)

N <- 2000

data.1 <- mlbench.cassini(N,relsize = c(1,1,0.25))


plot(data.1,main="Original Data with its true classification")
plot(data.1$x[,1],data.1$x[,2],main="Original Data without classification",xlab="",ylab="")



kmeans.3 <- cclust (data.1$x,3,iter.max=100,method="kmeans",dist="euclidean")
ch.3 <- clustIndex(kmeans.3,data.1$x, index="calinski")
for (i in 1:1000) {
  r <- cclust (data.1$x,3,iter.max=100,method="kmeans",dist="euclidean")
  s <- clustIndex(kmeans.3,data.1$x, index="calinski")
  if(s > ch.3){
    ch.3 <- s
    kmeans.3 <- r
  }
}


plot(data.1$x[,1],data.1$x[,2],col=(kmeans.3$cluster+1),main="K-means with k=3",xlab="",ylab="")
points(kmeans.3$centers,col=seq(1:kmeans.3$ncenters)+1,cex=2,pch=19)


#Apartat 2
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


plot(data.1$x[,1],data.1$x[,2],col=(kmeans.max$cluster+1),xlab="",ylab="",main="K-means with k=2")
points(kmeans.max$centers,col=seq(1:kmeans.max$ncenters)+1,cex=2,pch=19)


fammodel <- mixmodGaussianModel (family="diagonal", equal.proportions=FALSE)


z <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = 3)

summary(z)

means <- z@bestResult@parameters@mean
found.clusters <- z@bestResult@partition

plot(data.1$x[,1],data.1$x[,2],col=(found.clusters+1),xlab="",ylab="",main="E-M with k=3")


ks <- c(2,5,10,50)




par(mfrow=c(2,2))

for (k in ks) {
  
  z <- mixmodCluster (data.frame(data.1$x),models = fammodel, nbCluster = k)
  summary(z)
  means <- z@bestResult@parameters@mean
  found.clusters <- z@bestResult@partition
  plot(data.1$x[,1],data.1$x[,2],col=(found.clusters+1),xlab="",ylab="",main=paste(c("E-M with k =",k),sep=" "))
  
}

