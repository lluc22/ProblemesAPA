n <- 5
X <- matrix(c(rep(1,n), seq(n)),nrow=n)

t <- c( 7.97,10.2,14.2,16.0,21.2)

plot(X[,2],t,lwd=3)
dev.copy(png,'data.png')
dev.off()


#Comprovem rank de la matriu:
library(Matrix)
rankMatrix(X)

#és full rank per tant mínims qudrats té solució única
 
# 2. Solució amb pseudo-inversa
#w = (X^T X)^-1 X^T *t
(w <- solve(t(X) %*% X) %*% t(X) %*% t)
lines (X[,2], w[2,1]*X[,2]+w[1,1], type="l")
dev.copy(png,'model.png')
dev.off()
# 3. Solució amb svd

s <- svd(X)
#Podem reafirmar que el rank és dos gràcies a svd.


D <- diag(1/s$d)
(w <- s$v %*% D %*% t(s$u) %*% t)

#Dues solucions iguals:
 #w0 = e = 4.235
 #w1 = k = 3.226
