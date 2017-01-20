library(datasets)
library(nnet)

#Executem per una llavor en concret, per tenir els mateixos resultats.
set.seed(26212)

data(rock)

N <- nrow(rock)
rock.x <- data.frame(area = scale(rock$area), peri = scale(rock$peri),shape = scale(rock$shape))
rock.y <- log(rock$perm)
rocks <- data.frame(perm=rock.y,rock.x)

#Set de train
learn <- sample(1:N, round(3*N/4))

nlearn <- length(learn)
#set de test
ntest <- N - nlearn


#Calcula l'error quadràtic normalitzat
norm.mse <- function(x,pred){
  sum((x - pred)^2) / ((length(x)-1)*var(x))
}


library(caret)

#Probem sense regularitzar amb size 4
model.nnet <- nnet(perm ~., data = rocks, subset=learn, size=4 ,maxit=200, decay=0,linout=TRUE,trace=F)

#Error de training
norm.mse(rocks$perm[learn],predict(model.nnet))

#Error de test
testPred <- predict(model.nnet, newdata=rocks[-learn,]) 
norm.mse(rocks$perm[-learn],testPred)

#Ens quedem amb size 4 i regularitzem, provem amb diferents valors de lambda.
(decays <- 10^seq(-5,0,by=0.1))

trc <- trainControl (method="LOOCV")

model.LOOCV <- train (perm ~., data = rocks, subset=learn, method='nnet', maxit = 500, trace = FALSE,
                        tuneGrid = expand.grid(.size=4,.decay=decays), trControl=trc,linout=TRUE)

model.LOOCV$bestTune
lambda <- model.LOOCV$bestTune$decay


best <- nnet(perm ~., data = rocks, subset=learn, size=4 ,maxit=200, decay=lambda,linout=TRUE,trace=F)

#Calculem error
norm.mse(rocks$perm[learn],predict(best))
testPred <- predict(best, newdata=rocks[-learn,]) 
norm.mse(rocks$perm[-learn],testPred)

#Representem valor observat vs. predit
rockVsPred <- data.frame(perm=rocks$perm,pred=predict(model.nnet,newdata = rocks))
ggplot(rockVsPred,aes(perm,pred)) +
  geom_point() + 
  geom_abline(intercept=0,slope=1,colour='red') +
  xlab("Observed Value") +
  ylab("Predicted Value")

