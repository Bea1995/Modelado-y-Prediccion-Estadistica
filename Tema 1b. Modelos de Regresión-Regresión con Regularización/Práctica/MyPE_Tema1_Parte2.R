##################################################
#MASTER UNIVERSITARIO EN MATEMATICAS             #
#Universidad de Sevilla                          # 
#MODELADO Y PREDICCIÓN ESTADÍSTICA               #
#Beatriz Coronado Sanz                           #
#TRABAJO REGRESIÓN CON REGULARIZACIÓN            #
##################################################

set.seed("13579")

#Carga del conjunto de datos
datos=read.csv("dataworkMASTER.csv",sep=";",header=TRUE) #read csv file
dim(datos)
names(datos)
str(datos)
datos=datos[,-1] #Eliminación de la variable código Cod_Id

#Eliminación 
#DNI=(51511971->1971)
a=1
b=9
c=7
d=1
#La variable 1 es "varobj" y las siguientes xi para i entre [1,40]
elim=c(11+a,21+a,31+a,11+b,21+b,31+b,11+c,31+c,31+d);elim
elim2=sample(unique(elim),4,replace=FALSE);elim2

datos=datos[,-elim2]

#Apartado 1
#Modelo de regresión multiple de la variable "varobj" frente al resto
datos.mco <- lm(varobj ~ ., datos)
summary(datos.mco)

#Analizar el problema de multicolinealidad

#Metodo 1. Con el determinandte de la matriz de covarianzas
#Como el determinante es casi 0 hay multicolinalidad
det(vcov(datos.mco))

#Metodo 2. Con el factor inflactor de la varianza (VIF)
# Si VIF es alto (superior a 10) se dice que hay multicolinealidad
library(usdm)
#VIF de todas las variables
vif(datos[,-1])
#Vemos las que tienen multicolinealidad
vif(datos[,-1])[vif(datos[,-1])[2]>10,]

#Apartado 2
#Preparar los datos en formato matriz
mx=as.matrix(datos[,2:37])
my=as.matrix(datos[,1])

library(glmnet)

#Técnica de regularización elasticnet
solElastic=glmnet(mx,my,alpha=0.5,lambda=seq(0, 0.1, 0.001))
summary(solElastic)

# Matriz de coeficientes estimados (para cada valor de lambda por filas)
# La primera columna almacena el valor de lambda
head(coef(solElastic)[,1:3])
head(coef(solElastic)[,99:101])

#Regularización LASSO
solLasso=glmnet(mx,my,alpha=1,lambda=seq(0, 0.1, 0.001))
summary(solLasso)

# Matriz de coeficientes estimados (para cada valor de lambda por filas)
# La primera columna almacena el valor de lambda
head(coef(solElastic)[,1:3])
head(coef(solElastic)[,99:101])

#Regularización RIDGE
solRidge=glmnet(mx,my,alpha=0,lambda=seq(0, 0.1, 0.001))
summary(solRidge)

# Matriz de coeficientes estimados (para cada valor de lambda por filas)
# La primera columna almacena el valor de lambda
head(coef(solRidge)[,1:3])
head(coef(solRidge)[,99:101])

#Apartado 3
#Análisis comparativo entre los cuatro modelos

#Creamos un conjunto de entrenamiento y un conjunto de test
train=sample(1:nrow(mx),nrow(mx)/2)
test=(-train)
ytest=my[test]

#MSE regresión múltiple
regmul <- lm(varobj ~ ., datos[train,])
summary(regmul)

regmul.pred=predict(regmul,newx=mx[test,])
MSEregMult=mean((regmul.pred-ytest)^2); MSEregMult

#MSE rigde
rigde2=glmnet(mx[train,],my[train,],alpha=0,lambda=seq(0,0.1,0.001))
plot(rigde2)

#Predicción para lambda muy alto
rigde.pred=predict(rigde2,s=1e10,newx=mx[test,])
mean((rigde.pred-ytest)^2)

#Predicción para lambda óptimo (usando validación cruzada)
cv.outRG=cv.glmnet(mx[train,],my[train],alpha=0,lambda=seq(0,0.1,0.001))
plot(cv.outRG)
bestlamRG=cv.outRG$lambda.min; bestlamRG

rigde.pred=predict(rigde2,s=bestlamRG,newx=mx[test,])
MSEridge=mean((rigde.pred-ytest)^2); MSEridge

predict(solRidge,type="coefficients",s=bestlamRG)

#MSE lasso
lasso2=glmnet(mx[train,],my[train,],alpha=1,lambda=seq(0,0.1,0.001))
plot(lasso2)

#Predicción para lambda muy alto
lasso.pred=predict(lasso2,s=1e10,newx=mx[test,])
mean((lasso.pred-ytest)^2)

#Predicción para lambda óptimo (usando validación cruzada)
cv.outLS=cv.glmnet(mx[train,],my[train],alpha=0,lambda=seq(0,0.1,0.001))
plot(cv.outLS)
bestlamLS=cv.outLS$lambda.min; bestlamLS

lasso.pred=predict(lasso2,s=bestlamLS,newx=mx[test,])
MSElasso=mean((lasso.pred-ytest)^2); MSElasso

predict(solLasso,type="coefficients",s=bestlamLS)

#MSE elastic
elastic2=glmnet(mx[train,],my[train,],alpha=0.5,lambda=seq(0,0.1,0.001))
plot(elastic2)

#Predicción para lambda muy alto
elastic.pred=predict(elastic2,s=1e10,newx=mx[test,])
mean((elastic.pred-ytest)^2)

#Predicción para lambda óptimo (usando validación cruzada)
cv.outEL=cv.glmnet(mx[train,],my[train],alpha=0,lambda=seq(0,0.1,0.001))
plot(cv.outEL)
bestlamEL=cv.outEL$lambda.min; bestlamEL

elastic.pred=predict(elastic2,s=bestlamEL,newx=mx[test,])
MSEelastic=mean((elastic.pred-ytest)^2); MSEelastic

predict(solElastic,type="coefficients",s=bestlamEL)

#Comparacion de los metodos utilizados
Metodo=c("Ridge", "Lasso", "Elastic", "RegMult")
Lambda=c(bestlamRG,bestlamLS,bestlamEL,"")
MSE=c(MSEridge,MSElasso,MSEelastic,MSEregMult)

errores=cbind(Metodo,MSE,Lambda);errores
