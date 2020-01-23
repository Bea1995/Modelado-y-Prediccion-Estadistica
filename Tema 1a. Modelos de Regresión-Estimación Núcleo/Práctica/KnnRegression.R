#####################################
#MÁSTER UNIVERSITARIO EN MATEMÁTICAS#
#MODELADO Y PREDICCIÓN ESTADÍSTICA  #
#ESTIMACIÓN K-NN                    #
#Beatriz Coronado Sanz              #
#Andrea Prieto García               #
#####################################

#Cargamos los datos de "vino.csv", csv que recoge las variables de diferentes muestras de vino. 
#Dada una muestra, vamos a predecir el grado de alcohol que tiene.
Datos<-read.csv("vino.csv",header=TRUE)
n=nrow(Datos)

#Vamos a dividir los datos que tenemos en Entrenamiento/Test (70%-30%)
indices<- 1:n

inditest= sample(indices,ceiling(n*0.30)) #coge una muestra del 30% de los posibles valores de indices
indient= setdiff(indices,inditest) #coge el 70% restante

#install.packages("caret")
library(caret)

#Entrenamos con el método knn usando cv para obtener el número de vecinos
#Queremos predecir la columna alcohol a partir del conjunto de entrenamiento
entrenam_knn=train(alcohol~.,data=Datos,subset=indient,
                   method="knn",
                   preProcess = c("center", "scale"), #centramos los datos
                   tuneLength = 10,
                   trControl = trainControl(method = "cv")) #cross-validation

#Resumen del modelo calculado
entrenam_knn
plot(entrenam_knn)

#Predecimos el conjunto test con el modelo que hemos calculado
predKNN_test <- predict(entrenam_knn, newdata = Datos[inditest,]) #prediccion de los datos
str(predKNN_test)

#Definimos una función para ver los errores de la predicción y una representación gráfica de los mismos
Ajuste<- function(y,pred,titulo)
{
  residuos=y-pred
  plot(y,pred,main=titulo,ylab="Prediccion indice alcohol", xlab="Valor real indice alcohol",xlim=c(min(y,pred),max(y,pred)),ylim=c(min(y,pred),max(y,pred)))
  abline(a=0,b=1,col="green",lwd=2) #si la prediccion fuese exacta
  grid()
  MSE= mean(residuos^2) 
  RMSE= sqrt(MSE)
  R2= cor(y,pred)^2 #correlacion
  return(list(MSE=MSE,RMSE=RMSE,R2=R2))
}
Ajuste(Datos[inditest,12],predKNN_test,"Test, KNN")
