##################################################
#MASTER UNIVERSITARIO EN MATEMATICAS             #
#Universidad de Sevilla                          # 
#MINERIA ESTADISTICA DE DATOS                    #
#Beatriz Coronado Sanz                           #
#TRABAJO MyPE Tema 2                             #
##################################################

# Carga del paquete survival
library(survival)
set.seed(51511)

# Algunas opciones de los c�lculos
options(digits = 4) # n�mero de decimales en las salidas
options(columns = 40) # numero de columnas en las salidas

## Carga de datos
datos = read.csv("uissurv.csv",header=TRUE, sep=";")  # read csv file 

#Muestra del conjunto de entrenamiento
n<- nrow(datos);
nent<- ceiling(0.95*n);
indient=sample(n,nent);

datos=datos[indient,]

dim(datos)
names(datos)
summary(datos)

#Modelo de Cox
time.cox2<- coxph(Surv(time, censor) ~ age + ivhx + ndrugtx + race + site, 
                  data=datos)
summary(time.cox2)

#Interpretaci�n de los coeficientes

# Sobre la funci�n de riesgo h(t;x), el aumento de 1 unidad de la variable,
# provoca un cambio en dicha funci�n de:  h(t;x+1)=exp(coef) h(t;x)

#1. Para la variable edad
coefficients(time.cox2)[1]
exp(coefficients(time.cox2)[1])

#Observamos que el valor de exp(coef-Edad)=0.9736, lo que quiere decir
#que una persona un a�o m�s mayor que otra hace disminuir un poco el 
#valor de la funci�n de riesgo

#2. Para la variable ivhx
coefficients(time.cox2)[2:3]
exp(coefficients(time.cox2)[2:3])

#En este caso, la variable ivhx se define con dos variables dummy 
#(ivhxPrevious) e ivhxRecent) sobre la variable referencia ivhxNever. 
#Los coeficientes #indican el cambio que se produce en la funci�n de 
#riesgo frente a una persona que nunca tomo medicamentos IV.

#Haber tomado medicamentos IV recientemente implica que la funci�n de 
#riesgo se multiplica por exp(coef_Recent)=1.391 frente a no haber 
#tomado medicamentos IV nunca.

#Haber tomado medicamentos IV previamente implica que la funci�n de 
#riesgo se multiplica por exp(coef_Previous)=1.164 frente a no haber 
#tomado medicamentos IV nunca.

#Por lo tanto, haber tomado medicamentos IV hace poco tiene m�s riesgo
#de recaer en las drogas que haber tomado medicamentos IV previamente, 
# lo que tiene m�s riesgo de recaen en las drogas que no haber tomado 
#medicamentos IV nunca.

#Test de hip�tesis de riesgos proporcionales con los residuos de Schoenfeld
cox.zph(time.cox2,transform = "rank")
# Ninguna variable viola la hip�tesis de riesgos proporcionales.

#Comportamiento gr�fico para la variable edad
datos.fin <- with(datos, data.frame(
  id=c(1001,1002), age  = c(32,33),     
  beck = rep(mean(datos$beck,na.rm=TRUE), 2), 
  ndrugtx = rep(mean(datos$ndrugtx,na.rm=TRUE), 2), 
  hercoc  = c("Heroin & Cocain", "Heroin & Cocain"), 
  ivhx = c("Previous","Previous"), 
  race = c("White","White"), 
  treat = c("Long","Long"), site = c("B", "B"), los = rep(mean(datos$los), 2)))
datos.fin

# Plot de la funci�n de supervivencia
plot(survfit(time.cox2, newdata=datos.fin), conf.int=FALSE,
     lty=c(1, 2), ylim=c(0.0, 1), xlab="semanas",
     ylab="Proporci�n de no reincidentes",col = c(2,3),
     main="Estimaci�n de la supervivencia")
legend("topright", legend=c("Age = 32", "Age = 33"), lty=c(1 ,2), col = c(2,3),inset=0.02)

# Plot de la funci�n de riesgo acumulado
plot(survfit(time.cox2, newdata=datos.fin),conf.int=FALSE,
     lty=c(1, 2), 
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),lwd=2,fun="cumhaz", col = c(2,3),
     main="Estimaci�n de la tasa de fallo acumulada") 
legend("bottomright", legend=c("Age = 32", "Age = 33"), lty=c(1 ,2), col = c(2,3),inset=0.02)
abline(h=0) 
grid()

#Comportamiento gr�fico para la variable ivhx
datos.fin <- with(datos, data.frame(
  id=c(1001,1002,1003), age  = rep(mean(datos$age,na.rm=TRUE), 3),     
  beck = rep(mean(datos$beck,na.rm=TRUE), 3), 
  ndrugtx = rep(mean(datos$ndrugtx,na.rm=TRUE), 3), 
  hercoc  = c("Heroin & Cocain", "Heroin & Cocain","Heroin & Cocain"), 
  ivhx = c("Never","Previous","Recent"), 
  race = c("White","White","White"), 
  treat = c("Long","Long","Long"), site = c("B", "B","B"), los = rep(mean(datos$los), 3)))
datos.fin

# Plot de la funci�n de supervivencia
plot(survfit(time.cox2, newdata=datos.fin), conf.int=FALSE,
     lty=c(1, 2,3), ylim=c(0.0, 1), xlab="semanas",
     ylab="Proporci�n de no reincidentes",col = c(2,3,4),
     main="Estimaci�n de la supervivencia")
legend("topright", legend=c("Ivhx = never", "Ivhx = previous","Ihvx=recent"), lty=c(1 ,2,3), col = c(2,3,4),inset=0.02)

# Plot de la funci�n de riesgo acumulado
plot(survfit(time.cox2, newdata=datos.fin),conf.int=FALSE,
     lty=c(1, 2,3), 
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),lwd=2,fun="cumhaz", col = c(2,3,4),
     main="Estimaci�n de la tasa de fallo acumulada") 
legend("bottomright", legend=c("Ivhx = never", "Ivhx = previous","Ihvx=recent"), lty=c(1 ,2,3), col = c(2,3,4),inset=0.02)
abline(h=0) 
grid()
