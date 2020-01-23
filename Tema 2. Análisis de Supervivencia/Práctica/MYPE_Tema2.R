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

# Algunas opciones de los cálculos
options(digits = 4) # número de decimales en las salidas
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

#Interpretación de los coeficientes

# Sobre la función de riesgo h(t;x), el aumento de 1 unidad de la variable,
# provoca un cambio en dicha función de:  h(t;x+1)=exp(coef) h(t;x)

#1. Para la variable edad
coefficients(time.cox2)[1]
exp(coefficients(time.cox2)[1])

#Observamos que el valor de exp(coef-Edad)=0.9736, lo que quiere decir
#que una persona un año más mayor que otra hace disminuir un poco el 
#valor de la función de riesgo

#2. Para la variable ivhx
coefficients(time.cox2)[2:3]
exp(coefficients(time.cox2)[2:3])

#En este caso, la variable ivhx se define con dos variables dummy 
#(ivhxPrevious) e ivhxRecent) sobre la variable referencia ivhxNever. 
#Los coeficientes #indican el cambio que se produce en la función de 
#riesgo frente a una persona que nunca tomo medicamentos IV.

#Haber tomado medicamentos IV recientemente implica que la función de 
#riesgo se multiplica por exp(coef_Recent)=1.391 frente a no haber 
#tomado medicamentos IV nunca.

#Haber tomado medicamentos IV previamente implica que la función de 
#riesgo se multiplica por exp(coef_Previous)=1.164 frente a no haber 
#tomado medicamentos IV nunca.

#Por lo tanto, haber tomado medicamentos IV hace poco tiene más riesgo
#de recaer en las drogas que haber tomado medicamentos IV previamente, 
# lo que tiene más riesgo de recaen en las drogas que no haber tomado 
#medicamentos IV nunca.

#Test de hipótesis de riesgos proporcionales con los residuos de Schoenfeld
cox.zph(time.cox2,transform = "rank")
# Ninguna variable viola la hipótesis de riesgos proporcionales.

#Comportamiento gráfico para la variable edad
datos.fin <- with(datos, data.frame(
  id=c(1001,1002), age  = c(32,33),     
  beck = rep(mean(datos$beck,na.rm=TRUE), 2), 
  ndrugtx = rep(mean(datos$ndrugtx,na.rm=TRUE), 2), 
  hercoc  = c("Heroin & Cocain", "Heroin & Cocain"), 
  ivhx = c("Previous","Previous"), 
  race = c("White","White"), 
  treat = c("Long","Long"), site = c("B", "B"), los = rep(mean(datos$los), 2)))
datos.fin

# Plot de la función de supervivencia
plot(survfit(time.cox2, newdata=datos.fin), conf.int=FALSE,
     lty=c(1, 2), ylim=c(0.0, 1), xlab="semanas",
     ylab="Proporción de no reincidentes",col = c(2,3),
     main="Estimación de la supervivencia")
legend("topright", legend=c("Age = 32", "Age = 33"), lty=c(1 ,2), col = c(2,3),inset=0.02)

# Plot de la función de riesgo acumulado
plot(survfit(time.cox2, newdata=datos.fin),conf.int=FALSE,
     lty=c(1, 2), 
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),lwd=2,fun="cumhaz", col = c(2,3),
     main="Estimación de la tasa de fallo acumulada") 
legend("bottomright", legend=c("Age = 32", "Age = 33"), lty=c(1 ,2), col = c(2,3),inset=0.02)
abline(h=0) 
grid()

#Comportamiento gráfico para la variable ivhx
datos.fin <- with(datos, data.frame(
  id=c(1001,1002,1003), age  = rep(mean(datos$age,na.rm=TRUE), 3),     
  beck = rep(mean(datos$beck,na.rm=TRUE), 3), 
  ndrugtx = rep(mean(datos$ndrugtx,na.rm=TRUE), 3), 
  hercoc  = c("Heroin & Cocain", "Heroin & Cocain","Heroin & Cocain"), 
  ivhx = c("Never","Previous","Recent"), 
  race = c("White","White","White"), 
  treat = c("Long","Long","Long"), site = c("B", "B","B"), los = rep(mean(datos$los), 3)))
datos.fin

# Plot de la función de supervivencia
plot(survfit(time.cox2, newdata=datos.fin), conf.int=FALSE,
     lty=c(1, 2,3), ylim=c(0.0, 1), xlab="semanas",
     ylab="Proporción de no reincidentes",col = c(2,3,4),
     main="Estimación de la supervivencia")
legend("topright", legend=c("Ivhx = never", "Ivhx = previous","Ihvx=recent"), lty=c(1 ,2,3), col = c(2,3,4),inset=0.02)

# Plot de la función de riesgo acumulado
plot(survfit(time.cox2, newdata=datos.fin),conf.int=FALSE,
     lty=c(1, 2,3), 
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),lwd=2,fun="cumhaz", col = c(2,3,4),
     main="Estimación de la tasa de fallo acumulada") 
legend("bottomright", legend=c("Ivhx = never", "Ivhx = previous","Ihvx=recent"), lty=c(1 ,2,3), col = c(2,3,4),inset=0.02)
abline(h=0) 
grid()
