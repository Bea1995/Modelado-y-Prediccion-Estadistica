##################################################
#MASTER UNIVERSITARIO EN MATEMATICAS             #
#Universidad de Sevilla                          # 
#MINERIA ESTADISTICA DE DATOS                    #
#Beatriz Coronado Sanz                           #
#TRABAJO MyPE Tema 3                             #
##################################################

#----------------------- CARGA DEL PAQUETE GSTAT
library(sp)
library(graphics)
library(lattice)
library(xts)
library(gstat)
#------------------------

#----------------------- LECTURA DE DATOS
coal <- read.table("coal.dat", header=TRUE, sep=",")
names(coal)
dim(coal)
#------------------------

#-------- Asignación de coordenadas referencias geográficas
coordinates(coal) = ~X+Y
class(coal)
names(coal)
#------------------------

#---------------- PREPARACIÓN DEL "GRID" REJILLA
coal.grid <- read.table("coalgrid.csv", header=TRUE, sep=";")
names(coal.grid)
dim(coal.grid)
coordinates(coal.grid) = ~X+Y

gridded(coal.grid) = TRUE #Tipo rejilla
class(coal.grid)
names(coal.grid)
#------------------------  

#_____________________________________________________________________________
#  ACCIÓN 1. KRIGING DERIVA EXTERNA
#
# 1.1. CONSTRUCCIÓN DEL VARIOGRAMA (DE LO RESIDUOS) MUESTRAL Y AJUSTE MODELO TEÓRICO.
# 1.2. APLICACIÓN DEL KRIGING DERIVA EXTERNA
# 1.3. COMPARACIÓN ENTRE KRIGING DERIVA EXTERNA Y KRIGING ORDINARIO
# 1.4. COMPARACIÓN UTILIZANDO VALIDACIÓN CRUZADA 
#_________________________________________________________________________________

# 1.1. CONSTRUCCIÓN DEL VARIOGRAMA (DE LO RESIDUOS) MUESTRAL Y AJUSTE MODELO TEÓRICO.
#_________________________________________________________________________________

#------ Calcular el variograma de los residuos del modelo  y compararlo con el original

# Calculamos el variograma muestral
CalorV.vgm = variogram(CalorV~1, coal)
names(CalorV.vgm)
CalorV.vgm

# Calculo del variograma de los residuos del modelo lineal de la variable objetivo
# frente a las coordedanas (grado 1) y Elevat
CalorV.res.vgm = variogram(CalorV~X+Y+Elevat, coal)
CalorV.res.vgm

# Para comprobar que el variograma corresponde a los residuos del modelo lineal, se procede a
# - calcular el modelo lineal
# - obtener directamente el variograma muestral de los residuos
mco.fit <- lm(CalorV~X+Y+Elevat,coal)  # Modelo lineal
names(mco.fit)
summary(mco.fit)

coalres=data.frame(coal,mco.fit$residuals) # Inclusión de los residuos MCO en el data.frame
coalres[1:10,]
names(coalres)
coordinates(coalres) = ~X+Y
res.mco.vgm = variogram(mco.fit.residuals~1, coalres) # Semivariograma muestral delos residuos

# Puede comprobarse la igualdad de ambos variogramas nube
res.mco.vgm
CalorV.res.vgm

# Comparación de los variogramas nube de los residuos con el asociociado a los datos
# originales
comparar.vgm <- data.frame(np = CalorV.vgm$np, dist = CalorV.vgm$dist, 
                           gamma.ok = CalorV.vgm$gamma,   gamma.uk = CalorV.res.vgm$gamma, 
                           gamma.dif = CalorV.vgm$gamma - CalorV.res.vgm$gamma)
comparar.vgm

#Azul: semivariograma muestral de los datos originales
#Verde: semivariograma muestral de los residuos
plot(comparar.vgm$gamma.ok ~ comparar.vgm$dist, pch=20, col="blue", type="b", 
     xlab="Distancia", ylab="Gamma (semivariograma)", 
     ylim=c(0,max(comparar.vgm$gamma.ok, comparar.vgm$gamma.uk)), 
     main = " Variograma, CalorV", sub="OK:azul, UK:verde")
points(comparar.vgm$gamma.uk  ~ comparar.vgm$dist, pch=20, col="green", type="b") 

#Si los dos fueran practicamente igual espuedo utilizar el modelo ordinario. Como se separan,
#utilizamos el de los residuos.

rm(comparar.vgm) #eliminar comparar.vgm

#--------------------------------------
# 	Ajuste a un modelo teórico para el variograma de los residuos 		
#--------------------------------------

attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Sph", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Pen", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Gau", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Cir", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Exp", 600, 0.5)))$SSErr
# Óptimo: Gaussiano

CalorV.res.fit <- fit.variogram(CalorV.res.vgm, model = vgm(6, "Gau", 600, 0.5))
CalorV.res.fit
plot(CalorV.res.vgm, CalorV.res.fit)

#____________________________________________________________________________________
# 1.2. APLICACIÓN DEL KRIGING DERIVA EXTERNA
#_____________________________________________________________________________________

#Kirging ordinario 
#Cogemos el modelo Gaussiano, que ya hemos visto que es el óptimo en clase
CalorV.fit = fit.variogram(CalorV.vgm, model = vgm(8, "Gau", 600, 0.1))
CalorV.kriged = krige(CalorV~1, coal, coal.grid, model = CalorV.fit)

#Kriging deriva externa
CalorV.ukriged = krige(CalorV~X+Y+Elevat, coal, coal.grid, model = CalorV.res.fit)

# Plots de predicciones del Kriging ordinario y el Kriging deriva externa
spplot(CalorV.ukriged["var1.pred"],main="Predicción con kriging deriva externa")
spplot(CalorV.kriged["var1.pred"], main="Predicción con kriging ordinario")

# Plots del Kriging deriva externa
spplot(CalorV.ukriged, zcol="var1.pred", pretty=T, contour=T, col.regions=bpy.colors(64), main="Predicción con kriging deriva externa",
       xlab="X", ylab="Y", scales=list(draw=T))

spplot(CalorV.ukriged, zcol="var1.pred", pretty=T, contour=T, col.regions=bpy.colors(64), main="Predicción con kriging deriva externa", 
       xlab="X", ylab="Y", scales=list(draw=T), cuts=8)
#_____________________________________________________________

#_____________________________________________________________________________________
# 1.3. COMPARACIÓN ENTRE KRIGING DERIVA EXTERNA Y KRIGING ORDINARIO
#_____________________________________________________________________________________

# Cálculo y representación de diferencias entre kriging ordinario y deriva externa  
#	en predicciones y varianzas				   
plot(CalorV.kriged$var1.pred,CalorV.ukriged$var1.pred,xlab="Predición kriging ordinario",
     ylab="Predicción deriva externa")
abline(0,1,col="red")

dif.uk.ok <- data.frame(dif.pred = CalorV.ukriged$var1.pred - CalorV.kriged$var1.pred, 
                        dif.var = CalorV.ukriged$var1.var  - CalorV.kriged$var1.var)
summary(dif.uk.ok)

coordinates(dif.uk.ok) <- coordinates(CalorV.kriged)

spplot(dif.uk.ok, zcol="dif.pred", pretty=T, cuts=12, col.regions = cm.colors(64), main="Predicciones Krig. Der. Ext. - Predicciones Krig. Ordin.") 

spplot(dif.uk.ok, zcol="dif.var", pretty=T, col.regions = cm.colors(64), main="Varianzas Krig. Der. Ext. - Varianzas Krig. Ordin.") 

rm(dif.uk.ok)
#_____________________________________________________________________________________

#_________________________________________________________________________________
# 1.4. COMPARACIÓN UTILIZANDO VALIDACIÓN CRUZADA 
#_________________________________________________________________________________

#---- Estudio de los residuos a través de Validación Cruzada
CalorV.kvc <- krige.cv(CalorV~1, coal, model = CalorV.fit, nfold=96)
CalorV.ukriged.vc = krige.cv(CalorV~X+Y+Elevat, coal, model = CalorV.res.fit,nfold=96)
summary(CalorV.ukriged.vc) 

# Comparación con los residuos obtenidos por validación cruzada en el Kriging Ordinario
summary(CalorV.ukriged.vc$residual) # Residuos CV del Kriging Universal
summary(CalorV.kvc$residual)        # Residuos CV del Kriging Ordinario
summary(CalorV.ukriged.vc$zscore)   # Zscore CV del Kriging Universal
summary(CalorV.kvc$zscore)          # Zscore CV del Kriging Ordinario

###  Histograma de errores estandarizados
par(mfrow = c(2, 1))
hist(CalorV.ukriged.vc$zscore, breaks = seq(-3.5, 3.5, by = 0.5),freq=FALSE,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging Deriva Externa")
hist(CalorV.kvc$zscore, breaks = seq(-3.5, 3.5, by = 0.5), freq=FALSE ,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging Ordinario")
par(mfrow = c(1, 1))
#________________________________________________________________________________________
