#__________________________________________________________________________________
#                       Modelado y Predicci�n Estad�stica.                                      M�ster en Matem�ticas. Universidad de Sevilla                ###
#                           Juan Manuel Mu�oz Pichardo                        
#__________________________________________________________________________________
#__________________________________________________________________________________
#                 AN�LISIS ESTRUCTURAL, KRIGING ORDINARIO,      
#                 KRIGING UNIVERSAL, KRIGING CON DERIVA EXTERNA   
#                 Conjunto de datos: Coal + Coal.grid             
#__________________________________________________________________________________
   
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

#-------- NOTA SOBRE EL FICHERO DE DATOS coal.dat 
#
#  Un conjunto simulado de datos basado en una veta de carb�n real 
#  en el Sudeste de �frica. En los pozos perforados en la veta de carb�n se 
#  mide las variables: espesor de carb�n (metros), el contenido energ�tico o 
#  "valor calor�fico 'del carb�n (Megajulios por tonelada), el contenido de 
#  cenizas (%) y el contenido de azufre (%). 
#  Adem�s se consideran las coordenadas de la localizaci�n (X,Y) en metros
#  as� como la elevaci�n (en metros).
#------------------------

#-------- Asignaci�n de coordenadas referencias geogr�ficas
coordinates(coal) = ~X+Y
class(coal)
names(coal)
#------------------------

#-------- RESUMEN Y REPRESENTACION GR�FICA (BURBUJAS)
summary(coal)

# histogramas
hist(coal$Thickness, col = "lightblue",border = "red", main = "Espesor de carb�n")
rug(coal$Thickness)
hist(coal$CalorV, col = "lightblue",border = "red", main = "Valor Calor�fico")
rug(coal$CalorV)

# Plots de burbujas
bubble(coal, c("Thickness"),col=c("#00ff0088", "#00ff0088"), main = "Espesor del carbon")
bubble(coal, c("CalorV"),col=c("green","blue"), main = "Valor calor�fico")
bubble(coal, c("Ash"),  col=c("green","red"), main = "Porcentaje de cenizas")
bubble(coal, c("Sulphur"),  col=c("green","red"), main = "Porcentaje de azufre")
#------------------------

#---------------- PREPARACI�N DEL "GRID" REJILLA
coal.grid <- read.table("coalgrid.csv", header=TRUE, sep=";")
names(coal.grid)
dim(coal.grid)
coordinates(coal.grid) = ~X+Y

gridded(coal.grid) = TRUE #Tipo rejilla
class(coal.grid)
names(coal.grid)

#-------------   NOTA SOBRE EL FICHERO DE DATOS coal.grid 
# Las observaciones est�n georreferenciadas en coordenadas UTM (x e y)
# En cada localizaci�n se han recogido:
#   - x     : Coordenada X 
#   - y     : Coordenada Y 
#------------------------  

spplot(coal.grid, c("Elevat"), col.regions=terrain.colors(25),
       main="Elevat", scales=list(draw=TRUE),xlab="X",ylab="Y")

#_____________________________________________________________________________
#
# ACCI�N 1. AN�LISIS ESTRUCTURAL, KRIGING ORDINARIO
# 1.1. CONSTRUCCI�N DEL VARIOGRAMA MUESTRAL
# 1.2. AJUSTE DE MODELOS TEORICOS DE VARIOGRAMA
# 1.3. KRIGING ORDINARIO
# 1.4. M�TODO DE VALIDACI�N CRUZADA SAOBRE EL KRIGING ORDINARIO
# 1.5. PREDICCI�N SOBRE UN PUNTO ESPACIAL
#_____________________________________________________________________________

#____________________________________________________________________________
#1.1. CONSTRUCCI�N DEL VARIOGRAMA MUESTRAL
#____________________________________________________________________________

#------------- Funci�n "variogram"
### Calcula el variograma muestral, con opciones para variograma direccional 
### y robusto.
### 
### variogram(object, data, locations, cutoff, Cressie,...)
###    object = objeto de la clase "gstat". Se indica "variable ~ regresor1+regresor2+..."         
###             En caso de ausencia de regresores, utilizar "variable ~ 1"
###    data   = Conjunto de datos
###    cloud  = (l�gico), si es TRUE, c�lculo del demivariograma nube
###    locations = ubicaciones de datos espaciales
###    cutoff = valor de distancia de separaci�n m�xima para incluir los 
###             pares  de puntos en las estimaciones  de la semivarianza;
###             Por defecto, 1/3 de la longitud de la diagonal del marco que envuelve los datos 
###   Cressie = (l�gico) TRUE (Estimaci�n robusta de Cressie del variograma); 
###                      FALSE (m�todo cl�sico de los momentos)
#----------------------------------

#Variograma nube
CalorV.cloud <- variogram(CalorV ~ 1, coal, cloud = T) #~1: frente a nada
head(CalorV.cloud)
plot(CalorV.cloud,col="blue",main="Semivariograma nube CalorV") 

#Variograma muestral
CalorV.vgm = variogram(CalorV~1, coal)
CalorV.vgm
plot(CalorV.vgm, col="black",main="Semivariograma experim. CalorV")

CalorV.vgm

#____________________________________________________________________________
# 1.2. AJUSTE DE MODELOS TEORICOS DE VARIOGRAMA
#____________________________________________________________________________

#------------- Funci�n "vgm": vgm(s, "Sph", r, n)
###
### Genera un modelo de variograma seg�n el modelo Esf�rico (Sph), con 
###       umbral=s; rango=r; nugget=n
###
### Haciendo la llamada "vgm()" relaciona los distintos modelos
###    Modelos= Nug (nugget); Exp (exponential); Sph (spherical); Gau (gaussian)
###             Cir (circular); Lin (linear); Pen (pentaspherical);  Pow (power); ...
###
###                        Funci�n "fit.variogram"
### 
### Ajusta rango y/o sill a partir de un modelo de variograma simple o anidado 
### a un variograma muestral.
###     fit.variogram(object, model, fit.method = 7)
###
###         object= variograma muestral, resultado de la funci�n "variogram"
###         model = modelo variograma, resultado de la funci�n "vgm")
###         La funci�n tiene un  argumento opcional ("fit.method"), que especifica la ponderaci�n
###         de los puntos del variograma emp�rico para el ajuste de m�nimos cuadrados. 
###         Por defecto, "fit.method=7" pesos proporcionales a la cantidad de pares de puntos e 
###         inversamente proporcional al cuadrado de la distancia de separaci�n: $N_h/h^2$. 
###          
###         1 : pesos=$N_h$ ; 2 : pesos=$N_h/(gamma(h))^2$ (Propuesto por Cressie)
###         5 : M�xima verosimilitud restringida ; 6 : sin pesos (MCO).
#---------------------------------------------------------------------- 


#-------------- BONDAD DE AJUSTE DEL MODELO TEOR�CO DE VARIOGRAMA 
# Suma de cuadrados de los residuos, ponderados seg�n el modelo ajustado
# Se desea realizar diversos ajustes para elegir el m�s adecuado

# Ajuste con pesos proporcionales a la cantidad de pares de puntos e 
# inversamente proporcional al cuadrado de la distancia de separaci�n: $N_h/h^2$.

attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Sph", 700, 0.1)))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Pen", 700, 0.1)))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Gau", 700, 0.1)))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Cir", 700, 0.1)))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Exp", 700, 0.1)))$SSErr

# Ajuste con pesos=$N_h/(gamma(h))^2$ (Propuesto por Cressie)
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Sph", 700, 0.1),fit.method=2))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Pen", 700, 0.1),fit.method=2))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Gau", 700, 0.1),fit.method=2))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Cir", 700, 0.1),fit.method=2))$SSErr
attributes(fit.variogram(CalorV.vgm, model=vgm(5, "Exp", 700, 0.1),fit.method=2))$SSErr

# Seleccionamos el modelo Gaussiano
CalorV.fit = fit.variogram(CalorV.vgm, model = vgm(8, "Gau", 600, 0.1))
CalorV2.fit = fit.variogram(CalorV.vgm, model = vgm(8, "Gau", 600, 0.1),fit.method=2)
plot(CalorV.vgm, CalorV.fit,main="Ajuste variograma")
plot(CalorV.vgm, CalorV2.fit,main="Ajuste variograma (Cressie)")
CalorV.fit
CalorV2.fit

# Pr�cticamente coinciden los variogramas estimados. Mantenemos el primero CalorV.fit 
print(plot(CalorV.vgm, plot.numbers = F, pch =18, col = "darkblue", model = CalorV.fit))

#______________________________________________________________________________________
# 1.3. KRIGING ORDINARIO
#______________________________________________________________________________________

#--------------- Funci�n "krige" 
### Realiza las predicciones kriging sobre los puntos definidos en el grid
###
###             Orden: krige(variable~1, datos, datos.grid, model = modelo.fit)
###
###  datos.grid : conjunto de datos sobre el que se desea predecir
###  modelo.fit : modelo te�rico de variograma ajustado a trav�s de fit.variogram
###  
###  Genera un objeto con dos variables: 
###            "var1.pred" = Predicciones
###					   "var1.var" = Varianzas de las predicciones
#--------------------------------------------------------

CalorV.kriged = krige(CalorV~1, coal, coal.grid, model = CalorV.fit)
names(CalorV.kriged)
CalorV.kriged$var1.pred[1:5] # Predicciones de los cinco primeros casos
CalorV.kriged$var1.var[1:5]  # Varianzas de las predicciones de los cinco primeros casos

head(CalorV.kriged@data)

#------------- 	REPRESENTACIONES GR�FICAS

# Plot espacial de la predicci�n con graduaci�n de colores 
spplot(CalorV.kriged["var1.pred"],main="Plot espacial de la predicci�n" )

spplot(CalorV.kriged["var1.var"],main="Plot espacial de la predicci�n" )

# Plot espacial de la predicci�n y la varianza de la predicci�n con graduaci�n de colores
spplot(CalorV.kriged, main="Plots espaciales de la predicci�n y la varianza de la predicci�n" )

#---------
# NOTA:   En los spplots se puede cambiar la gama de colores a trav�s 
#         del argumento: col.regions=....
#---------

# Plot espacial de la predicci�n con graduaci�n de colores y l�neas de contorno 
spplot(CalorV.kriged, zcol="var1.pred", pretty=T, contour=T, col.regions=bpy.colors(64), 
       main="Predicci�n por kriging ordinario sobre Valor Calor�fico", 
       xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

# Plot espacial de la varianza de la estimaci�n con graduaci�n de colores y l�neas de contorno
spplot(CalorV.kriged, zcol="var1.var", pretty=T, contour=T, col.regions=terrain.colors(64),
       main="Varianza de la predicci�n por kriging ordinario sobre Valor Calor�fico", 
       xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))


# Plot espacial de la varianza de la estimaci�n con graduaci�n de colores sin l�neas de contorno
spplot(CalorV.kriged, zcol="var1.var", pretty=T, contour=F, col.regions=bpy.colors(10), 
       main="Varianza de la predicci�n por kriging ordinario sobre Valor Calor�fico", 
       xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

#		Representaci�n de curva de nivel de las predicciones con los puntos observados
contour(CalorV.kriged)
points(coordinates(coal))
#____________________________________________________________________________________

#____________________________________________________________________________________
# 1.4. M�TODO DE VALIDACI�N CRUZADA SOBRE EL KRIGING ORDINARIO
#____________________________________________________________________________________

##----------- Funci�n krige.cv
##
##		El argumento principal es "nfold=k", con k entero mayor que 1
## 		El procedimiento validaci�n cruzada de k iteraciones o K-fold cross-validation 
## 		los datos de muestra se dividen en k subconjuntos. Uno de los subconjuntos se utiliza 
## 		como datos de prueba y el resto (k-1) como datos de entrenamiento. El proceso se repite 
## 		k iteraciones, con cada uno de los posibles subconjuntos de datos de prueba. 
## 		Finalmente se realiza la media aritm�tica de los resultados de cada iteraci�n para obtener 
## 		un �nico resultado.  
## 		Para k=n�mero de observaciones, aplica el m�todo de validaci�n cruzada dejando uno fuera
#-----------------------------------------------------------
dim(coal)
CalorV.kvc <- krige.cv(CalorV~1, coal, model = CalorV.fit, nfold=96)

# 	Esta orden crea un objeto con las siguientes variables
#  		var1.pred (predicciones); var1.var (varianza de la predicci�n)
#  		observed (valores observados) ; residual (residuos)
#  		zscore (residuos estandarizados; fold (n�mero de grupo o submuestra)
dim(CalorV.kvc)
names(CalorV.kvc)

#----------   Criterios para evaluar el ajuste del variograma a trav�s 
#             de la validaci�n cruzada
#
#-- Medidas:  Medias de los residuos

summary(CalorV.kvc)
mean(CalorV.kvc$residual)
mean(CalorV.kvc$zscore^2)

#-- Representaciones gr�ficas

# Mapa de la ubicaci�n de los puntos con marcas 
# en funci�n del tama�o de los errores estandarizados
bubble(CalorV.kvc, "residual", main = "Cross Validaci�n: Residuos")
bubble(CalorV.kvc, "zscore"  , main = "Cross Validaci�n: Residuos Estandarizados")

#  Diagrama de dispersi�n de los valores estimados frente a los observados
plot(CalorV.kvc$observed,CalorV.kvc$var1.pred)
abline(0,1,col="red")

# Diagrama de dispersi�n de los valores estimados 
# frente a los residuos estandarizados
plot(CalorV.kvc$var1.pred, CalorV.kvc$zscore)
abline(0,0,col="red")

#  Histograma de errores estandarizados
hist(CalorV.kvc$zscore, breaks = seq(-3.5, 3.5, by = 0.5), 
     col = "lightblue",border = "red", main = "Errores estandarizados")
#_____________________________________________________________________________

#_____________________________________________________________________________
# 1.5. PREDICCI�N SOBRE UN PUNTO ESPACIAL
#_____________________________________________________________________________

#  Predicci�n sobre el punto de coordenadas X=10500, Y=13500 
punto <- SpatialPoints(data.frame(X=10500, Y=13500))

punto.kriged <- krige(CalorV~1, coal, punto, model = CalorV.fit)
names(punto.kriged)
punto.kriged
#_________________________________________________________________________________

#_____________________________________________________________________________
#  ACCI�N 2. KRIGING UNIVERSAL 
#
# 2.1. CONSTRUCCI�N DEL VARIOGRAMA (DE LO RESIDUOS) MUESTRAL Y AJUSTE MODELO TE�RICO.
# 2.2. KRIGING UNIVERSAL
# 2.3. COMPARACI�N ENTRE KRIGING UNIVERSAL Y KRIGING ORDINARIO
# 2.4. COMPARACI�N UTILIZANDO VALIDACI�N CRUZADA 
#_________________________________________________________________________________

# 2.1. CONSTRUCCI�N DEL VARIOGRAMA (DE LO RESIDUOS) MUESTRAL Y AJUSTE MODELO TE�RICO.
#_________________________________________________________________________________

#------  Calcular el variograma de los residuos del modelo  y compararlo con el original
#----------------------------------------

# El variograma muestral de CalorV era: "CalorV.vgm
names(CalorV.vgm)
CalorV.vgm

# Calculo del variograma de los residuos del modelo lineal de la variable objetivo
# frente a las coordedanas (grado 1)
CalorV.res.vgm = variogram(CalorV~X+Y, coal)
CalorV.res.vgm


# Para comprobar que el variograma corresponde a los residuos del modelo lineal, se procede a
# - calcular el modelo lineal
# - obtener directamente el variograma muestral de los residuos
mco.fit <- lm(CalorV~X+Y,coal)  # Modelo lineal
names(mco.fit)
summary(mco.fit)

coalres=data.frame(coal,mco.fit$residuals) # Inclusi�n de los residuos MCO en el data.frame
coalres[1:10,]
names(coalres)
coordinates(coalres) = ~X+Y
res.mco.vgm = variogram(mco.fit.residuals~1, coalres) # Semivariograma muestral delos residuos

# Puede comprobarse la igualdad de ambos variogramas nube
res.mco.vgm
CalorV.res.vgm

# Comparaci�n de los variogramas nube de los residuos con el asociociado a los datos
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

rm(comparar.vgm) #eliminar comparar.vgm

#--------------------------------------
# 	Ajuste a un modelo te�rico para el variograma de los residuos 		
#--------------------------------------

attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Sph", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Pen", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Gau", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Cir", 600, 0.5)))$SSErr
attributes(fit.variogram(CalorV.res.vgm, model=vgm(6, "Exp", 600, 0.5)))$SSErr
# �ptimo: Gaussiano

CalorV.res.fit <- fit.variogram(CalorV.res.vgm, model = vgm(6, "Gau", 600, 0.5))
CalorV.res.fit
plot(CalorV.res.vgm, CalorV.res.fit)

#--------------  APLICACI�N DEL KRIGING UNIVERSAL
#--------------------------------------------------

CalorV.ukriged = krige(CalorV~X+Y, coal, coal.grid, model = CalorV.res.fit)

# Plots de predicciones del Kriging ordinario y el Kriging Universal
spplot(CalorV.ukriged["var1.pred"],main="Predicci�n con kriging universal")
spplot(CalorV.kriged["var1.pred"], main="Predicci�n con kriging ordinario")

spplot(CalorV.ukriged, zcol="var1.pred", pretty=T, contour=T, col.regions=bpy.colors(64), main="Predicci�n con kriging universal",
       xlab="X", ylab="Y", scales=list(draw=T))

spplot(CalorV.ukriged, zcol="var1.pred", pretty=T, contour=T, col.regions=bpy.colors(64), main="Predicci�n con kriging universal", 
       xlab="X", ylab="Y", scales=list(draw=T), cuts=8)
#_____________________________________________________________

#_____________________________________________________________________________________
# 2.3. COMPARACI�N ENTRE KRIGING UNIVERSAL Y KRIGING ORDINARIO
#_____________________________________________________________________________________

# 	C�lculo y representaci�n de diferencias entre kriging ordinario y universal  
#		en predicciones y varianzas				   

plot(CalorV.kriged$var1.pred,CalorV.ukriged$var1.pred,xlab="Predici�n kriging ordinario",
     ylab="Predicci�n kriging universal")
abline(0,1,col="red")

dif.uk.ok <- data.frame(dif.pred = CalorV.ukriged$var1.pred - CalorV.kriged$var1.pred, 
                         dif.var = CalorV.ukriged$var1.var  - CalorV.kriged$var1.var)
summary(dif.uk.ok)

coordinates(dif.uk.ok) <- coordinates(CalorV.kriged)

spplot(dif.uk.ok, zcol="dif.pred", pretty=T, cuts=12, col.regions = cm.colors(64), main="Predicciones Krig. Univ. - Predicciones Krig. Ordin.") 

spplot(dif.uk.ok, zcol="dif.var", pretty=T, col.regions = cm.colors(64), main="Varianzas Krig. Univ. - Varianzas Krig. Ordin.") 

rm(dif.uk.ok)
#_____________________________________________________________________________________

#_________________________________________________________________________________
# 2.4. COMPARACI�N UTILIZANDO VALIDACI�N CRUZADA 
#_________________________________________________________________________________

#---- Estudio de los residuos a trav�s de Validaci�n Cruzada
#------------------------------------------------------------
CalorV.ukriged.vc = krige.cv(CalorV~X+Y, coal, model = CalorV.res.fit,nfold=96)
summary(CalorV.ukriged.vc) 

# Comparaci�n con los residuos obtenidos por validaci�n cruzada en el Kriging Ordinario
summary(CalorV.ukriged.vc$residual) # Residuos CV del Kriging Universal
summary(CalorV.kvc$residual)        # Residuos CV del Kriging Ordinario
summary(CalorV.ukriged.vc$zscore)   # Zscore CV del Kriging Universal
summary(CalorV.kvc$zscore)          # Zscore CV del Kriging Ordinario

###  Histograma de errores estandarizados
par(mfrow = c(2, 1))
hist(CalorV.ukriged.vc$zscore, breaks = seq(-3.5, 3.5, by = 0.5),freq=FALSE,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging Universal")
hist(CalorV.kvc$zscore, breaks = seq(-3.5, 3.5, by = 0.5), freq=FALSE ,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging Ordinario")
par(mfrow = c(1, 1))
#________________________________________________________________________________________

#________________________________________________________________________________________
# ACCI�N 3. KRIGING CON DERIVA EXTERNA (Variable explicativa: Elevat) 
# 3.1. KRIGING CON DERIVA EXTERNA
# 3.2. COMPARACI�N ENTRE KRIGING CON DERIVA EXTERNA Y KRIGING ORDINARIO
#________________________________________________________________________________________

# Deriva externa considerando la variable Elevat
CalorV.dekriged = krige(CalorV~ Elevat, coal, coal.grid, model = CalorV.fit)

plot.pred.dekriged <- spplot(CalorV.dekriged, zcol="var1.pred", pretty=T, contour=T, 
                      col.regions=bpy.colors(64), xlab="Etq X", ylab="Etiq Y", scales=list(draw=T),
                      main="Predicciones Krig. Deriva Externa (Elevat)")

plot.var.dekriged <- spplot(CalorV.dekriged, zcol="var1.var", pretty=T, contour=T, 
                      col.regions=bpy.colors(64), xlab="Etq X", ylab="Etiq Y", scales=list(draw=T), 
                      main="Varianzas Krig. Deriva Externa (Elevat)")

print(plot.pred.dekriged, split=c(1,1,1,2), more =T)

print(plot.var.dekriged, split=c(1,2,1,2), more =F)

#_____________________________________________________________________________________
# 3.2. COMPARACI�N ENTRE KRIGING CON DERIVA EXTERNA Y KRIGING ORDINARIO
#_____________________________________________________________________________________

#   C�lculo y representaci�n de diferencias entre kriging ordinario y con Deriva Externa  
#		en predicciones y varianzas				   
plot(CalorV.kriged$var1.pred,CalorV.dekriged$var1.pred,xlab="Predici�n kriging ordinario",
     ylab="Predicci�n kriging Deriva Externa")
abline(0,1,col="red")

dif.dek.ok <- data.frame(dif.pred = CalorV.dekriged$var1.pred - CalorV.kriged$var1.pred, 
                        dif.var  = CalorV.dekriged$var1.var  - CalorV.kriged$var1.var)
summary(dif.dek.ok)

coordinates(dif.dek.ok) <- coordinates(CalorV.kriged)

spplot(dif.dek.ok, zcol="dif.pred", pretty=T, cuts=12, col.regions = cm.colors(64),
       main="Predicciones Krig.Deriva Ext. - Predicciones Krig. Ordin.") 

spplot(dif.dek.ok, zcol="dif.var", pretty=T, col.regions = cm.colors(64), 
       main="Varianzas Krig. Deriva Ext. - Varianzas Krig. Ordin.") 

rm(dif.dek.ok)

#----------------------------------------------------------
# Estudio de los residuos a trav�s de Validaci�n Cruzada
#----------------------------------------------------------

CalorV.dekriged.vc = krige.cv(CalorV~Elevat, coal, model = CalorV.res.fit,nfold=96)
summary(CalorV.dekriged.vc) 

# Comparaci�n con los residuos obtenidos por validaci�n cruzada en el Kriging Ordinario
summary(CalorV.dekriged.vc$residual) # Residuos CV del Kriging con Deriva Externa
summary(CalorV.kvc$residual)        # Residuos CV del Kriging Ordinario
summary(CalorV.dekriged.vc$zscore)   # Zscore CV del Kriging con Deriva Externa
summary(CalorV.kvc$zscore)          # Zscore CV del Kriging Ordinario

###  Histograma de errores estandarizados
par(mfrow = c(2, 1))
hist(CalorV.dekriged.vc$zscore, breaks = seq(-3.5, 3.5, by = 0.5),freq=FALSE,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging con Deriva Ext")
hist(CalorV.kvc$zscore, breaks = seq(-3.5, 3.5, by = 0.5), freq=FALSE ,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging Ordinario")
par(mfrow = c(1, 1))

#________________________________________________________________________________________
#  ACCI�N 4. KRIGING RESIDUAL DIRECTO   (Variable explicativa: Elevat) 
#
# 4.1. Estimaci�n de los par�metros que determinan la deriva a trav�s Del m�todo de 
#      m�nimos cuadrados ordinarios (MCO)
# 4.2. Kriging ordinario sobre los residuos del m�todo de m�nimos cuadrados ordinarios.
# 4.3. C�lculo de la predicci�n del proceso como suma de la deriva ajustada y el ajuste 
#      de  los residuos realizados a trav�s del Kriging
# 4.4. Estudio de los residuos y comparaci�n con el Kriging ordinario
#________________________________________________________________________________________

#_____________________________________________________________________
# 4.1. Estimaci�n de los par�metros que determinan la deriva a trav�s  
#    del m�todo de m�nimos cuadrados ordinarios (MCO)  
#_____________________________________________________________________
#   Se aplica la funci�n "lm" (linear model) indicando                  
#    (var_objetivo~var_explicativas, conjunto_datos)                   

deriva=lm(CalorV~Elevat,coal)
summary(deriva)

# Se salvan los residuos ordinarios en una variable que se incluye en el
# conjunto de datos original (se hace copia par mantener el inicial) 

R0calorv<-residuals(deriva)
coal2=coal
coal2@data=cbind(coal2@data,R0calorv)

#_____________________________________________________________________
# 4.2. Kriging ordinario sobre los residuos del m�todo de m�nimos     
#   cuadrados ordinarios (MCO). 
#_____________________________________________________________________

#   Se inicia el proceso, estimando el modelo de semivariograma te�rico  
#    y se realiza el krigeado, considerando el conjunto de datos          
#    grid o rejilla inicial.  

R0.vgm = variogram(R0calorv~1,coal2)
R0.vgm
plot(R0.vgm)

attributes(fit.variogram(R0.vgm, model=vgm(5, "Sph", 600, 0.5)))$SSErr
attributes(fit.variogram(R0.vgm, model=vgm(5, "Pen", 600, 0.5)))$SSErr
attributes(fit.variogram(R0.vgm, model=vgm(5, "Gau", 600, 0.5)))$SSErr
attributes(fit.variogram(R0.vgm, model=vgm(5, "Cir", 600, 0.5)))$SSErr
attributes(fit.variogram(R0.vgm, model=vgm(5, "Exp", 600, 0.56)))$SSErr

R0.fit = fit.variogram(R0.vgm, model = vgm(5, "Gau", 600, 0.5))
R0.fit
plot(R0.vgm, R0.fit)

R0.kriged = krige(R0calorv~1, coal2, coal.grid, model = R0.fit)

#_____________________________________________________________________
# 4.3. C�lculo de la predicci�n del proceso como suma de la deriva      
#     ajustada y el ajuste de los residuos realizados a trav�s del Kriging 
#_____________________________________________________________________

#     La predicci�n se realiza sobre la rejilla, para lo cual se hace
#     previamente una copia y se incrustan dichas predicciones

Predderiva=predict(deriva,coal.grid)
Predresid=R0.kriged@data$var1.pred

Predfinal=Predderiva+Predresid
coal2.grid=coal.grid
coal2.grid@data=cbind(coal2.grid@data,Predfinal)
spplot(coal2.grid, zcol="Predfinal", pretty=T, contour=T, col.regions=bpy.colors(64), 
       main="Prediciones Kriging Residual", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

#------------------------------
## Estudio de los residuos 
#------------------------------

Predderiva2=predict(deriva,coal)
summary(Predderiva2)
coal2@data=cbind(coal2@data,Predderiva2,errorkrd=coal2$CalorV-Predderiva2)
m=mean(coal2$errorkrd)
s=sd(coal2$errorkrd)
coal2@data=cbind(coal2@data,zserrorkrd=(coal2$errorkrd-m)/s)     
names(coal2)

# Comparaci�n con los residuos obtenidos por validaci�n cruzada en el Kriging Ordinario
summary(coal2$errorkrd) # Residuos CV del Kriging residual directo
summary(CalorV.kvc$residual)        # Residuos CV del Kriging Ordinario
summary(coal2$zserrorkrd)   # Zscore CV del Kriging con Deriva Externa
summary(CalorV.kvc$zscore)          # Zscore CV del Kriging Ordinario

###  Histograma de errores estandarizados
par(mfrow = c(2, 1))
hist(coal2$zserrorkrd, breaks = seq(-3.5, 3.5, by = 0.5),freq=FALSE,
     col = "lightblue",border = "red", main = "Errores estand. Kriging con Deriva Ext (Residual Directo)")
hist(CalorV.kvc$zscore, breaks = seq(-3.5, 3.5, by = 0.5), freq=FALSE ,
     col = "lightblue",border = "red", main = "Errores estandarizados Kriging Ordinario")
par(mfrow = c(1, 1))

#------------------------------