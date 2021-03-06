###################################################################################
###################################################################################
####									                                                          ###
####	 	              MASTER UNIVERSITARIO EN MATEM�TICAS		                     ###        	
####									                                                          ###   
###             Asignatura: MODELADO Y PREDICCI�N ESTAD�STICA                   ###
###                                                                             ###
###                 TEMA 3. MODELOS PARA DATOS ESPACIALES                       ###
###                                                                             ###
###                                                                             ###
###                                                                             ###
###       A.M. MU�OZ-REYES - J.M. MU�OZ-PICHARDO (Universidad de Sevilla)       ###
###################################################################################

#########                SCRIPT 1 : KRIGING ORDINARIO                     #########
#########                  Conjunto de datos: Meuse                       ######### 
###################################################################################

#################### CARGA DEL PAQUETE GSTAT
library(sp)
library(lattice)
library(xts)
library(gstat)


#################### LECTURA DE DATOS
data(meuse)
class(meuse)
names(meuse)

## NOTA SOBRE EL FICHERO DE DATOS #####################################
## Las observaciones est�n georreferenciadas en coordenadas UTM (x e y)
## En cada localizaci�n se han recogido:
##   - las concentraciones de 4 metales pesados  (cadmio, cobre, plomo y zinc)
##   - variables geof�sicas: 
##       -- elevaci�n sobre el nivel de lecho del rio ("elev"; en metros)
##       -- distancia topom�trica al grid o rejilla(?) ("dist") 
##       -- contenido de materia org�nica (%) del suelo ("om") , 
##       -- freq.inundaci�n ("ffreq"): 1=1 vez en 2 a�os; 2=1 vez en 10 a.; 3=1 vez en 50 a.
##       -- tipo de suelo ("soil"): 1-caliza; 2-arcilla pesada; 3-arcilla limosa
##       -- clase de limo ("lime"): 0-ausente; 1-presente
##       -- uso de la parcela ("landuse")- diversas modalidades
##       -- distancia al rio en m. ("dist.m"), obtenida durante el trabajo de campo

###### Asignaci�n de coordenadas referencias geogr�ficas
coordinates(meuse) = ~x+y
class(meuse) #DataFrame de puntos espaciales
names(meuse)
plot(meuse)

####################################################################################
###
###	1. AN�LISIS ESTRUCTURAL
###
####################################################################################

#################### RESUMEN Y REPRESENTACION GR�FICA (BURBUJAS)
summary(meuse)

hist(meuse$zinc, breaks = seq(0, 1900, by = 100), col = "lightblue",border = "red", main = "Concentraci�n de zinc (peso: ppm)")
rug(meuse$zinc)

bubble(meuse, c("zinc"),col=c("#00ff0088", "#00ff0088"), main = "concentraci�n de Zinc (ppm)")


## Relaci�n lineal entre log(zinc) y raiz cuadrada(dist)
plot(log(zinc)~sqrt(dist), meuse)
abline(lm(log(zinc)~sqrt(dist), meuse))

summary(log(meuse$zinc))

hist(log(meuse$zinc), breaks = seq(4, 8, by = 0.3), col = "lightblue",border = "red", main = "Log-Concentraci�n de Zinc (peso: ppm)")
rug(log(meuse$zinc))

#################### PREPARACI�N DEL "GRID" REJILLA   
data(meuse.grid)

## NOTA SOBRE EL FICHERO DE DATOS mesure.grid ################################################
## Las observaciones est�n georreferenciadas en coordenadas UTM (x e y)
## En cada localizaci�n se han recogido:
##   - x     : Coordenada X 
##   - y     : Coordenada Y 
##   - dist  : distancia al borde de r�o Meuse, espacial), normalizada a [0,1]
##   - ffreq : freq.inundaci�n ("ffreq"): 1=1 vez en 2 a�os; 2=1 vez en 10 a.; 3=1 vez en 50 a.
##   - soil  : tipo de suelo ("soil"): 1-caliza; 2-arcilla pesada; 3-arcilla limosa
##   - part.a : divisi�n arbitraria del �rea en dos zonas, Zona A
##   - part.b : divisi�n arbitraria del �rea en dos zonas, Zona B
#############################################################################################

coordinates(meuse.grid) = ~x+y

gridded(meuse.grid) = TRUE
class(meuse.grid) #Base de datos de p�xeles
names(meuse.grid)

spplot(meuse.grid, c("dist"), col.regions=terrain.colors(20),main="DISTANCIA AL RIO MEUSE", scales=list(draw=TRUE),xlab="X",ylab="Y")

spplot(meuse.grid, c("part.a","part.b"), col.regions=terrain.colors(20),main="Partes A y B", scales=list(draw=TRUE),xlab="X",ylab="Y")

####################################################################################
###
###	2. CONSTRUCCI�N DEL VARIOGRAMA MUESTRAL 
###
####################################################################################

### 	Funci�n "variogram": variogram(object, data, locations, cutoff, Cressie,...)
###
### Calcula el variograma muestral, con opciones para variograma direccional y robusto.    
###    
###    object : objeto de la clase "gstat". Se indica "variable ~ regresor1+regresor2+..."         
###             En caso de ausencia de regresores, utilizar "variable ~ 1"
###    data   : Conjunto de datos
###    cloud  : (l�gico), si es TRUE, c�lculo del semivariograma nube
###    locations : ubicaciones de datos espaciales
###    cutoff : valor de distancia de separaci�n m�xima para incluir los 
###             pares  de puntos en las estimaciones  de la semivarianza;
###             Por defecto, 1/3 de la longitud de la diagonal del marco que envuelve los datos 
###   cressie = (l�gico) TRUE (Estimaci�n robusta de Cressie del variograma); 
###                      FALSE (m�todo cl�sico de los momentos) [---> Valor por defecto]

##########################################################
###        SEMIVARIOGAMA NUBE
##########################################################
#Se usa para detectar outliers
#Kriging sin ninguna variable auxiliar (indicado con el 1)
lzn.cloud <- variogram(log(zinc) ~ 1, meuse, cloud = T)
head(lzn.cloud)
plot(lzn.cloud,col="blue",main="Semivariograma nube Log(Zinc)")

##########################################################
###       SEMIVARIOGAMA EXPERIMENTAL (MET. MOMENTOS)
##########################################################
lzn.vgm = variogram(log(zinc)~1, meuse)
lzn.vgm
plot(lzn.vgm, col="black",main="Semivariograma experim. Log(Zinc)")

##########################################################
###       ESTIMACI�N ROBUSTA DE CRESSIE
##########################################################
#M�s robusta frente a determinadas variaciones de las hipotesis que se plantean
lzn.Cressie = variogram(log(zinc)~1, meuse, cressie=TRUE)
lzn.Cressie
plot(lzn.Cressie, col="black",main="Semivariograma de cressie. Log(Zinc)")

####################################################################################
###
###	3. AJUSTE DE MODELOS TEORICOS DE VARIOGRAMA
###
####################################################################################

###                   Funci�n "vgm": vgm(s, "Sph", r, n)
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
######################################################################## 

vgm() #Genera una funci�n de variograma

lzn.fit = fit.variogram(lzn.vgm, model = vgm(0.5, "Sph", 900, 0.1))
lzn2.fit = fit.variogram(lzn.vgm, model = vgm(0.5, "Sph", 900, 0.1),fit.method=2)
plot(lzn.vgm, lzn.fit,main="ajuste variograma")
plot(lzn.vgm, lzn2.fit,main="ajuste variograma (Cressie)")
lzn.fit
lzn2.fit

############### BONDAD DE AJUSTE DEL MODELO TEOR�CO DE VARIOGRAMA 
### Suma de cuadrados de los residuos, ponderados seg�n el modelo ajustado
attributes(lzn.fit)$SSErr

### Si se desea realizar diversos ajustes para elegir el m�s adecuado
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1)))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Pen", 900, 1)))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Gau", 900, 1)))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Cir", 900, 1)))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Exp", 900, 1)))$SSErr

attributes(fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1),fit.method=2))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Pen", 900, 1),fit.method=2))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Gau", 900, 1),fit.method=2))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Cir", 900, 1),fit.method=2))$SSErr
attributes(fit.variogram(lzn.vgm, model=vgm(1, "Exp", 900, 1),fit.method=2))$SSErr

print(plot(lzn.vgm, plot.numbers = F, pch =18, col = "darkblue", model = lzn.fit))

####################################################################################
###
###	4. KRIGING ORDINARIO
###
####################################################################################
###
###                   Funci�n "krige": krige(object, data, newdata, model,...)
###
###     Realiza las predicciones kriging sobre los puntos definidos en el grid (o rejilla)
#####       Genera un un marco de datos con dos variables: 
#####           "var1.pred" = Predicciones
#####		"var1.var" = Varianzas de las predicciones
#####
#####     object : objeto de la clase "gstat". Se indica "variable ~ regresor1+regresor2+..."         
###             En caso de ausencia de regresores, utilizar "variable ~ 1"
###       data   : Conjunto de datos
###       newdata  :  Grid donde realizar las predicciones
###       model  :  modelo de variograma obtenido por la ejecuci�n de "vgm" o "fit.variogram"
############################################################################################
lzn.kriged = krige(log(zinc)~1, meuse, meuse.grid, model = lzn.fit)

###################  Para que aparezcan en pantalla los valores de la variable var1.pred
lzn.kriged$var1.pred

###################  Para que aparezcan en pantalla los valores de todas las variables
lzn.kriged@data

###################   Construcci�n de tabla de resultados (a�adiendo coordenadas de los puntos) 
tabla<-cbind(coordinates(lzn.kriged),lzn.kriged$var1.pred,lzn.kriged$var1.var)
colnames(tabla)<-c("Coord X","Coord Y","Predicci�n","Varianza Pred.")
tabla

############## 	REPRESENTACIONES GR�FICAS

### Plot espacial de la predicci�n con graduaci�n de colores 
spplot(lzn.kriged["var1.pred"],main="Plot espacial de la predicci�n" )

### Plot espacial de la predicci�n y la varianza de la predicci�n con graduaci�n de colores
spplot(lzn.kriged, main="Plots espaciales de la predicci�n y la varianza de la predicci�n" )

#####################################################################################################
### NOTA:   En los spplots se puede cambiar la gama de colores a trav�s del argumento: col.regions=....
#####################################################################################################

### Plot espacial de la predicci�n con otra graduaci�n de colores
spplot(lzn.kriged, zcol="var1.pred", pretty=T, col.regions=bpy.colors(64), main="Titulo deseado", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

### .... y con l�neas de contorno 
spplot(lzn.kriged, zcol="var1.pred", pretty=T, contour=T, col.regions=bpy.colors(64), main="Titulo deseado", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

### Plot espacial de la varianza de la estimaci�n (con otra graduaci�n de colores y l�neas de contorno)
spplot(lzn.kriged, zcol="var1.var", pretty=T, contour=T, col.regions=terrain.colors(64), main="titulo deseado 3", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

### Plot espacial de la varianza de la estimaci�n (con otra graduaci�n de colores sin l�neas de contorno)
spplot(lzn.kriged, zcol="var1.var", pretty=T, contour=F, col.regions=bpy.colors(10), main="titulo deseado", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T))

###############################################################################################################################
##########
##########		Seleccionamos los 2 primeros puntos del meuse.grid, y 6 puntos del meuse
##########			para obtener el estimador kriging sobre estos puntos (y analizar la salida con debug.level=32)
##########			Enviar la salida de pantalla a un fichero con la orden sink
################################################################################################################################
meuse.selec <- meuse[1:6,]
selec.grid <- meuse.grid[1:2, ]

sink("Resultados.txt")

lzn.krigsel = krige(log(zinc)~1, meuse.selec, selec.grid, model = lzn.fit, debug.level=32)

##############      Notas sobre DEBUG.LEVEL   ###########################################################
####
### Valores para el argumento debug.level: 
####
### 	0: Suprime todas las salidas, salvo mensajes de aviso y error; 
###	1: salida normal (valor por defecto): informe breve sobre datos, tipo de kriging usado, progreso del programa en %, tiempo de ejecuci�n total;
### 	2: imprime modelo de variograma utilizado, ficheros leidos y escritos, incluye nombre del fuchero fuente y n�mero de linea en mensaje de error;
### 	4: imprime  diagnosticos de ajuste MCO (M�nimos Cuadsrados Ordinario) and MCP (M�n. Cuad. Ponderado); 
### 	8:  imprime todos los datos despu�s de leerlos; 
### 	16: imprime el entorno para cada predicci�n local; 
### 	32: imprime matrices de covarianzas, matrices del dise�o, soluciones, pesos kriging, etc.;
### 	64: imprime diagn�sticos del ajuste al  variograma (n�mero de iteraciones y modelo de variogra<ma en cada paso de iteraci�n);
### Para combinar salidas, sumar los respectivos valores. 
### Cualquiera de los valores anteriores en negativo, a�ade a la salida un contador del proceso ejecutado
######################################################################################################################################
sink()

######################################################################################################################################
#######################             C�LCULO DE LA PREDICCI�N DEL ZINC  (deshaciendo la transformaci�n del logaritmo neperiano)
zn.kriged.pred <- exp(lzn.kriged$var1.pred)

#################   Para guardar en un nuevo objeto (zn.kriged) las coordenadas junto con las predicciones del zinc: 
#################   Se crea el objeto zn.kriged (copiando lzn.kriged) y se le a�ade la variable calculada antes   "zn.kriged.pred"
#################    
zn.kriged <- lzn.kriged
zn.kriged@data = cbind(zn.kriged@data,zn.kriged.pred)

names(zn.kriged)

spplot(zn.kriged["zn.kriged.pred"],main="Plot espacial de la predicci�n del Zinc",col.regions=bpy.colors(10) )

##########################################################################################################################################
###########	 	KRIGING ORDINARIO TOMANDO LAS OBSERVACIONES M�S CERCANAS 			################################## 
##########################################################################################################################################

######################		FIJANDO EL N�MERO M�XIMO DE OBSERVACIONES			########################################
lzn.kriged.nmax = krige(log(zinc)~1, meuse, meuse.grid, model = lzn.fit, nmax = 30)

spplot(lzn.kriged.nmax, zcol="var1.pred", contour=T, col.regions=terrain.colors(64), main="Predicciones KRIGING (ORDINARIO) CON 30 OBSERVACIONES")

######################		FIJANDO LA DISTANCIA M�XIMA					########################################
lzn.kriged.maxdist = krige(log(zinc)~1, meuse, meuse.grid, model = lzn.fit, maxdist = 900, nmin = 12)

spplot(lzn.kriged.maxdist , zcol="var1.pred", pretty=T, contour=T, col.regions=terrain.colors(64),main = "Predicciones KRIGING (ORDINARIO) a un radio de 900m.")

#################		Comparaci�n gr�fica con el kriging ordinario  (predicciones y varianzas)
#################			(Para hacerlo a la misma escala se utiliza el argumento "at")
range(lzn.kriged.maxdist$var1.pred, lzn.kriged$var1.pred, na.rm=TRUE)

range(lzn.kriged.maxdist$var1.var, lzn.kriged$var1.var,  na.rm=TRUE)

at.pred <- 4:8
at.var <- seq(0.05, 0.6, by=0.1)

plot.1 <-spplot(lzn.kriged.maxdist , zcol="var1.pred", pretty=T, contour=F, col.regions=bpy.colors(64),main = "Predicciones KRIGING (ORDINARIO) - a un radio de 900m.", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T), at=at.pred)

plot.2 <- spplot(lzn.kriged, zcol="var1.pred", pretty=T, contour=F, col.regions=bpy.colors(64), main="Predicciones KRIGING (ORDINARIO)", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T), at=at.pred)

print(plot.1, split=c(1,1,2,1), more =T)

print(plot.2, split=c(2,1,2,1), more =F)

plot.1 <-spplot(lzn.kriged.maxdist , zcol="var1.var", pretty=T, contour=F, col.regions=bpy.colors(64),main = "Varianzas KRIGING (ORDINARIO) - a un radio de 900m.", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T), at=at.var)

plot.2 <- spplot(lzn.kriged, zcol="var1.var", pretty=T, contour=F, col.regions=bpy.colors(64), main="Varianzas KRIGING (ORDINARIO)", xlab="Etq X", ylab="Etiq Y", scales=list(draw=T), at=at.var)

print(plot.1, split=c(1,1,2,1), more =T)

print(plot.2, split=c(2,1,2,1), more =F)

####################################################################################
###
###	5. M�TODO DE VALIDACI�N CRUZADA 
###
####################################################################################
####
#### Funci�n krige.cv
####
####		El argumento principal es "nfold=k", con k entero mayor que 1
#### 		El procedimiento validaci�n cruzada de k iteraciones o K-fold cross-validation 
#### 		los datos de muestra se dividen en k subconjuntos. Uno de los subconjuntos se utiliza 
#### 		como datos de prueba y el resto (k-1) como datos de entrenamiento. El proceso se repite 
#### 		k iteraciones, con cada uno de los posibles subconjuntos de datos de prueba. 
#### 		Finalmente se realiza la media aritm�tica de los resultados de cada iteraci�n para obtener 
#### 		un �nico resultado.  
#### 		Para k=n�mero de observaciones, aplica el m�todo de validaci�n cruzada dejando uno fuera
##########################################################################################################
lzn.kvc <- krige.cv(log(zinc)~1, meuse, model = lzn.fit, nfold=155)

### 	Esta orden crea un objeto con las siguientes variables
###  		var1.pred (predicciones); var1.var (varianza de la predicci�n)
###  		observed (valores observados) ; residual (residuos)
###  		zscore (residuos estandarizados; fold (n�mero de grupo o submuestra)

############################################################################################################
#####   Criterios para evaluar el ajuste del variograma a traves de la validaci�n cruzada	############

############ Medidas:  Medias de los residuos
summary(lzn.kvc)
mean(lzn.kvc$residual)
mean(lzn.kvc$zscore^2)

### Porcentaje de datos con un error absoluto superior a epsilon 
epsilon <- 1
errorabs.sup <- abs(lzn.kvc$residual)> epsilon
table(errorabs.sup)/length(errorabs.sup)

############ Representaciones gr�ficas

####  Mapa de la ubicaci�n de los puntos con marcas en funci�n del tama�o de los errores estandarizados
bubble(lzn.kvc, "residual", main = "Cross Validaci�n: Residuos")
bubble(lzn.kvc, "zscore", main = "Cross Validaci�n: Residuos Estandarizados")

###  Diagrama de dispersi�n de los valores estimados frente a los observados
plot(lzn.kvc$observed,lzn.kvc$var1.pred)

### Diagrama de dispersi�n de los valores estimados frente a los residuos estandarizados
plot(lzn.kvc$var1.pred, lzn.kvc$zscore)

###  Histograma de errores estandarizados
hist(lzn.kvc$zscore, breaks = seq(-2.5, 3.5, by = 0.25), col = "lightblue",border = "red", main = "Errores estandarizados")
