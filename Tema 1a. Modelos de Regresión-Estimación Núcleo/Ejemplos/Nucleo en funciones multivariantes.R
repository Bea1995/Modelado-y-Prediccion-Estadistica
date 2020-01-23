#Vamos a calcular el nucleo para funciones de dimension d>=2
install.packages("ks")
library("ks")

#Matriz con datos en la filas y dimensión en las columnas
#Calculamos datos para una distribución normal y lo dividimos en dos columnas
#para tener una muestra bivariante: tenedremos X_1,...,X_500 iid N_2(0,I_2)
x=matrix(rnorm(1000),ncol=2)

#Regla de Silverman (exacto)
(Hn=Hns(x)) #Ejecuta y pinta: sentencia entre paréntesis

#Cross-validation (la version diagonal tarda menos)
(Hv=Hlscv(x)); (Hv=Hlscv.diag(x))

#Plug-in en 2 pasos
(Hp=Hpi(x)); (Hp=Hpi.diag(x))

#Nucleo para nuestra distribucion bivariante
fthat=kde(x)
plot(fthat) #Curva de nivel (el exacto son círculos cocentricos)

#La primera curva de nivel deja por debajo el 25% de la muestra, la segunda el 50% y la tercera el 75%

#Para pintar 10 curvas de nivel (cada 10%)
v=1:9*10
plot(fthat,cont=v)

#Tengo una rejilla fina en X e Y y estima la densidad en las Z
x1=fthat$eval.points[[1]]
x2=fthat$eval.points[[2]]
x3=fthat$estimate
#Dibujamos la perspectiva
persp(x1,x2,x3)
persp(x1,x2,x3,theta=45) #Gira la perspectiva 45 grados
persp(x1,x2,x3,theta=45,phi=45)

#Ejemplo con datos reales: faithful
data=faithful
head(data)
plot(data)

#Funcion de densidad de los datos
fthat=kde(data)
plot(fthat)

#Calculamos la perspectiva de los datos
x1=fthat$eval.points[[1]]
x2=fthat$eval.points[[2]]
x3=fthat$estimate

persp(x1,x2,x3)
