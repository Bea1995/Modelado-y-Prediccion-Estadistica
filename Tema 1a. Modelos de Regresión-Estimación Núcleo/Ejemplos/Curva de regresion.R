#install.packages("KernSmooth")
library(KernSmooth) #suavizado Kernel

set.seed(1234567)

x=runif(1000) #datos de una ley U(0,1)
x=sort(x)     #ordeno los valores en orden creciente
m.true=sin(6*x)
y=m.true+0.75*rnorm(1000) #Meto ruido a m.true
plot(x,y)
lines(x,m.true, col="blue", lwd=2) #lwd=grosor de la linea
h=dpill(x,y);h #utiliza el método plug in
m.hat.1=locpoly(x,y,bandwidth=h)
#Donde la curva se aproxima bien por una recta, el sesgo es despreciable
#Donde la curvatura es grande, existen discrepancias (hay sesgo)
#En los extremos, el ajuste es malo (no tenemos información de los dos lados)
lines(m.hat.1, col="red", lwd=3)

#Para que ajuste en un rango de las x
m.hat.1.rango=locpoly(x,y,bandwidth=h, range.x=c(0.2,0.6))
plot(x,y)
lines(x,m.true, col="blue", lwd=2)
lines(m.hat.1.rango, col="red", lwd=2)
#En los extremos siguen funcionando mal, aunque sea un trozo de curvatura 0
#Hay que tomar un h distintos en los extremos y en el centro para solucionar el problema

#Estimador de Nadarata-Watson, p=0
m.hat.0=locpoly(x,y,degree=0,bandwidth=h)
#Estimador con p=3
m.hat.3=locpoly(x,y,degree=3,bandwidth=h)
#Dibujamos
plot(x,y)
lines(m.hat.1,col="red",lwd=3)
lines(m.hat.0,col="green",lwd=2) #Es peor en los extremos
lines(m.hat.3,col="blue",lwd=2) #Sale abollado

#Podemos estimar las derivadas
m.prima.true=6*cos(6*x)
m.prima.hat=locpoly(x,y,drv=1,bandwidth=h)
#Dibujamos
plot(x,m.prima.true,col="green",type="l", lwd=2)
lines(m.prima.hat, col="red", lwd=2)
#Sale muy rugoso porque hemos usado la h de m
m.prima.hat=locpoly(x,y,drv=1,bandwidth=2*h)
plot(x,m.prima.true,col="green",type="l", lwd=2)
lines(m.prima.hat, col="red", lwd=2)
#En las lineas rectas coincide, en los extremos lo hace mal
#Tocará el valle con muchos más elementos de muestra

library(np)
h
h2=npregbw(y~x);h2

