#Histograma de los datos
set.seed(1234567)
n=1000
x=rnorm(n)

#hist(x) nos daría un grafico con la frecuencia observada
#Para la frecuencia relativa:
hist(x,freq=F)

#Para cambiar la anchura (poner los intervalos que queramos):
hist(x,breaks=c(-8:8)/2,freq=F)
hist(x,breaks=c(-12:12)/3,freq=F) #24 intervalos
#Nos sale muy rugoso

hist(x,breaks=c(-8:8)/2,freq=F)
curve(dnorm,add=T,col="red")

hist(x,breaks=c(-12:12)/3,freq=F)
curve(dnorm,add=T,col="red")

hist(x,breaks=c(-16:16)/4,freq=F)
curve(dnorm,add=T,col="red")

#Vamos a mejorarlo con el estimador núcleo