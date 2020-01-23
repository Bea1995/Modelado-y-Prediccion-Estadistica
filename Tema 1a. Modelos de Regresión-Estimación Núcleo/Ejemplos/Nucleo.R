#Calculamos el estimador nucleo
d=density(x)
plot(d)
curve(dnorm,add=T,col="red")
hist(x,breaks=c(-8:8)/2,freq=F,add=T,border="green")

#Vamos a probar con varios nucleos
layout(matrix(1:6,ncol=3))

#Gaussiano
d=density(x,kernel="gaussian")
plot(d,main="nucleo gaussiano")

#Epanechnikov
d=density(x,kernel="epanechnikov")
plot(d,main="nucleo epanechnikov")

#Rectangular
d=density(x,kernel="rectangular")
plot(d,main="nucleo rectangular")

#Triangular
d=density(x,kernel="triangular")
plot(d,main="nucleo triangular")

#Biweight
d=density(x,kernel="biweight")
plot(d,main="nucleo biweight")

#Cosine
d=density(x,kernel="cosine")
plot(d,main="nucleo cosine")

layout(matrix(1:6,ncol=3))

#Podemos observar que casi no se aprencian diferencias entre núcleos

#Vamos a cambiar el parámetro ventana

#Silverman (exacta)
bw.nrd(x)
#Cross-validation
bw.ucv(x)
#Plug-in en 4 pasos
bw.SJ(x)

#Cross-validation y plug-in convergen muy lentamente

#Veamos el efecto del parámetro ventana en nuestros datos
layout(matrix(1:4,ncol=2))

h=0.265 #Silverman
d=density(x,bw=h)
plot(d,main="h=0.265")

h=0.14 #Cross-Validation
d=density(x,bw=h)
plot(d,main="h=0.14")

h=0.07
d=density(x,bw=h)
plot(d,main="h=0.07")

h=0.035
d=density(x,bw=h)
plot(d,main="h=0.035")

#Con un h pequeño, sale muy rugosa la densidad

#Un h muy grande oculta caracteristicas de la distribucion
#Un h muy pequeño conlleva a una densidad muy arrugada e irregular

layout(matrix(1,ncol=1))

#Silverman
h1=bw.nrd(x); plot(density(x,bw=h1))
#Cross-validation
h2=bw.ucv(x); lines(density(x,bw=h2),col="blue")
#Plug-in en 4 pasos
h3=bw.SJ(x); lines(density(x,bw=h3),col="red")

