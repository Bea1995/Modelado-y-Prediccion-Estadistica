#Teorema de Glivenko-Cantelli: sup_x(|F_n(x)-F(x)|)->0 (c.s.)

#Fijamos la semilla
set.seed(1234567)

#Generamos datos de una normal
n=100

# rnorm:datos de una normal (pseudorandom data generation)
# dnorm:densidad de una normal (fdd ó fdProb)
# pnorm:funcion de distribucion de una normal (fdD)
# qnorm: cuantiles (Quantilefuncion, F^-1)

#Las letras r,d,p,q+distribución estan reservadas para estos conceptos
x=rnorm(n)
x
Fn=ecdf(x) #calcula la funcion de distribucion empirica asociada
plot(Fn)
#Añadimos al gráfico (add=T) la curva de la normal
curve(pnorm,add=T,col="red")

#Para 1000 datos
n=1000
x=rnorm(n)
Fn=ecdf(x)
plot(Fn)
curve(pnorm,add=T,col="red")

#Para 10000 datos
n=10000
x=rnorm(n)
Fn=ecdf(x)
plot(Fn)
curve(pnorm,add=T,col="red")

#Ponemos los tres graficos en una misma pantalla
layout(matrix(1:3, ncol=3))
enes=c(100,1000,10000)
for(i in 1:3)
{n=enes[i]
 x=rnorm(n)
 Fn=ecdf(x)
 plot(Fn, main=n)
 curve(pnorm,add=T,col="red")
}

#Volvemos a dejar una sola pantalla
layout(matrix(1, ncol=1))
