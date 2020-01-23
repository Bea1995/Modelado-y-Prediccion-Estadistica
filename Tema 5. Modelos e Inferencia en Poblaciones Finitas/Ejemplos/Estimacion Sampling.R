########################################################################################################################
###
###           ESTIMACIÓN con la LIBRERIA SAMPLING
###
###  Supongamos dado un conjunto de datos muestrales  
###        m = elementos seleccionados de la población
###        Ym = variable de interés en los elementos muestrales 
###        pikm = probabilidades de inclusión (Pi_i) en los elementos muestrales
###        pik2m = matriz n x n, probabilidades de inclusión de segundo orden (Pi_ij) en los elementos muestrales
########################################################################################################################
library(sampling)

###### ESTIMADOR DE HORVITZ-THOMPSON (del Total, T(Y) )
HTestimator(Ym,pikm)  ## Estimador
varHT(Ym,pik2m)       ## Estimación insesgada de la varianza

Zm=n*Ym/pikm
(1/n)*var(Zm)    ### Estimación de la varianza, sup. con reemplazamiento

####### ESTIMADOR DE HAJEK (de la media, media(Y) )
est.N=sum(1/pikm)
HTestimator(Ym,pikm)/est.N

Hajekestimator(Ym,pikm,type="mean")

est.HJ=Hajekestimator(Ym,pikm,type="mean")
Zm=(Ym-est.HJ)/(est.N)  ## Para utilizar aproximación lineal

### Estimacion insesgada de la aprox. lineal de la varianza
varHT(Zm,pik2m)   

Zm2=n*Zm/pikm

### Estimación de la varianza (aprox. lineal), sup. con reemplazamiento
(1/n)*var(Zm2)   

### La función vartaylor_ratio(Ym,Xm,pik2m) 
### obtiene el estimador(R)=estimador(T(Y))/estimador(T(X))
### y el estimador insesgado de la aprox. lineal de la varianza

## Haciendo X=1 --> est(R)=est(T(Y))/est(N)=estimador de Hajek de media(Y)
vartaylor_ratio(Ym,rep(1,n),pik2m) 
