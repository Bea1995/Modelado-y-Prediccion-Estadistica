########################################################################
###
###           ESTIMACI�N con la LIBRERIA SURVEY
###
#########################################################################
###  Hay que definir el dise�o muestral aplicado a traves de la funci�n 
###  svydesign  (procedimiento usado para elegir los elementos)
#########################################################################
library(survey)

dis.MAS=svydesign(ids=m, fpc=rep(N,n))  ### m=sample(N,n) 

#######################################
####  ESTIMACION EN DISE�OS COMPLEJOS
#######################################

## MC2-MAS-MAS  Unid.Prim.= cant�n (CT)  Unid.Secund.= municipios (Nom)

#### Selecci�n de la muestra a trav�s de la funci�n mstage
library(sampling)

data("swissmunicipalities")
datos=swissmunicipalities
N=nrow(datos)

m=mstage(datos,stage=list("cluster",""), varnames=list("CT","Nom"),
         size=list(6,rep(2,6)), method=list("srswor","srswor"))

datosmuestra=getdata(datos,m)[[2]]

Ym=datosmuestra$Surfacesbois

dis.MC2= svydesign(ids=~CT+1, probs=datosmuestra$Prob , data=datosmuestra )
confint(svymean(Ym,design=dis.MC2))

## Estimacion de total (HT) 
## Estima varianza sup. con reemplazamiento
t=svytotal(Ym,design=dis.MC2)
SE(t)
coef(t)

## Estimaci�n de la media (Hajek)  
## Estima varianza sup. con reemplazamiento
svymean(Ym,design=dis.MC2)   
