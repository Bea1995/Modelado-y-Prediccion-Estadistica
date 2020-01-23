library(sampling)

data("swissmunicipalities")
datos=swissmunicipalities
names(datos)

N=length(datos$CT)

######################################################################
###    1.SELECCIONAR UNA MUESTRA ALEATORIA SIMPLE DE TAMAÑO  n
######################################################################
###        sample(N,n): tamaño poblacional (N) y tamaño muestral (n)
###        srswor(n,N): simple random sample without replacement
#######################################################################
n=120

### Modo 1
m=sample(N,n) #Vector con los elementos seleccionados
datosmuestrales=datos[m,]

### Modo 2
s=srswor(n,N) #Vector de 0'y 1's
datosmuestrales2=datos[s==1, ]

################################################################
#####  2. SELECCIONAR UNA MUESTRA SISTEMÁTICA (Paso = 25)
################################################################

#Función para realizar un muestreo sistemático uniforme de paso k
#Genera un primer elemento aleatorio y luego coge todos los
#elementos siguientes con paso k
MS=function(N,k){
  m=c() #conjunto
  gamma=1+trunc(k*runif(1)) #elemento aleatorio menor que k
  while (gamma<=N){
    m=c(m,gamma)
    gamma=gamma+k}
  return(m)
}

K=25
m=MS(N,K)

muestra=datos[m, ]
nrow(muestra)

#########################################################################
### 3. SELECCIONAR UN PiPS (Muestreo probabilidades proporcionales a X)
########################################################################
######  Met. Brewer: UPbrewer(pik) (UP: probabilidades desiguales)
######  Met. Midzuno: UPmidzuno(pik)
######  Met. Madow: UPsystematic(pik): paso proporcional a X (variante SIS)
######  Met. Sampford: UPsampford(pik)
######  etc.etc.etc.
########################################################################
n=20

#Probabilidad que le quiero asignar a cada elemento de la población
#POPTOT: población total
#Cada municipio es seleccionado de forma proporcional a su población
pik=n*datos$POPTOT/sum(datos$POPTOT) #Vector de N valores
s=UPmidzuno(pik)
datosmuestrales=getdata(datos,s==1)

### Para cada metodo tenemos la función correspondiente que calcula las 
###  probabilidades de inclusión de segundo orden para toda la población
pik2=UPmidzunopi2(pik) #No se suele utilizar

############################################################################
####  4. MUESTREO ESTRATIFICADO, según variable "REG"
############################################################################
###       strata(data, stratanames=NULL, size, method=c("srswor","srswr",
###             "poisson","systematic"), pik,description=FALSE)
######################################################################

### Hay que ordenar los datos segun la variable de estratificacion
datos=datos[order(datos$REG), ] 

### Tamaño poblacional de cada estrato
N_estrat=as.vector(table(datos$REG))

### Afijación: Tamaño muestral de cada estrato
### Afijación proporcional
n=120
n_estrat=round(N_estrat*n/N)

#Los estratos de la muestra tienen que ser los de la población
#Al hacer un MAS en cada estrato tenemos la misma probabilidad por elemento
#Al tener afijación proporcional todos los estratos tienen la misma
#probabilidad salvo asuntos de redondeo
st=strata(datos, stratanames = "REG",size=n_estrat,method = "srswor")
muestra=getdata(datos,st)

######################################################################
####  5. MUESTREO POR CONGLOMERADOS
######################################################################
###  cluster(data, clustername, size, method=c("srswor","srswr","poisson",
##          "systematic"),pik,description=FALSE)
######################################################################

## MC1-MAS: Seleccion de 3 regiones (conglomerados) segun MAS
cl=cluster(datos,clustername="REG",size=3,method="srswor",description=TRUE)
mc1=getdata(datos,cl)

unique(mc1$REG) #Regiones seleccionadas
table(mc1$REG)  #Tabla de frecuencias para cada una de las regiones

## MC1-PiPS: Seleccion de 3 regiones (conglomerados) segun 
## Probabilidades Proporcionales a la Población 
tapply(datos$POPTOT,datos$REG,sum)

pik=3*as.vector(tapply(datos$POPTOT,datos$REG,sum))/sum(datos$POPTOT)
cl_prop=cluster(datos,clustername=c("REG"),size=3,method="systematic",pik)
table(cl_prop$REG)

#########################################################################
###  6. MUESTREO MULTIETÁPICO: muestreo por conglomerado en varias etapas
#########################################################################
###   mstage(data, stage=c("stratified","cluster",""), varnames, size, 
###          method=c("srswor","srswr","poisson","systematic"), pik, 
###           description=FALSE)
#########################################################################

##  MC2-MAS-MAS  Unid.Prim.= región (REG)  Unid.Secund.=cantón (CT)
m=mstage(datos,stage=list("cluster","cluster"), varnames=list("REG","CT"),
         size=list(4,c(1,1,1,1)), method=list("srswor","srswor"), description = TRUE)

class(m)
class(m[[1]])

table(m[[1]]$REG)
table(m[[2]]$CT)

datosmuestrales=getdata(datos,m)[[2]]
table(datosmuestrales$REG,datosmuestrales$CT)

##  MC2-MAS-MAS  Unid.Prim.= cantón (CT)  Unid.Secund.= municipios (Nom)
m=mstage(datos,stage=list("cluster",""), varnames=list("CT","Nom"),
         size=list(6,rep(2,6)), method=list("srswor","srswor"), description = TRUE)

table(m[[1]]$CT)
table(m[[2]]$Nom)
datosmuestrales=getdata(datos,m)[[2]]

###  Estratificado - MC1 (en cada estrato), selecciono 2 cantones de cada region segun MAS
datos=datos[order(datos$REG), ] 
N_estrat=as.vector(table(datos$REG))
m=mstage(datos, stage=list("stratified","cluster"), varnames=list("REG","CT"), 
         size=list(N_estrat,rep(1,7)),method=list("","srswor")) 
table(m[[1]]$REG)
table(m[[2]]$CT)
xx=getdata(datos,m)[[2]]
table(xx$REG,xx$CT)

###  Estratificado - MC2(MAS,MAS) 
###  (en cada estrato), selecciono 2 cantones de cada region segun MAS
###  y en cada cantón se selecciona un MAS de municipios (tam. muestral n)
###  tamaño muestral n, afijación proporcional al numero de municipios en cada estrato
datos=datos[order(datos$REG), ] 
N_estrat=as.vector(table(datos$REG))
# n=50
n=20
n_estrat=round(N_estrat*n/sum(N_estrat))

m=mstage(datos, stage=list("stratified","cluster",""), varnames=list("REG","CT","Nom"), 
         size=list(N_estrat,rep(1,7),n_estrat),method=list("","srswor","srswor")) 

table(m[[1]]$REG)
table(m[[2]]$CT)

MC=getdata(datos,m)[[2]]
xx=getdata(datos,m)[[3]]
table(xx$REG,xx$CT)
