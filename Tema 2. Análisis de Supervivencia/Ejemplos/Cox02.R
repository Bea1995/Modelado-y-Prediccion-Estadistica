#__________________________________________________________________________________
#                     Modelado y Predicci�n Estad�stica.                       
#                M�ster en Matem�ticas. Universidad de Sevilla                
#                           Juan Manuel Mu�oz Pichardo                        
#                                                                             
#__________________________________________________________________________________
#
#   AN�LISIS DE SUPERVIVENCIA     
#                SCRIPT   cox02                                  
#                Estimaci�n de Kaplan Meier                       
#                Regresi�n de Cox                                
#                Datos: uissurv.csv                              
#__________________________________________________________________________________


# CARGA DEL PAQUETE survival
library(survival)

# Algunas opciones de los c�lculos
options(digits = 4) # n�mero de decimales en las salidas
options(columns = 40) # numero de columnas en las salidas

## DESCRIPCI�N DE LOS DATOS ##
# Fuente:
#   University of Massachusetts AIDS Research Unit (UMARU) IMPACT Study (UIS).
#   Provided by Drs. Jane McCusker, Carol Bigelow and Anne Stoddard.
#   https://www.umass.edu/statdata/statdata/stat-survival.html
#
# REFERENCIA:
#   Hosmer, D.W. and Lemeshow, S. and May, S. (2008) 
#   Applied Survival Analysis: Regression Modeling of Time to Event Data. 2nd Edition, 
#   John Wiley and Sons Inc., New York.
#
#  El estudio UIS fue un proyecto de investigaci�n colaborativo, de 5 a�os (1989-1994),
#  (Benjamin F. Lewis, PI, del National Institute on Drug Abuse Grant #R18-DA06151) 
#  El objetivo del estudio fue comparar los programas de tratamiento de diferentes duraciones 
#  previstas dise�ados para reducir el abuso de drogas y prevenir conductas de alto riesgo de VIH.
# 
#  El UIS busc� determinar si los enfoques de tratamiento residenciales alternativos son variables 
#  en efectividad y si la eficacia depende de la duraci�n del programa dise�ado.
#  
# LISTADO DE VARIABLES
#   Name		Description				Codes/Values
# ********************************************************************************************************************
#  1 id		   Identification                       Code 1 - 628
#  2 age		   Age at Enrollment                    Years
#  3 beck      Beck Depression Score	              0.000 - 54.000
#  4 hercoc		 Heroin/Cocaine Use During	          1 = Heroin & Cocaine
#              3 Months Prior to Admission	        2 = Heroin Only 
#                                                   3 = Cocaine Only
#                                                   4 = Neither Heroin 
#                                                   5 = Neither 
#  5 ivhx		   IV Drug Use History at		            1 = Never
#              Admission		                        2 = Previous
#                                                   3 = Recent
#  6 ndrugtx	 Number of Prior Drug Treatments      0 - 40
#  7 race		   Subject's Race			                  0 = White
# 						                                      1 = Non-White
#  8 treat	  Treatment Randomization		            0 = Short
#   		      Assignment		                        1 = Long
# 9 site		  Treatment Site			                  0 = A , 1 = B
# 10 los		  Length of Stay in Treatment	          Days
# 		        (Admission Date to Exit Date)
# 11 time		  Time to Drug Relapse		              Days
# 		        (Measured from Admission Date)
# 12 censor		Event for Treating Lost to	          1 = Returned to Drugs or Lost to Follow-Up
# 		        Follow-Up as Returned to Drugs 	      0 = Otherwise
# 
# 
#__________________________________________________________________________________

##      ENTRADA DE DATOS
datos = read.csv("uissurv.csv",header=TRUE, sep=";")  # read csv file 

# Otra forma: datos <- read.table("uissurv.csv", header=TRUE, sep=";")

dim(datos)
names(datos)
summary(datos)

##-----------------------------------------------------------------------------------
##
## (a) Obtener la estimaci�n de Kaplan-Meier de la Funci�n de Supervivencia para el tiempo 
##     hasta la reacaida con su correspondiente error de estimaci�n
##     y un intervalo de confianza al 95%.
##     
##
##-----------------------------------------------------------------------------------


##
## Estimaci�n de Kaplan-Meier de la Funci�n de Supervivencia
## 
## Si se quiere que el IC se base solo en el ES, se calcule seg�n la f�rmula de Greenwood:
##    km.fit <- survfit(Surv(time,status)~1,type="kaplan-meier",data = aml1,conf.type="plain")
## Por defecto el IC se calcula mediante el m�todo "log" (ver p�g.30 de (1)), basado 
##    en el m�todo delta:
##       exp[ log(FS(t))+- 1.96 SE(H(t)) ]
##    donde H(t) es la estimaci�n de la funci�n tasa de fallo acumulado y 
##    SE(H(t)) su error de estimaci�n 
##    Tanto el proporcionado por la f�rmula de Greenwood como el m�todo "log" pueden proporcionar
##    intervalos fuera del intervalo  [0,1], lo cual es "il�gico" dado que son intervalos de confianza 
##    de la probabilidad de supervivencia. Para evitar este problema, se puede usar el m�todo log-log
##    propuesto en Kalbfleisch & Prentice (1980) : W = log(-log( FS(t)))
##          conf.type="log-log" 
##
## (1) Mara Tableman. "Survival Analysis Using S/R.". 
##     Notes of the first six chapters of the book 
##     Mara Tableman and Jong Sung Kimz. "Survival Analysis Using S: Analysis of
##     Time-to-Event Data". Chapman & Hall/CRC, Boca Raton, 2004
##
##-----------------------------------------------------------------------------------

km.fit <- survfit(Surv(time,censor)~1,type="kaplan-meier",data = datos)
km.fit
summary(km.fit)
attributes(km.fit)
summary(km.fit)$time #los tiempos no censurados
summary(km.fit)$n.censor

## Representaci�n gr�fica de FS estimada por K-M 

plot(km.fit,conf.int=FALSE,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Proporci�n sin reincidencia",lab=c(10, 10, 7),col="blue",lwd=2,
     main="Estimaci�n de S(t) datos uissurv.csv")
abline(h=0) 
grid()

## Tabla de estimaci�n de la FS e intervalos de confianza 

superv<-summary(km.fit)$surv  #F.Superv. estimada por K-M
inferior<- summary(km.fit)$lower 
superior<- summary(km.fit)$upper
cbind(inferior,superv,superior)

## Representaciones gr�ficas de FS y H estimadas por K-M con Intervalos de confianza
##     fun="cumhaz"  representa funci�n H estimada ( H(t)=-log(S(t))
##     fun="log"     representa curva log-supervivencia (eje marcado con log (S) 
##     fun="sqrt"    representa curva raiz-supervivencia
##     fun="event"   representa n�mero de eventos acumulados 1-S(t) 
##     fun="cloglog" representa log[H(t)]=log[-log(S(t)]

plot(km.fit,conf.int=T,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Proporci�n sin recaida",lab=c(10, 10, 7),col="blue",lwd=2,
     main="Estimaci�n K-M de S(t). Datos datos.txt")
abline(h=0) 
grid()

plot(km.fit,conf.int=T,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),col="blue",lwd=2,fun="cumhaz",
     main="Estimaci�n K-M de la tasa de fallo acumulada. Datos datos.txt") 
abline(h=0) 
grid()

##-----------------------------------------------------------------------------------
##
## (b)	Aplicando m�todos no param�tricos, "�Hay alguna diferencia en los tiempos de 
##      reca�da entre los dos grupos de tratamiento (treat)?  
##      �Qu� tratamiento resulta m�s adecuado?"
##
##-----------------------------------------------------------------------------------

##  Comparaci�n de curvas de supervivencia
##
##   Orden para realizar la comparaci�n:
##      survdiff(Surv(variable_tiempo,variable_estado)~var_grupo,data=fichero,rho=valor)
##
##   Par�metro "rho=valor" indica el test de comparaci�n que se ha de aplicar, perteneciente 
##   a la familia de tests de Fleming-Harrington con pesos en cada instante de fallo de 
##   [SKM(t)]^rho siendo SKM(t) el estimador de K-M de la supervivencia  
##         - rho=0 (por defecto): test de log-rangos o de mantel-Haenszel.
##         - rho=1 : Test de Peto (modificaci�n del test de Wilcoxon. 
##

## Ajuste de la supervivencia de dos grupos

table(datos$race)
km2.fit <- survfit(Surv(time,censor)~race,data=datos, type="kaplan-meier")
summary(km2.fit)
km2.fit 

## Plot de las curvas K-M para ambos grupos 

plot(km2.fit,conf.int=FALSE,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Proporci�n sin reincidencia",lab=c(10, 10, 7),lwd=2, col=c("red","blue"),
     main="Estimaci�n K-M de S(t). Datos uissurv",lty=1:2)
legend("topright",lty=1:2,lwd=2,col=c("red","blue"),legend=c("Non-White","White"))
abline(h=0) 
grid()

##  Test de comparaci�n de la supervivencia entre ambos grupos. 
survdiff(Surv(time,censor)~race,data=datos,rho=1)

##-----------------------------------------------------------------------------------
##
## (c)	Obtener el modelo de Cox para dicha variable tiempo de supervivencia frente 
##      a las restantes variables, es decir: "age", "beck", "hercoc" , "ivhx", "ndrugtx",
##      "race", "treat",  "site", "los". 
##      
##-----------------------------------------------------------------------------------

# Orden para el estudio de un modelo de supervivencia: 
#                      Riesgos proporcionales de Cox (coxph) 
#   function (formula, data, weights, subset, na.action, init, control, 
#              ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
#              robust = FALSE, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
#              ...) 
# formula: objeto de supervivencia (creado por Surv()) en funci�n de ("~") variables predictoras
# data:    conjunto de datos
# weights: vector de ponderaciones, en su caso.
# subset : subconjunto de datos que debe ser considerado en el ajuste (�til para datos de
#          entrenamiento y datos de prueba)
# ties   : cadena de caracteres que especifica el m�todo para el caso de empates.
#          Si no hay empates en los tiempos de fallo todos los m�todos son equivalentes. 
#          Aunque casi todos los programas de regresi�n de Cox usan por defecto el m�todo de Breslow,
#          no es as� en este caso. Esta librer�a usa por defecto la aproximaci�n Efron,
#          (m�s precisa y eficiente computacionalmente). 
#          Es apropiado cuando los tiempos son un peque�o conjunto de valores discretos.
#          Las opciones son: "efron" (por defecto),"breslow","exact".
# method : Opci�n id�ntica a "ties"
# singular.ok :  valor l�gico que indica qu� hacer en caso de colinealidad en la matriz de modelo. 
#                Si TRUE (por defecto), el programa no considera autom�ticamente las columnas de la 
#                matriz X que son combinaciones lineales de las restantes columnas. 
#                En este caso, los coeficientes para tales columnas ser�n "NA", y la matriz de varianzas 
#                contienen ceros. 
#                Para los c�lculos auxiliares (predictor lineal,...) los coeficientes que faltan se 
#                tratan como ceros.
# (M�s detalles en el manual)
##

args(coxph)

names(datos)
time.cox1<- coxph(Surv(time, censor) ~ age + beck +  hercoc + ivhx + ndrugtx + race + treat + site + los , data=datos)
time.cox1

summary(time.cox1)
#  Se obtiene:
# - Tabla con coeficientes estimados, errores de estimaci�n y significaci�n de cada uno (z es el test
#   de Wald, asint�ticamente normal bajo la hipotesis de nulidad del coeficiente)
# - Los "exp(coef) tienen interpretaci�n como efectos multiplicativos sobre la funci�n de riesgo h(t).
#   As�, por ejemplo, manteniendo las otras variables constantes, un a�o adicional de la edad reduce 
#   el riesgo semanal de fallo por un factor de exp(-0.0574)= 0.944, en promedio, es decir, un 
#   (1-0.944)100%=5.6% . An�logamente, cada condena previa aumenta el riesgo con un factor de 1.096, 
#   o 9,6 por ciento.
# - El test de raz�n de verosimilitudes, test de Wald y test score son asint�ticamente equivalentes, 
#   y corresponden al test "�mnibus" con hip�tesis nula de que todas los coef. son 0. 
#

names(time.cox1)  # Para ver todo lo que contiene el objeto creado

time.cox1$coefficients # Estimaci�n de los coeficientes del modelo
time.cox1$var          # Estimaci�n de la matriz de covarianzas de los estimadores
time.cox1$loglik       # Log-verosimilitud del modelo y del modelo bajo H0:beta=0
time.cox1$score        # Test score
time.cox1$wald.test    # Test de Wald
time.cox1$linear.predictors[1:5]  # Valores "centrados" ajustados del predictor lineal 
                                  # (s�lo las 5 primeras)
time.cox1$means        # Medias de las variables predictoras
time.cox1$n            # N�mero de observaciones
time.cox1$nevent       # N�mero de eventos (fallos, muertes,...)
time.cox1$iter         # N�mero de iteraciones
time.cox1$method       # M�todo utilizado
time.cox1$residuals[1:5]  # Residuos "martingala". Para m�s detalles v�ase
                          # la orden residuals.coxph() y la referencia siguiente:
                          # T. Therneau, P. Grambsch, and T. Fleming. 
                          # "Martingale based residuals for survival models", Biometrika, March 1990.

##-----------------------------------------------------------------------------------
##
## (d)	Bas�ndose en los residuos de Schoenfeld, realizar el test sobre la hip�tesis   
##      de riesgos proporcionales del modelo de Cox.
##
##-----------------------------------------------------------------------------------

## Obtenci�n de los residuos: 
## Funci�n:
##    residuals(object, type=c("martingale", "deviance", "score", "schoenfeld",
##            "dfbeta", "dfbetas", "scaledsch","partial"),  ...)
##            
##
## object : modelo ajustado por la funci�n coxph
## type   : tipos de residuos solicitados con todas las posibilidades indicadas arriba
## (Otras opciones: v�ase el manual)
##
## - Los residuos de Schoenfeld se pueden considerar como los valores observados de las 
##   covariantes o variables explicativas menos los valores esperados  en cada instante de fallo:
##       x[i,k] - ajust.x[i,k]
##    
##    x[i,k] - valor observado en la variable X[k] en el sujeto o caso que "falla" en el instante t[i].
##    ajust.x[i,k] - valor esperado, seg�n el modelo ajustado, en la variable X[k] en el instante t[i]. 
##    
##       ajust.x[i,k] = media ponderada de la variable X[k] en los casos en riesgo R(t[i]) 
##                      en el instante t[i], instante en el que "falla" el caso i-�simo,
##                      donde los pesos son los valores de la funci�n de riesgo 
##                      en dichos casos para ese instante. 
##                      En concreto, si el caso j-�simo est� en R(t[i]), su funci�n de riesgo es:
##                       h(t[i] | x(j))= hbase(t[i]) exp{ x[j,1] beta[1]+...+x[j,p] beta[p] } 
##                      Por tanto, los pesos son: 
##                           w[i,j]= exp{ x[j,1] beta[1]+...+ x[j,p] beta[p] }
##                           para el caso j-�simo si est� incluido en R(t[i])
##
##   - Un valor positivo del residuo muestra un valor de X[k] m�s alto de lo esperado en ese instante de fallo.
##   - Los residuos de Schoenfeld para cada variable predictora suman cero.
##   - Se dispone de una matriz de residuos: una columna por cada variable explicativa, una fila por cada "fallo".
##   - Si estos residuos mantienen un patr�n aleatorio, es decir, no sistem�tico, proporciona una evidencia de 
##     que el efecto de la covariable no cambia respecto del tiempo, algo que presupone el modelo de Cox. 
##     Si hay alg�n tipo de patr�n sistem�tico, sugiere que el efecto de la covariable cambia a lo largo del tiempo.
##     As�, si es cierta la propiedad de riesgos proporcionales, los residuos no mostrar�n tendencias temporales 
##     y en el plot de los residuos frente al tiempo, la pendiente debe ser nula.
#

residuals(time.cox1,type=c("schoenfeld"))
res=residuals(time.cox1,type=c("schoenfeld"))

## Test sobre la hip�tesis de riesgos proporcionales del modelo de Cox
##  
##  Bas�ndose en los residuos de Schoenfeld, se consideran:
##  U[i]= rango de t[i]  en la colecci�n de tiempos de fallo
##  Se realiza el test de igualdad a cero del coeficiente de correlaci�n lineal
##  entre los rangos y los residuos asociados a cada una de las variables explicativas.
##  Si se rechaza la hip�teis nula, se debe concluir que se viola la hip�tesis de 
##  riesgos proporcionales.
##
##  Orden:
##          cox.zph(modeloajustado, transform="rank", global=TRUE)
##
##  - modeloajustado:   modelo de regresi�n de Cox ajustado a trav�s de la funci�n coxph
##  - transform     :   existen otras posibilidades, v�ase el manual.
##  - global        :   porporciona un test global conjunto.
##

cox.zph(time.cox1,transform = "rank")

# La variable "los" viola la hip�tesis de riesgos proporcionales.

## Test Gr�fico de riesgos proporcionales: plot.cox.zph 
##
## Gr�fico de los residuos de Schoenfeld (re-escalados) con el ajuste de un curva suavizada
## a trav�s de splines.
##
## Orden:  plot(x, resid=TRUE, se=TRUE, df=4, nsmo=40, var, ...)
##
## Argumentos:
## 
##   x     : modelo ajustado a trav�s de la funci�n cox.zph
##   resid : valor l�gico, TRUE se incluyen los residuos y la curva ajustada.
##   se    : valor l�gico, TRUE se incluyen bandas de confianza con 2 errores estandar 2se.
##   df    : grados de libertad de los splines naturales ajustados (df=2 ajuste lineal).
##   var   : conjunto de variables para las cuales se realizan los plots.

plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=1)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=2)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=3)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=4)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=5)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=6)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=7)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=8)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=9)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=10)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=11)
plot(cox.zph(time.cox1,transform = "rank"),df=2,col="red",var=12)

##-----------------------------------------------------------------------------------
##
## (e)	Ajustar de nuevo el modelo, eliminando las no significativas y las que provocan 
##      la violaci�n de la hip�tesis de riesgos proporcionales y evaluar nuevamente la 
##      hip�tesis de proporcionalidad. 
##
##-----------------------------------------------------------------------------------

time.cox2<- coxph(Surv(time, censor) ~ age + ivhx + ndrugtx + race + site, 
                 data=datos)
time.cox2

summary(time.cox2)
cox.zph(time.cox2,transform = "rank")
summary(datos$ivhx)

# C�lculo de la funci�n de supervivencia para un modelo de riesgos proporcionales de Cox
#
# Orden: survfit(objeto, newdata, se.fit=TRUE, conf.int=.95, individual=FALSE, type,vartype,
#                conf.type=c("log","log-log","plain","none"), censor=TRUE, id, 
#               newstrata, na.action=na.pass, ...)
#
# objeto      : objeto creado con la orden coxph()
# newdata     : data.frame con los mismos nombres de variables que aparecen en la f�rmula coxph. 
#               La curva(s) producida ser� representantiva de una cohorte cuyas covariables corresponden 
#               a los valores en newdata. 
#               El valor predeterminado es la media de las covariables utilizados en el ajuste coxph.
# conf.int    : Nivel para el intervalo de confianza para las curvas de supervivencia
#               Default is 0.95.
# se.fit      : Valor logico inidcando si los errores estandar deben  ser obtenidos o no
#               Por defecto, TRUE.
# conf.type   : Una de las opciones: "none", "plain", "log" (por defecto), o "log-log".
#               La opci�n "none" no calcula los intervalos de confianza. 
#               La segunda proporciona la curva de intervalos est�ndar +- k * se(curva), donde 
#               k se determina a partir conf.int. 
#               La opci�n "log" calcula intervalos basados en el riesgo acumulado o log(supervivencia).
#               La opci�n "log-log" determina intervalos basasdos en la transformaci�n
#               logar�tmica de la funci�n de riesgo o log(-log(supervivencia)).
# type,vartype : Cadena de caracteres especificando el tipo de curva de supervivencia. 
#               Los valores posibles son "aalen", "efron", o "kalbfleish-aprendiz" 
#               (s�lo los dos primeros caracteres son necesarios: aa, ef, ka). 
#               El valor predeterminado lo hace coincidir con el c�lculo utilizado en el modelo 
#               de Cox en caso de empates. Es decir,  
#                 - La estimaci�n de Nelson-Aalen-Breslow para " ties =  'Breslow' ", 
#                 - La estimaci�n de Efron para " ties = 'efron' "
#                 - La estimaci�n de Kalbfleisch-Prentice para un modelo de tiempo discreto
#                   con " ties = 'exact' ".
#               Estimaciones de la varianza son: Aalen-Link-Tsiatis, Efron, y Greenwood. 
#               El valor por defecto ser� la estimaci�n Efron para "ties = 'efron' " y 
#               la estimaci�n de Aalen en caso contrario.

# funci�n de supervivencia para valores medios de las predictoras
fsuperv.time2 = survfit(time.cox2)
names(fsuperv.time2)

fsuperv.time2$n          # tama�o muestral
fsuperv.time2$time       # tiempos observados
fsuperv.time2$n.risk     # casos en riesgos en cada tiempo
fsuperv.time2$n.event    # n�mero de eventos en cada instante de tiempo
fsuperv.time2$n.censor   # n�mero de censurados en cada instante de tiempo
fsuperv.time2$surv       # funci�n de supervivencia en cda instante de tiempo para 
                         # los casos con las covariantes indicadas.
fsuperv.time2$type       # tipo de censura (derecha, izquierda,..)
fsuperv.time2$std.err    # Error estandar de estimaci�n
fsuperv.time2$conf.int   # Nivel de confianza de los intervalos de confianza
fsuperv.time2$conf.type  # Aproximaci�n utilizada en el calculo del intervalo de confianza
fsuperv.time2$upper      # Limite superior del intervalo de confianza
fsuperv.time2$lower      # Limite inferior del intervalo de confianza

# Tabla de valores de supervivencias estimados e intervalos de confianza
t = cbind ("Dia"=fsuperv.time2$time,"Int. Inf"=fsuperv.time2$lower,"Func. Sup."=fsuperv.time2$surv, "Int. Sup"=fsuperv.time2$upper)
t

# Plot de la funci�n de supervivencia para valores medios de las predictoras
plot(survfit(time.cox2), xlab="d�as", ylab="Proporci�n de no reincidentes")

##-----------------------------------------------------------------------------------
##
##  (f)	�Qu� diferencias existen en el comportamiento de la variable tiempo hasta 
##     la reca�da entre dos individuos, uno con raza ="White" y otro con raza "Non-White"? 
##     Interpretar lo indicado por el modelo al respecto, seg�n las estimaciones de los
##     par�metros. Para ello, considerar dos casos con valores medios en "age", "beck", 
##     "ndrugtx" y "los", con hercoc="Heroin & Cocain, ivhx="Previous", treat="Long", site="B". 
##      
##-----------------------------------------------------------------------------------
names(datos)
mean(datos$los)

datos.fin <- with(datos, data.frame(
                  id=c(1001,1002), age  = rep(time.cox1$means[1], 2),     
                  beck = rep(time.cox1$means[2], 2), 
                  ndrugtx = rep(time.cox1$means[8], 2), 
                  hercoc  = c("Heroin & Cocain", "Heroin & Cocain"), 
                  ivhx = c("Previous","Previous"), 
                  race = c("White","Non-White"), 
                  treat = c("Long","Long"), site = c("B", "B"), los = rep(time.cox1$means[12], 2)))
datos.fin

# Plot de la funci�n de supervivencia
plot(survfit(time.cox2, newdata=datos.fin), conf.int=FALSE,
     lty=c(1, 2), ylim=c(0.0, 1), xlab="semanas",
     ylab="Proporci�n de no reincidentes",col = c(2,3))
legend("topright", legend=c("race = White", "race = Non-White"), lty=c(1 ,2), col = c(2,3),inset=0.02)

# Plot de la funci�n de riesgo acumulado
plot(survfit(time.cox2, newdata=datos.fin),conf.int=FALSE,
     lty=c(1, 2), 
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),lwd=2,fun="cumhaz", col = c(2,3),
     main="Estimaci�n de la tasa de fallo acumulada") 
legend("bottomright", legend=c("race = White", "race = Non-White"), lty=c(1 ,2), col = c(2,3),inset=0.02)
abline(h=0) 
grid()

##-----------------------------------------------------------------------------------
##
##  (g)	Interpretar lo indicado por el modelo, seg�n las estimaciones de los par�metros, 
##      sobre la variable ndrugtx n�mero de tratamientos previos. 
##-----------------------------------------------------------------------------------

summary(time.cox2)
names(time.cox2)
coefficients(time.cox2)[4]   # Coeficiente de la variable
exp(coefficients(time.cox2)[4]) # 
# Sobre la funci�n de riesgo h(t;x), el aumento de 1 unidad de la variable,
# provoca un cambio en dicha funci�n de:
#  h(t;x+1)=exp(coef) h(t;x)

##-----------------------------------------------------------------------------------
##
##  (h)	Obtener una estimaci�n de las funciones bases del modelo.
##-----------------------------------------------------------------------------------

# C�lculo de la funci�n de riesgo "base" ACUMULADA del modelo de Cox
# Orden 
#          basehaz(fit, centered = TRUE)
#  fit      : objeto ajustado por coxph
#  centered : Si es TRUE, la curva resultante es para un sujeto hipot�tico cuyos valores 
#             en las variables predictoras son las medias correspondientes de los datos 
#             originales, de lo contrario para un sujeto hipot�tico con valor nulo

basehaz(time.cox2, centered = TRUE)
plot(basehaz(time.cox2, centered = TRUE),type="l",col="2",
     ylab="Tiempo hasta reincidencia (en semanas)",
     main="Funci�n de riesgo ACUMULADA")
abline(h=0) 
grid()

plot(basehaz(time.cox2, centered = TRUE)$time,
     basehaz(time.cox2, centered = TRUE)$hazard, type="l",col="2",
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Riesgo base",
     main="Funci�n de riesgo ACUMULADA")
abline(h=0) 
grid()

plot(basehaz(time.cox2, centered = FALSE )$time,
     basehaz(time.cox2, centered = FALSE)$hazard, type="l",col="2",
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Riesgo base",
     main="Funci�n de riesgo base ACUMULADA. Datos datos.txt")
abline(h=0) 
grid()
