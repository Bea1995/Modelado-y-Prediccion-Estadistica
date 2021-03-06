#__________________________________________________________________________________
#                     Modelado y Predicci�n Estad�stica.                       
#                M�ster en Matem�ticas. Universidad de Sevilla                
#                           Juan Manuel Mu�oz Pichardo                        
#                                                                             
#__________________________________________________________________________________
#
#                AN�LISIS DE SUPERVIVENCIA      
#                SCRIPT   cox01                                   
#                Estimaci�n de Kaplan Meier                       
#                Regresi�n de Cox Orden:                          
#                Datos: Rossi.txt                                 
#__________________________________________________________________________________

# CARGA DEL PAQUETE survival
library(survival)

# Algunas opciones de los c�lculos
options(digits = 4) # n�mero de decimales en las salidas
options(columns = 40) # numero de columnas en las salidas

# LECTURA DE DATOS
# El archivo Rossi.txt contiene los datos de un estudio experimental de la reincidencia 
# de 432 reclusos varones, que fueron observados durante un a�o despu�s de haber sido
# liberados de la prisi�n. 
#      Rossi, P. H., Berk, R. A., and Lenihan, K. J. (1980). Money, Work and Crime: 
#      Some Experimental Results. Academic Press, New York.
# La ilustraci�n se obtiene de: 
#      J, Fox & S. Weisberg (2011). Cox Proportional-Hazards Regression for Survival 
#      Data in R. An Appendix to "An R Companion to Applied Regression", Second Edition
# Las siguientes variables se incluyen en los datos.
#
## Variables
#
# week :    semana del primer arresto despu�s de la liberaci�n, o tiempo censurado.
# arrest :  Indicador de evento, igual a 1 para los detenidos durante el per�odo del estudio 
#            y 0 para los que no fueron detenidos (censurados).
# fin  :    Un factor, con niveles de SI, si el individuo recibi� ayuda financiera despu�s 
#           de salir de la c�rcel, y NO en caso contrario; la ayuda econ�mica fue un factor 
#           asignado aleatoriamente por los investigadores.
# age :     Edad en los a�os en el momento de la liberaci�n.
# race      Raza con niveles: Raza Negra y Otras razas.
# wexp:     Factor con los niveles de SI, si el individuo ten�a experiencia de trabajo a
#           tiempo completo antes de su encarcelamiento y No en caso contrario.
# mar:      Factor con los niveles "married" si el individuo estaba casado en el 
#           momento de la liberaci�n y "no.married" en caso contrario.
# paro :    Factor codificado Si, si el individuo fue puesto en libertad condicional y No en c.c.
# prio:     n�mero de condenas anteriores.
# educ:     Educaci�n, una variable categ�rica codificada num�ricamente, con los c�digos 
#           2 (grado 6 o menos), 3 (grados 6 a 9), 4 (de los grados 10 y 11), 5 (grado 12), 
#           � 6 (algunos postsecundaria )
# emp1-emp52 : Factores: SI, si el individuo estaba empleado en la semana correspondiente del 
#            estudio y NO, en caso contrario
#-------------------------------------------------------------

Rossi <- read.table("Rossi.txt", header=TRUE)
dim(Rossi)
Rossi[1:5,]
summary(Rossi)

#______________________________________________________________________
#
#  ACCI�N 1. ESTIMACI�N DE FUNCION DE SUPERVIVENCIA DE KAPLAN-MEIER. 
#            COMPARACI�N ENTRE GRUPOS 
#______________________________________________________________________

##----------------------------------------------------------------------------
#  1.1. Obtener la estimaci�n de Kaplan-Meier de la Funci�n de Supervivencia 
#  para el tiempo hasta el primer arresto despu�s de la liberaci�n con su
#  correspondiente error de estimaci�n y un intervalo de confianza al 95%.
##-----------------------------------------------------------------------------

# F.Superv.Emp�rica (notaci�n fse.fit)
#     Con  Surv(THR,status) se crea un objeto "de supervivencia" para esta librer�a
#     Con  survfit se ajusta la funci�n de supervivencia emp�rica 
#     Los Errores de Estimaci�n (ES) se obtienen seg�n la f�rmula de Greenwood
#        
Rossi$nocens=1  # Creamos una variable ficticia de censura para incluir todos
                # los datos como no censurados en el c�lculo
Rossi[1:5,]     # Comprobaci�n de la creaci�n de la variable en los 5 primeros casos
fse.fit <- survfit(Surv(week,nocens)~1,data=Rossi)
attributes(Surv(Rossi$week,Rossi$arrest))
attributes(fse.fit)
summary(fse.fit)     

## Representaci�n gr�fica de FSE

plot(fse.fit,conf.int=FALSE,xlab="Tiempo (en semanas)",
     ylab="Proporci�n sin reincidencia",lab=c(10,10,7),lwd=2)
mtext("Funci�n de supervivencia emp�rica para Tiempo hasta reincidencia",3,-1)
legend(5,.20,c("Datos: Rossi.txt","Sin tener en cuenta la censura"))
abline(h=0)
grid()

#_____________________________________________________________________________
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
 
km.fit <- survfit(Surv(week,arrest)~1,type="kaplan-meier",data = Rossi)
km.fit
summary(km.fit)
attributes(km.fit)
summary(km.fit)$time #los tiempos no censurados
summary(km.fit)$n.censor
tiempos=summary(km.fit)$time

## Representaci�n gr�fica de FS estimada por K-M frente a la FSE (emp�rica)
plot(km.fit,conf.int=FALSE,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Proporci�n sin reincidencia",lab=c(10, 10, 7),col="blue",lwd=2,
     main="Estimaci�n de S(t) datos Rossi.txt")
abline(h=0) 
grid()
lines(fse.fit,conf.int=FALSE,col="red",lty=2,lwd=2)
legend("topright",lty=1:2,col=c("blue","red"),lwd=2, legend=c("K-M","FSE"))

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
     main="Estimaci�n K-M de S(t). Datos Rossi.txt")
abline(h=0) 
grid()

plot(km.fit,conf.int=T,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),col="blue",lwd=2,fun="cumhaz",
     main="Estimaci�n K-M de la tasa de fallo acumulada. Datos Rossi.txt") 
abline(h=0) 
grid()

#-------------------------------------------------------------------------------------------
#
#  1.2. Comparar las curvas de supervivencia correspondientes a los grupos
#  seg�n "race".
#
##-------------------------------------------------------------------------------------------

##
##  Comparaci�n de curvas de supervivencia
##
##  Se pretende comparar las curvas de supervivencia de los grupos definidos por Race 
##  Black /No Black.
##  
##  Orden para realizar la comparaci�n:
##      survdiff(Surv(variable_tiempo,variable_estado)~var_grupo,data=fichero,rho=valor)
##
##   Par�metro "rho=valor" indica el test de comparaci�n que se ha de aplicar, perteneciente 
##   a la familia de tests de Fleming-Harrington con pesos en cada instante de fallo de 
##   [SKM(t)]^rho siendo SKM(t) el estimador de K-M de la supervivencia  
##         - rho=0 (por defecto): test de log-rangos o de mantel-Haenszel.
##         - rho=1 : Test de Peto (modificaci�n del test de Wilcoxon. 
##

## Ajuste de la supervivencia de dos grupos (seg�n raza)
##
##  Se pretende comparar las curvas de supervivencia de los grupos definidos por Race 
##  Black /No Black (1/0)

summary(Rossi$race)  # Est� considerada como num�rica
table(Rossi$race)

km2.fit <- survfit(Surv(week,arrest)~race,data=Rossi, type="kaplan-meier")
summary(km2.fit)
km2.fit 

## Plot de las curvas K-M para ambos grupos 
plot(km2.fit,conf.int=FALSE,xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Proporci�n sin reincidencia",lab=c(10, 10, 7),lwd=2, col=c("red","blue"),
     main="Estimaci�n K-M de S(t). Datos Rossi",lty=1:2)
legend("topright",lty=1:2,lwd=2,col=c("red","blue"),legend=c("No black","Black"))
abline(h=0) 
grid()

##  Test de comparaci�n de la supervivencia entre ambos grupos. 
survdiff(Surv(week,arrest)~race,data=Rossi,rho=1)

## ADICIONAL: AJUSTES SEG�N LOS CUANTILES SOBRE DIVERSOS MODELOS  
##-------------------------------------------------------------------------

# DISTRIBUCION EXPONENCIAL: 
#  (log(-log(S(t)) , log(t)) recta con pendiente 1 e intercept -ln(lambda)
# DISTRIBUCION WEIBULL:
#  (log(-log(S(t)) , log(t)) recta con pendiente 1/alpha e intercept -ln(lambda)

summary( lm(log(tiempos) ~ log(-log(superv)) ))
intercepto= coef(lm(log(tiempos) ~ log(-log(superv))))[1]
pendiente = coef(lm(log(tiempos) ~ log(-log(superv))))[2]  
pendiente
confint(lm(log(tiempos) ~ log(-log(superv)) ), level=0.95)
plot(log(-log(superv)),log(tiempos),main="Plot para el ajuste de modelos Exponencial o Weibull ")
abline(a=intercepto,b=pendiente,col=2)

# DISTRIBUCION tipo GUMBPEL: 
#  (log(-log(S(t)) , t) recta con pendiente -sigma e intercept mu
summary( lm(tiempos ~ log(-log(superv)) ))
intercepto= coef(lm(tiempos ~ log(-log(superv))))[1]
pendiente = coef(lm(tiempos ~ log(-log(superv))))[2]  
confint(lm(tiempos ~ log(-log(superv)) ), level=0.95)
plot(log(-log(superv)),tiempos,main="Plot para el ajuste de modelos Gumbel ")
abline(a=intercepto,b=pendiente,col=2)

# DISTRIBUCION LOG-LOG�STICA:
#  (-log[ S(t)/ (1-S(t))] , log(t)) recta con pendiente 1/alpha e intercept -ln(lambda)
logcoc=-log( superv / (1-superv))
summary( lm( log(tiempos) ~ logcoc ))
intercepto= coef(lm(log(tiempos) ~ logcoc))[1]
pendiente = coef(lm(log(tiempos) ~ logcoc))[2]  
confint(lm(log(tiempos) ~ logcoc ), level=0.95)
plot(logcoc,log(tiempos),main="Plot para el ajuste de modelo Log-Logistica ")
abline(a=intercepto,b=pendiente,col=2)

#__________________________________________________________________________________
#
# ACCI�N 2. MODELO DE RIESGOS PROPORCIONALES DE COX
#
#__________________________________________________________________________________
#
# 2.1. Analizar el modelo de Cox para dicha variable tiempo de supervivencia frente a
# las variables: 
#           fin, age, race, wexp, mar, paro y prio.
#
#__________________________________________________________________________________

args(coxph) # Argumentos de la funci�n coxph

# Orden para el estudio de un modelo de supervivencia: Riesgos proporcionales de Cox
#
#   function (formula, data, weights, subset, na.action, init, control, 
#              ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
#              robust = FALSE, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
#              ...) 
#
# formula: objeto de supervivencia (creado por Surv()) en funci�n de ("~") variables 
#          predictoras
# data:    conjunto de datos
# weights: vector de ponderaciones, en su caso.
#
# subset : subconjunto de datos que debe ser considerado en el ajuste (�til para datos 
#          de entrenamiento y datos de prueba)
#
# ties   : cadena de caracteres que especifica el m�todo para el caso de empates.
#          Si no hay empates en los tiempos de fallo todos los m�todos son equivalentes. 
#          Aunque casi todos los programas de regresi�n de Cox usan por defecto 
#          el m�todo de Breslow, no es as� en este caso. Esta librer�a usa por 
#          defecto la aproximaci�n Efron (m�s precisa y eficiente computacionalmente). 
#          Es apropiado cuando los tiempos son un peque�o conjunto de valores discretos.
#          Las opciones son: "efron" (por defecto),"breslow","exact".
#
# method : Opci�n id�ntica a "ties"
#
# singular.ok :  valor l�gico que indica qu� hacer en caso de colinealidad 
#                en la matriz del modelo. 
#                Si TRUE (por defecto), el programa no considera autom�ticamente las 
#                columnas de la matriz X que son combinaciones lineales
#                de las restantes columnas. 
#                En este caso, los coeficientes para tales columnas ser�n "NA", y 
#                la matriz de varianzas contiene ceros. 
#                Para los c�lculos auxiliares (predictor lineal,...) los coeficientes  
#                que faltan se tratan como ceros.
# (M�s detalles en el manual)
#---------------------------------------------------------------------------------------------------------

week.cox1<- coxph(Surv(week, arrest) ~ fin + age + race + wexp + mar + paro + prio, data=Rossi)
week.cox1

summary(week.cox1)
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

names(week.cox1)  # Para ver todo lo que contiene el objeto creado

week.cox1$coefficients # Estimaci�n de los coeficientes del modelo
week.cox1$var          # Estimaci�n de la matriz de covarianzas de los estimadores
week.cox1$loglik       # Log-verosimilitud del modelo y del modelo bajo H0:beta=0
week.cox1$score        # Test score
week.cox1$wald.test    # Test de Wald
week.cox1$linear.predictors[1:5]  # Valores "centrados" ajustados del predictor lineal 
                                  # (s�lo las 5 primeras)
week.cox1$means        # Medias de las variables predictoras
week.cox1$n            # N�mero de observaciones
week.cox1$nevent       # N�mero de eventos (fallos, muertes,...)
week.cox1$iter         # N�mero de iteraciones
week.cox1$method       # M�todo utilizado
week.cox1$residuals[1:5]  # Residuos "martingala". Para m�s detalles v�ase
                          # la orden residuals.coxph() y la referencia siguiente:
                          # T. Therneau, P. Grambsch, and T. Fleming. 
                          # "Martingale based residuals for survival models", Biometrika, March 1990.

##--------------- Obtenci�n de los residuos: 
##
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
##--------------------------------------------------------------------------------

residuals(week.cox1,type=c("schoenfeld"))
res=residuals(week.cox1,type=c("schoenfeld"))

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

cox.zph(week.cox1,transform = "rank")

# Las variables "age" y "wexp" violan la hip�tesis de riesgos proporcionales.

#--------------   Test Gr�fico de riesgos proporcionales: 
##    plot.cox.zph(...)
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

plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=1)
plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=2)
plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=3)
plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=4)
plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=5)
plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=6)
plot(cox.zph(week.cox1,transform = "rank"),df=2,col="red",var=7)

# Cambio de modelo: se eliminan las variables no significativas, 
#                   Evaluar nuevamente la hip�tesis de proporcionalidad
week.cox2<- coxph(Surv(week, arrest) ~ fin + prio, data=Rossi)
week.cox2

summary(week.cox2)
cox.zph(week.cox2,transform = "rank")

# Adem�s del ajuste, es interesante examinar la distribuci�n estimada de los tiempos de 
# supervivencia. Las estimaciones de la funci�n de supervivencia S(t), por defecto en 
# los valores medios de las covariables.

#--------------------------------------------------------------
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
#------------------------------------------------------------------------------------

# funci�n de supervivencia para valores medios de las predictoras
fsuperv.week1 = survfit(week.cox2)
names(fsuperv.week1)

fsuperv.week1$n          # tama�o muestral
fsuperv.week1$time       # tiempos observados
fsuperv.week1$n.risk     # casos en riesgos en cada tiempo
fsuperv.week1$n.event    # n�mero de eventos en cada instante de tiempo
fsuperv.week1$n.censor   # n�mero de censurados en cada instante de tiempo
fsuperv.week1$surv       # funci�n de supervivencia en cda instante de tiempo para 
                         # los casos con las covariantes indicadas.
fsuperv.week1$type       # tipo de censura (derecha, izquierda,..)
fsuperv.week1$std.err    # Error estandar de estimaci�n
fsuperv.week1$conf.int   # Nivel de confianza de los intervalos de confianza
fsuperv.week1$conf.type  # Aproximaci�n utilizada en el calculo del intervalo de confianza
fsuperv.week1$upper      # Limite superior del intervalo de confianza
fsuperv.week1$lower      # Limite inferior del intervalo de confianza

# Tabla de valores de supervivencias estimados e intervalos de confianza
t = cbind (fsuperv.week1$lower,fsuperv.week1$surv, fsuperv.week1$upper)
t

# Plot de la funci�n de supervivencia para valores medios de las predictoras
plot(survfit(week.cox1), ylim=c(0.7, 1), xlab="semanas", 
     ylab="Proporci�n de no reincidentes")

##_____________________________________________________________________________________
##
##     2.2.En dicho modelo, realizar un an�lisis particular de la variable "fin", 
##     representado las funciones de supervivencia estimadas para los casos con
##     financiaci�n y sin financiaci�n, as� como la funci�n de riesgo acumulada.
##
##_____________________________________________________________________________________

#Dos individuos iguales en todos salvo que uno recibe financiaci�n y otro no
Rossi.fin <- with(Rossi, data.frame(fin=c(0, 1),
                                    age=rep(mean(age), 2), race=rep(mean(race == "other"), 2),
                                    wexp=rep(mean(wexp == "yes"), 2), mar=rep(mean(mar == "not married"), 2),
                                    paro=rep(mean(paro == "yes"), 2), prio=rep(mean(prio), 2)))
Rossi.fin
week.cox2$means

# Plot de la funci�n de supervivencia
plot(survfit(week.cox2, newdata=Rossi.fin), conf.int=FALSE,
     lty=c(1, 2), ylim=c(0.6, 1), xlab="semanas",
     ylab="Proporci�n de no reincidentes",col = c(2,3))
legend("bottomleft", legend=c("fin = No", "fin = Si"), lty=c(1 ,2), col = c(2,3),inset=0.02)

# Plot de la funci�n de riesgo acumulado
plot(survfit(week.cox1, newdata=Rossi.fin),conf.int=FALSE,
     lty=c(1, 2), 
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="H(t)",lab=c(10, 10, 7),lwd=2,fun="cumhaz", col = c(2,3),
     main="Estimaci�n de la tasa de fallo acumulada. Datos Rossi.txt") 
legend("topleft", legend=c("fin = No", "fin = Si"), lty=c(1 ,2), col = c(2,3),inset=0.02)
abline(h=0) 
grid()

##_____________________________________________________________________________________
##
##  2.3 Obtener una estimaci�n de las funci�n base del modelo
##
##_____________________________________________________________________________________

# C�lculo de la funci�n de riesgo "base" ACUMULADA del modelo de Cox
# Orden 
#          basehaz(fit, centered = TRUE)
#  fit      : objeto ajustado por coxph
#  centered : Si es TRUE, la curva resultante es para un sujeto hipot�tico cuyos valores 
#             en las variables predictoras son las medias correspondientes de los datos 
#             originales, de lo contrario para un sujeto hipot�tico con valor nulo

basehaz(week.cox2, centered = TRUE)

plot(basehaz(week.cox2, centered = TRUE),type="l",col="4",
     ylab="Tiempo hasta reincidencia (en semanas)",
     main="Func.riesgo base ACUMULADA para casos promedios. Datos Rossi.txt")
abline(h=0) 
grid()
# Con la orden anterior, realiza la representaci�n tiempo(eje vertical) frente
# a funci�n riesgo (eje horizontal)

# Para obtenerlo en forma transpuesta:
plot(basehaz(week.cox2, centered = TRUE)$time,
     basehaz(week.cox2, centered = TRUE)$hazard, type="l",col="2",
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Riesgo base",
     main="Funci�n de riesgo ACUMULADA para casos promedios. Datos Rossi.txt")
abline(h=0) 
grid()

plot(basehaz(week.cox2, centered = FALSE )$time,
     basehaz(week.cox2, centered = FALSE)$hazard, type="l",col="4",
     xlab="Tiempo hasta reincidencia (en semanas)",
     ylab="Riesgo base",
     main="Funci�n de riesgo base ACUMULADA. Datos Rossi.txt")
abline(h=0) 
grid()
