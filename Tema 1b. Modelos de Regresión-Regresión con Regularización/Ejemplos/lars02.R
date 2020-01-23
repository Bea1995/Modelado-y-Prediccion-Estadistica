#__________________________________________________________________________________
###                     Modelado y Predicción Estadística.                      ### 
###                Máster en Matemáticas. Universidad de Sevilla                ###
###                           Juan Manuel Muñoz Pichardo                        ###
###                                                                             ###
#__________________________________________________________________________________

#########                SCRIPT  lars02.R                                 #########
#########                Regresión Ridge y Lars                           #########
#########                Orden: lm.ridge; lars                            #########
#########                Datos: diabetes (con formato csv)                #########
#__________________________________________________________________________________

######### CARGA DEL PAQUETE lars
library(MASS)
install.packages("lars")
library(lars)

######### Algunas opciones de los cálculos
options(digits = 4) # número de decimales en las salidas
options(columns = 40) # numero de columnas en las salidas
#__________________________________________________________________________________
#                           LECTURA DE DATOS "diabetes"
#
# Un dataframe con 442 casos y 11 "columnas". 
# Datos utilizados en el artículo:
# Efron, Hastie, Johnstone and Tibshirani (2003). "Least Angle Regression" (with discussion)
# Annals of Statistics
#
# El data frame contiene las siguientes columnas: 
# x : una matriz con 10 columnas (variables iniciales)
# y : un vector numeríco 
# x2: una matriz con 64 columnas (diversas interacciones de variablesd iniciales)
#
# Son el resultado de un estudio sobre 442 pacientes de diabetes a los que se les midieron:
# - diez variables iniciales,  incluidas en x. 
#   ** Edad, Sexo, Índice de masa corporal, Presión sanguínea media, 
#      y seis medidas de suero sanguíneos (tc,ldl, hdl,tch,ltg,glu).
# - Una medida de la progresión de la enfermedad un año después, incluida en y, considerada
#   como medida respuesta.
#__________________________________________________________________________________

diabetes = read.csv("diabetes.txt",sep=";",header=TRUE)  # read csv file 
dim(diabetes)
names(diabetes)
diabetes[1:5,]

#__________________________________________________________________________________
# Estandarización a norma 2 igual a 1 (dividir por raiz(n-1)*sd)
#    La función "scale" tipifica, es decir, transforma a media cero y varianza 1
#    por lo que divide por sd

diabet <- data.frame(scale(diabetes[,-c(2,11)])) # Eliminamos la variable SEX y la Y
# Comprobación de la función "scale"
diabet[1:3,]
summary(diabet)
sd(diabet$AGE)
sd(diabet$BMI)

n=dim(diabetes)[1]
diabet <- data.frame(scale(diabetes[,-c(2,11)])/sqrt(n-1))
# Comprobación del calculo
diabet[1:3,]
sd(diabet$AGE)
sd(diabet$BMI)

# Inclusión de la variable "Y"
y = diabetes$Y
diabet <- data.frame(diabet,y)
diabet[1:5,]
names(diabet)
dim(diabet)

#__________________________________________________________________________________
#                AJUSTE POR MÍNIMOS CUADRADOS ORDINARIOS
#
#   Ajuste por MCO de variable y frente al resto
#__________________________________________________________________________________
diabetes.mco <- lm(y ~ ., diabet)
summary(diabetes.mco)

#__________________________________________________________________________________
#                       MULTICOLINEALIDAD
#
# Estudiar la multicolinealidad con las correlaciones bivariadas, con el
# factor de inflación de la varianza (librería usdm) y con la matriz de
# correlaciones parciales entre las variables predictoras (paquete ggm,
# orden parcor(S), siendo S la matriz de covarianzas muestrales de las variables.
cor(diabet)

# FACTOR INFLACTOR DE LA VARIANZA
#
# El coeficiente VIF de una variable predictora X(j) es funcion del cuadrado del 
# coeficiente de correlación múltiple de X(j) respecto al resto de predictoras: R2(j).
# En concreto:
#          VIF[X(j)] = 1 / (1-R2(j))
#
# Así, si VIF es alto (superior a 10) se dice que hay multicolinealidad
install.packages("usdm")
library(usdm)
vif(diabet[,-10])

#__________________________________________________________________________________
#                           AJUSTE POR REGRESIÓN RIDGE 
#
# Ajuste por Regresión Ridge con valores de lambda variando 
# entre 0 y 0.1 (0.001, 0,002,...,0.100). 
# Para k=0 coincide con MCO 
diabet.rr <- lm.ridge(y ~ ., diabet,lambda = seq(0, 50, 0.01))
summary(diabet.rr)
names(diabet.rr)

# Matriz de coeficientes estimados (para cada valor de lambda por filas)
# La primera columna almacena el valor de lambda
coef(diabet.rr)[1:3, ]
coef(diabet.rr)[99:101, ]

# Media de la variable objetivo
diabet.rr$ym

# Media de las variables predictoras
diabet.rr$xm

# Vector de las "escalas" utilizadas en la estadarización de las variables
# predictoras (sus desviaciones típicas)
diabet.rr$scales

# Vector de valores de lambda
diabet.rr$lambda[1:3]

# Vector de Valores de errores por validación cruzada generalizada (GCV values)
diabet.rr$GCV[1:6]

# Valor del criterio HKM propuesto por Hoerl et al. (1975)  
diabet.rr$kHKB
# Valor del criterio HKM propuesto por Lawless and Wang.(1976) 
diabet.rr$kLW

# La función "select" aplicada al objeto que proporciona lm.ridge devuelve
# los valores óptimos de tres de los siguientes criterios:
#   - Criterio HKM propuesto por Hoerl et al. (1975) y tiene una justificación bayesiana. 
#   - Criterio LW propuesto en Lawless and Wang.(1976)
#   - Criterio de validación cruzada generalizada (GCV).
select(diabet.rr)

# Selección del lambda óptimo de acuerdo al criterio GCV y su regresión ridge
nGCV <- which.min(diabet.rr$GCV)  # número para el cual se alcanza el mínimo
nGCV
lGCV <- diabet.rr$lambda[nGCV]  # valor de lambda para el cual se alcanza el mínimo
lGCV

diabet.ridge <- lm.ridge(y ~ ., diabet, lambda = lGCV)
summary(diabet.ridge)
names(diabet.ridge)

coef(diabet.ridge)
diabet.ridge$scales
diabet.ridge$lambda
diabet.ridge$ym
diabet.ridge$xm
diabet.ridge$GCV

# Gráficos de la traza ridge y las medidas del ECM en función de lambda
matplot(diabet.rr$lambda, t(diabet.rr$coef), type = "l", xlab = expression(lambda),
         ylab = expression(beta[i]))
abline(v = lGCV)
mtext(expression(lambda[GCV]), side = 3, at = lGCV)
title(main = "Trazas ridge")

plot(diabet.rr$lambda, diabet.rr$GCV, type = "l", xlab = expression(lambda),
       ylab = "GCV", main = "Criterio GCV",col=2)
abline(v = lGCV)
mtext(expression(lambda[GCV]), side = 3, at = lGCV)
abline(v = diabet.rr$kLW)
mtext(expression(lambda[LW]), side = 3, at = diabet.rr$kLW)

#__________________________________________________________________________________
#                         AJUSTE POR REGRESIÓN LASSO
# 
# ORDEN:  lars(x, y, type = c("lasso", "lar", "forward.stagewise", "stepwise"),
#              trace = FALSE, normalize = TRUE, intercept = TRUE, max.steps, ...)  
#             
#
# x : MATRIZ DE PREDICTORES
# y : vector de respuestas
# type: Tipo de técnica de regresión LASSO, LAR, forward.stagewise, stepwise.
#       Por defecto lasso
# trace : Si TRUE guarda todo el proceso del algoritmo. Por defecto FALSE
# normalize : Si TRUE estandariza para pasar a norma L2 igual a 1 (dividir por raiz(n-1)*sd)
#             por defecto, TRUE
# intercept : Si TRUE (por defecto) incluye intercept en el modelo.
# max.steps : Número máximo de pasos. Por defecto 8*min(num.var , num.obs - 1)
#__________________________________________________________________________________

x <- as.matrix(diabet[,1:9])
y <- as.vector(diabet[,10])
diabet.las = lars(x,y)

names(diabet.las)
diabet.las$type # tipo de regresión realizada
diabet.las$df   # grados de libertad efectivo o estimado
diabet.las$lambda # Valores de lambda
diabet.las$R2 # Coeficiente R2
diabet.las$RSS # Suma de cuadrados del error
diabet.las$Cp  # Estadistico Cp
               # Estadistico Cp de Mallows, como ajuste del modelo. 
               # Se define como: Cp=[SCE(p) / CME ] - (n-2p)
               #  SCE(p) suma de cuadrados de los residuos del modelo con p parámetros
               #  CME : Error cuadrático Medio del modelo con todas las variables
diabet.las$beta # coeficientes
diabet.las$mu   # media de la variable objetivo
diabet.las$meanx   # media de las var. predictoras

summary.lars(diabet.las)

# Gráficos de la traza
# La orden 
#  plot(objeto, xvar= c("norm", "df", "step"), 
#       breaks = TRUE, plottype = c("coefficients","Cp"), 
#       omit.zeros = TRUE, eps = 1e-10, ...)
#
# Produce un plot de ajuste lars. El valor predeterminado es un plot de las trazas
# de los coeficientes 
# 
# Parámetros:
# objeto : objeto lars
# plottype : "coefficients" evolución de los coeficientes (en el eje vertical)
#            "Cp" evolución del estadístico Cp (en el eje vertical)
# xvar   : Parámetro o variable que determina el eje horizontal del plot.
#   xvar = norm (por defecto) Norma L1 del vector t=|| beta||, es decir el valor t
#                de la restricción del problema de optimización.
#   xvar = step  número de pasos (que es esencialmente grados de libertad efectivo para LAR, 
#                no para LASSO o Forward Stagewise).
#   xvar = df    df estimado
plot.lars(diabet.las)
plot.lars(diabet.las,plottype="Cp",xvar="step")
abline(h=min(diabet.las$Cp),col=4,lty=3,lwd=3)
mtext(expression("Cpmin"), side = 2, at = min(diabet.las$Cp))
abline(v = 7,col=2)

# Plot
plot(diabet.las$R2, col=4, lwd=3, type="l",main="Plot de los Coeficiente R2 ",ylab="R2")
abline(h=max(diabet.las$R2))
mtext(expression("R2max"), side = 2, at = max(diabet.las$R2))
abline(v = 7)

# Seleccionamos el paso 6 (índice 7)
diabet.las$R2[7]
max(diabet.las$R2)

plot.lars(diabet.las,breaks = FALSE, xvar="step")
abline(v = 6,col=4)
mtext(expression("con R2=0.5006"), side = 3, at = 6,col=4)
abline(v = 12,col=3)
mtext(expression("con R2=0.5177"), side = 3, at = 12,col=3)

coef(lm.ridge(y ~ ., diabet, lambda = lGCV)) #coeficientes de la regresión ridge
diabet.las$beta[7,] # coeficientes
