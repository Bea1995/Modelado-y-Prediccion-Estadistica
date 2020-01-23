#Modificamos las columnas 11 y 12 de wine.csv
Date<-read.csv("wine.csv",header=TRUE)
aux<-Date[,11]
Date[,11]<-Date[,12]
Date[,12]<-aux
write.csv(Date, file="vino.csv",row.names = FALSE)