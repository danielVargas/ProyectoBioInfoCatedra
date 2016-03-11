##### QUITAR COMENTARIOS PARA INSTALAR LIBRERIAS
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("Biobase")
#biocLite("plyr")
#biocLite("siggenes")
#biocLite("limma")
####################
# Importación de librerias

library(Biobase)
library(GEOquery)
library(plyr)
library(siggenes)
library(limma)
library(amap)

####################


# CARGADO DE DATOS
db <- getGEO('GDS2881', destdir=".")
#db<-getGEO(filename='GDS2881.soft.gz')

# CARGADO DE METADADOS DEL DATASET
Meta(db)$dataset_id
Meta(db)$platform_organism
Meta(db)$feature_count
Meta(db)$sample_count

colnames(Table(db))


Table(db)[1:5,1:7]

matriz_log<- GDS2eSet(db, do.log2=TRUE)


muestras<-sampleNames(matriz_log)
genes<-featureNames(matriz_log)
expresion <- exprs(matriz_log)
matriz_expresion<-exprs(matriz_log[,])


### SE NOMBRAR FILAS Y COLUMNAS CON NOMBRES DE MUESTRAS Y GENES
colnames(matriz_expresion) <- pData(matriz_log)$"disease.state"
rownames(matriz_expresion) <- featureData(matriz_log)$"Gene symbol"
###########################################################################
### SE ELIMINAN VALORES NA Y ESPACIOS VACIOS
matriz_expresion <- na.omit(matriz_expresion) #NA
matriz_expresion <- matriz_expresion[rownames(matriz_expresion) != "",] #vacios
matriz_expresion[1:10,1:6]
##########################################################################
#Encuentra nombre de genes con / y reemplaza / por -
#los nombres de variables no pueden llevar / para evitar problemas de programación
rownames(matriz_expresion)<-gsub("/","-",rownames(matriz_expresion))
#promediar genes repetidos
matriz_expresion_promedio<-aggregate(matriz_expresion, by=list(rownames(matriz_expresion)),
                                     FUN=mean)
#Al promediar los nombres de genes se convierten en columnas
#por eso se deben reasignar como nombres de filas
matriz_expresion_promedio<- data.frame(matriz_expresion_promedio[,-1],
                                       row.names=matriz_expresion_promedio[,1])
colnames(matriz_expresion_promedio) <- pData(matriz_log)$"disease.state"
#rownames(matriz_expresion) <- featureData(matriz_log)$"Gene Symbol"
matriz_expresion_promedio[1:10,1:6]
#################################################################################

#########################################################################################################
#Hipótesis
#(H0): No hay capacidad de discriminación entre clase normal y clase cáncer.
#(H1): Hay capacidad de discriminación entre clase normal y clase cáncer.
#alfa=0,05
#Si p =< alfa se rechaza la hipótesis nula (H0).
#Si p > alfa no se rechaza la hipótesis nula (H0).
#Obtengo p.value desde Test t, variables: NORMAL(muestra 14-21) y CANCER (muestra 1-13)
#El 1, indica que se aplica a nivel de filas.
pvals=apply(matriz_expresion_promedio,1,function(x) {t.test(x[14:21],x[1:13])$p.value})
#Asignar p.value al gen correspondiente
matriz_expresion_promedio$p <- pvals
#Ordenar de - a + la matriz de acuerdo al p.value de cada gen.
ordenarValoresP<-matriz_expresion_promedio[order(pvals),]
matriz_expresion_discriminativos<-subset(ordenarValoresP, p<=0.05)
head(matriz_expresion_discriminativos)
# SE NOMRAN COLUMNAS NUEVAMENTE AHORA CON "P.VALUE"
colnames(matriz_expresion_discriminativos) <- c(array(pData(matriz_log)$"disease.state"),"p.value")
########################################################################################################


### LIMMMAZ
#El paquete limma contiene funciones para el uso de t-test y ANOVA para identificar 
#la expresión diferencial de los datos de microarrays. Estas funciones se pueden utilizar
#para las plataformas de toda gama y trabajan incluso para los datos de microarrays con 
#diseños complejos de múltiples muestras. La idea central se basa en ajustar un modelo lineal a 
#los datos de expresión de cada gen . Los datos de expresión pueden ser log-ratios o log-intensities. 
#Limma está diseñado para ser utilizado en combinación con el paquete affy y es compatible con muchos
#otros paquetes

## Limma necesita la información de phenodata que esta en la matriz_log extraida desde el dataset

ph = matriz_log@phenoData
ph@data[ ,2] = c(array(pData(matriz_log)$"disease.state"))
colnames(ph@data)[2]="source"
ph@data
groups = ph@data$source
f = factor(groups,levels=c("normal","stage I cRCC","stage II cRCC"))
## Se crea el modelo necesario para la implementación del método
design = model.matrix(~ 0 + f)
colnames(design) = c("normal","stageIcRCC","stageIIcRCC")
rownames(expresion) <- featureData(matriz_log)$"Gene symbol"
data.fit = lmFit(expresion,design)
data.fit$coefficients[1:10,]
## Se crea un segundo factor de prueba 
f2 = factor(groups,levels=c("normal","stageIcRCC","stageIIcRCC"))
design2 = model.matrix(~ 0 + f2) # y un segundo modelo de prueba
colnames(design2) = c("normal","stageIcRCC","stageIIcRCC")
#Se crea la matrix de contraste
contrast.matrix = makeContrasts(normal - stageIcRCC - stageIIcRCC,levels=design2)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
#Mediante el metodo eBayes se obtienen los genes diferencialmente expresados
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients[1:10,]

data.fit.eb$t[1:10,]

## Se trabajan loS datos para poder mostrarlos más ordenandamente
data.fit.eb$p.value[1:10,]
temp = data.matrix(data.fit.eb$p.value)
temp = temp[order(temp)]
temp[1:10]
temp = data.matrix(temp)
colnames(temp) <-  c("p-value")
temp[1:10]
summary(temp)
## Finalmene se seleciconan los genes más diferencialmente expresados.
temp =topTable(data.fit.eb, adjust = "fdr")
results <- decideTests(data.fit.eb)

## Y se arma una matriz para ingrearla al heatmap
matrix.genes.dif = c()
for (variable in temp$ID) {
  matrix.genes.dif = rbind(matrix.genes.dif, matriz_expresion_promedio[variable,])
  # print(matriz_expresion_promedio[variable,])
}
matrix.genes.dif<-na.omit(matrix.genes.dif[,1:21])

t = as.matrix(matrix.genes.dif)

colnames(t) <- c(groups,"p")

heat = as.matrix(matriz_expresion_discriminativos)
colnames(heat) <- c(groups,"p")
heat<-heat[ ,-21]
heatmap(heat, main = "Heatmap de genes")
t<-t[ ,-21]
t <- t[ , sort(groups)]
colnames(t) <- sort(groups)
heatmap(t, main = "Heatmap de genes")

##############################################################################################3
#### ARREGLOS SOLICITADOS POR EL PROFESOR
## y se aplica un algortimo de klustering para el análisis de los datos.

matriz_expresion_discriminativoskluster = as.matrix(matriz_expresion_discriminativos)
algoritmoKg <- Kmeans(x = matriz_expresion_discriminativoskluster, centers=4, method="kendall",iter.max = 100)
# using package ade4
library(ade4)
pca    <-prcomp(matriz_expresion_discriminativoskluster, scale.=T, retx=T)  # principal components analysis
plot.df <- cbind(pca$x[,1], pca$x[,2]) # first and second PC
coul <- c("black", "red", "green", "blue")
s.class(plot.df, factor(algoritmoKg$cluster),col = coul)


matrizCluster = as.matrix(algoritmoKg$cluster)
grupos = rownames(matrizCluster)
grupo1 = data.frame(c())
grupo1Names = c()
grupo2 = data.frame(c())
grupo2Names = c()
grupo3 = data.frame(c())
grupo3Names = c()
grupo4 = data.frame(c())
grupo4Names = c()
i=1
for (variable in matrizCluster){
   if (variable == 1){
      grupo1= rbind(grupo1, matriz_expresion_discriminativos[grupos[i],])
      grupo1Names = c(grupo1Names,grupos[i])
   }
  if (variable == 2){
    grupo2= rbind(grupo2, matriz_expresion_discriminativos[grupos[i],])
    grupo2Names = c(grupo2Names,grupos[i])
  }
  if (variable == 3){
    grupo3= rbind(grupo3, matriz_expresion_discriminativos[grupos[i],])
    grupo3Names = c(grupo3Names,grupos[i])
  }
  if (variable == 4){
    grupo4= rbind(grupo4, matriz_expresion_discriminativos[grupos[i],])
    grupo4Names = c(grupo4Names,grupos[i])
  }
  i= i +1
}

grupo1<-grupo1[ ,-21]
heat = as.matrix(grupo1)
heatmap(heat, main = "Heatmap de Grupo1")
grupo2<-grupo2[ ,-21]
heat = as.matrix(grupo2)
heatmap(heat, main = "Heatmap de Grupo2")
grupo3<-grupo3[ ,-21]
heat = as.matrix(grupo3)
heatmap(heat, main = "Heatmap de Grupo3")
grupo4<-grupo4[ ,-21]
heat = as.matrix(grupo4)
heatmap(heat, main = "Heatmap de Grupo4")

### NORMALIDAD
require(nortest)
testgrupo1 <- grupo1
plot(density(as.matrix(testgrupo1)), main="Grupo1")
testgrupo2 <- grupo2
plot(density(as.matrix(testgrupo2)), main="Grupo2")
testgrupo3 <- grupo3
plot(density(as.matrix(testgrupo3)), main="Grupo3")
testgrupo4 <- grupo4
plot(density(as.matrix(testgrupo4)), main="Grupo4")


########################################################################################################
#LIMMA GRUPO1
###########
groups = colnames(matriz_expresion_discriminativoskluster)
f = factor(groups,levels=c("normal","stage I cRCC","stage II cRCC"))
## Se crea el modelo necesario para la implementaci�n del m�toMCF7Areprodo
design = model.matrix(~ 0 + f)
colnames(design) = c("normal","stageIcRCC","stageIIcRCC")
colnames(grupo1) <- c(groups)

data.fit = lmFit(grupo1,design)
data.fit$coefficients[1:10,]

#####

## Se crea un segundo factor de prueba 
f2 = factor(groups,levels=c("normal","stageIcRCC","stageIIcRCC"))
design2 = model.matrix(~ 0 + f2) # y un segundo modelo de prueba
colnames(design2) = c("normal","stageIcRCC","stageIIcRCC")
#Se crea la matrix de contraste
contrast.matrix = makeContrasts(normal - stageIcRCC - stageIIcRCC,levels=design2)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
#Mediante el metodo eBayes se obtienen los genes diferencialmente expresados
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients[1:10,]

data.fit.eb$t[1:10,]


## Se trabajan loS datos para poder mostrarlos más ordenandamente
data.fit.eb$p.value[1:10,]
temp = data.matrix(data.fit.eb$p.value)
temp = temp[order(temp)]
temp[1:10]
temp = data.matrix(temp)
colnames(temp) <-  c("p-value")
temp[1:10]
summary(temp)
## Finalmene se seleciconan los genes más diferencialmente expresados.
temp =topTable(data.fit.eb, adjust = "fdr")
results <- decideTests(data.fit.eb)

#########################################################################################################

########################################################################################################
#LIMMA GRUPO2
###########
f = factor(groups,levels=c("normal","stage I cRCC","stage II cRCC"))
## Se crea el modelo necesario para la implementaci�n del m�toMCF7Areprodo
design = model.matrix(~ 0 + f)
colnames(design) = c("normal","stageIcRCC","stageIIcRCC")
colnames(grupo2) <- c(groups)

data.fit = lmFit(grupo2,design)
data.fit$coefficients[1:10,]

#####

## Se crea un segundo factor de prueba 
f2 = factor(groups,levels=c("normal","stageIcRCC","stageIIcRCC"))
design2 = model.matrix(~ 0 + f2) # y un segundo modelo de prueba
colnames(design2) = c("normal","stageIcRCC","stageIIcRCC")
#Se crea la matrix de contraste
contrast.matrix = makeContrasts(normal - stageIcRCC - stageIIcRCC,levels=design2)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
#Mediante el metodo eBayes se obtienen los genes diferencialmente expresados
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients[1:10,]

data.fit.eb$t[1:10,]


## Se trabajan loS datos para poder mostrarlos más ordenandamente
data.fit.eb$p.value[1:10,]
temp = data.matrix(data.fit.eb$p.value)
temp = temp[order(temp)]
temp[1:10]
temp = data.matrix(temp)
colnames(temp) <-  c("p-value")
temp[1:10]
summary(temp)
## Finalmene se seleciconan los genes más diferencialmente expresados.
temp =topTable(data.fit.eb, adjust = "fdr")
results <- decideTests(data.fit.eb)


#########################################################################################################

########################################################################################################
#LIMMA GRUPO3
###########
f = factor(groups,levels=c("normal","stage I cRCC","stage II cRCC"))
## Se crea el modelo necesario para la implementaci�n del m�toMCF7Areprodo
design = model.matrix(~ 0 + f)
colnames(design) = c("normal","stageIcRCC","stageIIcRCC")
colnames(grupo3) <- c(groups)

data.fit = lmFit(grupo3,design)
data.fit$coefficients[1:10,]

#####

## Se crea un segundo factor de prueba 
f2 = factor(groups,levels=c("normal","stageIcRCC","stageIIcRCC"))
design2 = model.matrix(~ 0 + f2) # y un segundo modelo de prueba
colnames(design2) = c("normal","stageIcRCC","stageIIcRCC")
#Se crea la matrix de contraste
contrast.matrix = makeContrasts(normal - stageIcRCC - stageIIcRCC,levels=design2)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
#Mediante el metodo eBayes se obtienen los genes diferencialmente expresados
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients[1:10,]

data.fit.eb$t[1:10,]


## Se trabajan loS datos para poder mostrarlos más ordenandamente
data.fit.eb$p.value[1:10,]
temp = data.matrix(data.fit.eb$p.value)
temp = temp[order(temp)]
temp[1:10]
temp = data.matrix(temp)
colnames(temp) <-  c("p-value")
temp[1:10]
summary(temp)
## Finalmene se seleciconan los genes más diferencialmente expresados.
temp =topTable(data.fit.eb, adjust = "fdr")
results <- decideTests(data.fit.eb)


#########################################################################################################

########################################################################################################
#LIMMA GRUPO4
###########
f = factor(groups,levels=c("normal","stage I cRCC","stage II cRCC"))
## Se crea el modelo necesario para la implementaci�n del m�toMCF7Areprodo
design = model.matrix(~ 0 + f)
colnames(design) = c("normal","stageIcRCC","stageIIcRCC")
colnames(grupo4) <- c(groups)

data.fit = lmFit(grupo4,design)
data.fit$coefficients[1:10,]

#####

## Se crea un segundo factor de prueba 
f2 = factor(groups,levels=c("normal","stageIcRCC","stageIIcRCC"))
design2 = model.matrix(~ 0 + f2) # y un segundo modelo de prueba
colnames(design2) = c("normal","stageIcRCC","stageIIcRCC")
#Se crea la matrix de contraste
contrast.matrix = makeContrasts(normal - stageIcRCC - stageIIcRCC,levels=design2)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
#Mediante el metodo eBayes se obtienen los genes diferencialmente expresados
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients[1:10,]

data.fit.eb$t[1:10,]


## Se trabajan loS datos para poder mostrarlos más ordenandamente
data.fit.eb$p.value[1:10,]
temp = data.matrix(data.fit.eb$p.value)
temp = temp[order(temp)]
temp[1:10]
temp = data.matrix(temp)
colnames(temp) <-  c("p-value")
temp[1:10]
summary(temp)
## Finalmene se seleciconan los genes más diferencialmente expresados.
temp =topTable(data.fit.eb, adjust = "fdr")
results <- decideTests(data.fit.eb)


#########################################################################################################



##############################################################################################3
#### ARREGLOS SOLICITADOS POR EL PROFESOR
## y se aplica un algortimo de klustering para el análisis de los datos.
algoritmoKgnuevo <- Kmeans(x = matriz_expresion_discriminativoskluster, centers=4, method="kendall",iter.max = 100)
par(mfrow=c(1,2))
grafico2<-plot(matriz_expresion_discriminativoskluster, col = algoritmoKgnuevo$cluster,
              type='n', main="K-means en genes")
points(algoritmoKgnuevo$centers, col = c("green","blue"), pch = 15, cex = 1)
text(matriz_expresion_discriminativoskluster, labels=rownames(matriz_expresion_discriminativoskluster),
     col=algoritmoKgnuevo$cluster)
