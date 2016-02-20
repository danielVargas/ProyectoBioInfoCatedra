#source("http://www.bioconductor.org/biocLite.R")
#biocLite("ISLR")
######################
#biocLite("GEOquery")
#biocLite("Biobase")
#biocLite("plyr")
#biocLite("siggenes")
#biocLite("limma")
# http://www.johnverostek.com/wp-content/uploads/2014/02/ISLR-First-Printing_CH10-Clustering.pdf
library("ISLR")
#################
library(Biobase)
library(GEOquery)
library(plyr)
library(siggenes)
library(limma)
library(amap)
############
a = data.matrix(NCI60$data)
b = data.matrix(NCI60$labs)

rownames(a)<-b
renal = data.matrix(a[rownames(a)=="RENAL", ])

######### PCA

pr.out =prcomp (a , scale = TRUE)

Cols= function (vec ){
  cols= rainbow (length ( unique (vec )))
  return (cols[as.numeric (as.factor ( vec ))]) 
}

par ( mfrow =c(1 ,2) )
plot(pr.out$x [ ,1:2] , col =Cols(b), pch =19,
       xlab ="Z1", ylab =" Z2")
plot(pr.out$x [,c(1 ,3) ], col =Cols(b), pch =19,
     xlab ="Z1", ylab =" Z3")

summary(pr.out)


plot(pr.out)

####################


pve =100* pr.out$sdev ^2/ sum (pr.out$sdev ^2)
par ( mfrow =c(1 ,2) )
plot(pve , type ="o", ylab =" PVE ", xlab =" Principal Component ",
       col =" blue ")
plot( cumsum (pve ), type ="o", ylab =" Cumulative PVE ", xlab ="
        Principal Component ", col =" brown3 ")


############## CLUSTERING


sd.data= scale (a)
par ( mfrow =c(1 ,3) )
data.dist =dist(sd.data)
plot( hclust ( data.dist ), labels =b , main =" Complete
        Linkage ", xlab ="", sub ="", ylab ="")
plot( hclust ( data.dist , method ="average"), labels =b,
        main =" Average Linkage ", xlab ="", sub ="", ylab ="")
plot( hclust ( data.dist , method ="single") , labels =b ,
        main =" Single Linkage ", xlab ="", sub ="", ylab ="")



hc.out =hclust (dist (sd.data))
hc.clusters = cutree (hc.out ,4)
table (hc.clusters ,b)

par (mfrow =c(1 ,1) )
plot(hc.out , labels =b)
abline (h=139 , col =" red ")

####################################

# CARGADO DE DATOS PARA LOS NOMBRES DE LAS MUESTRAS DESDE EL GEN OMNIBUS CON BASE DE DATOS DE AHI
db <- getGEO('GDS1761', destdir=".")

# CARGADO DE METADADOS DEL DATASET
Meta(db)$dataset_id
Meta(db)$platform_organism
Meta(db)$feature_count
Meta(db)$sample_count

colnames(Table(db))


Table(db)[1:5,1:7]

matriz_log<- GDS2eSet(db, do.log2=TRUE)

NCI60.genes=featureData(matriz_log)$"Gene symbol"
###############################################################
### INICIO DE BÚSQUEDA DE GENES DIFERENCIALMENTE EXPRESADOS
#########################################################3
NEC60.matrizexpresion = t(sd.data)

###########################################################################
### SE ELIMINAN VALORES NA Y ESPACIOS VACIOS
NEC60.matrizexpresionomit <- na.omit(NEC60.matrizexpresion) #NA
NEC60.matrizexpresionomit <- NEC60.matrizexpresion[rownames(NEC60.matrizexpresionomit) != "",] #vacios
NEC60.matrizexpresionomit[1:10,1:6]
###############################
rownames(NEC60.matrizexpresionomit)<-gsub("/","-",rownames(NEC60.matrizexpresionomit))
#promediar genes repetidos
NEC60.matriz_expresion_promedio<-aggregate(NEC60.matrizexpresionomit, by=list(rownames(NEC60.matrizexpresionomit)),
                                     FUN=mean)
#Al promediar los nombres de genes se convierten en columnas
#por eso se deben reasignar como nombres de filas
NEC60.matriz_expresion_promedio<- data.frame(NEC60.matriz_expresion_promedio[,-1],
                                       row.names=NEC60.matriz_expresion_promedio[,1])

NEC60.matriz_expresion_promedio[1:10,1:6]
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
pvals=apply(NEC60.matriz_expresion_promedio,1,function(x) {t.test(x[14:21],x[1:13])$p.value})
#Asignar p.value al gen correspondiente
NEC60.matriz_expresion_promedio$p <- pvals
#Ordenar de - a + la matriz de acuerdo al p.value de cada gen.
ordenarValoresP<-NEC60.matriz_expresion_promedio[order(pvals),]
NEC60.matriz_expresion_discriminativos<-subset(ordenarValoresP, p<=0.05)
head(NEC60.matriz_expresion_discriminativos)
# SE NOMRAN COLUMNAS NUEVAMENTE AHORA CON "P.VALUE"
colnames(NEC60.matriz_expresion_discriminativos) <- c(b,"p.value")
########################################################################################################
#LIMMA
###########
groups = b
f = factor(groups,levels=c("RENAL","CNS","BREAST","NSCLC","UNKNOWN","OVARIAN","MELANOMA","PROSTATE","LEUKEMIA","K562B-repro","K562A-repro","COLON","MCF7A-repro","MCF7D-repro"))
## Se crea el modelo necesario para la implementación del métoMCF7Areprodo
design = model.matrix(~ 0 + f)
colnames(design) =c("RENAL","CNS","BREAST","NSCLC","UNKNOWN","OVARIAN","MELANOMA","PROSTATE","LEUKEMIA","K562Brepro","K562Arepro","COLON","MCF7Arepro","MCF7Drepro")
NEC60.matriz_expresion_discriminativos_sinp<-NEC60.matriz_expresion_discriminativos[ ,-65]
colnames(NEC60.matriz_expresion_discriminativos_sinp) <- c(groups)

data.fit = lmFit(NEC60.matriz_expresion_discriminativos_sinp,design)
data.fit$coefficients[1:10,]

#####

## Se crea un segundo factor de prueba 
f2 = factor(groups,c("RENAL","CNS","BREAST","NSCLC","UNKNOWN","OVARIAN","MELANOMA","PROSTATE","LEUKEMIA","K562Brepro","K562Arepro","COLON","MCF7Arepro","MCF7Drepro"))
design2 = model.matrix(~ 0 + f2) # y un segundo modelo de prueba
colnames(design2) = c("RENAL","CNS","BREAST","NSCLC","UNKNOWN","OVARIAN","MELANOMA","PROSTATE","LEUKEMIA","K562Brepro","K562Arepro","COLON","MCF7Arepro","MCF7Drepro")
#Se crea la matrix de contraste
contrast.matrix = makeContrasts(RENAL - CNS - BREAST - NSCLC - UNKNOWN- OVARIAN- MELANOMA -PROSTATE- LEUKEMIA -K562Brepro -K562Arepro -COLON-MCF7Arepro-MCF7Drepro,levels=design2)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
#Mediante el metodo eBayes se obtienen los genes diferencialmente expresados
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

data.fit.eb$coefficients[1:10,]

data.fit.eb$t[1:10,]


## Se trabajan loS datos para poder mostrarlos mÃ¡s ordenandamente
data.fit.eb$p.value[1:10,]
temp = data.matrix(data.fit.eb$p.value)
temp = temp[order(temp)]
temp[1:10]
temp = data.matrix(temp)
colnames(temp) <-  c("p-value")
temp[1:10]
summary(temp)
## Finalmene se seleciconan los genes mÃ¡s diferencialmente expresados.
temp =topTable(data.fit.eb, adjust = "fdr")
results <- decideTests(data.fit.eb)


## Y se arma una matriz para ingrearla al heatmap
matrix.genes.dif = c()
for (variable in rownames(temp)) {
  matrix.genes.dif = rbind(matrix.genes.dif, NEC60.matriz_expresion_promedio[variable,])
  #print(NEC60.matriz_expresion_promedio[variable,])
}
matrix.genes.dif<-na.omit(matrix.genes.dif[,1:21])

t = as.matrix(matrix.genes.dif)

heatmap(t, main = "Heatmap de genes")