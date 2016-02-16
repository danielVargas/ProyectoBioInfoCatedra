#source("http://www.bioconductor.org/biocLite.R")
#biocLite("ISLR")

library("ISLR")

a = data.matrix(NCI60$data)
b = data.matrix(NCI60$labs)

rownames(a)<-b

renal = data.matrix(a[rownames(a)=="RENAL", ])