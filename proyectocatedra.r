#source("http://www.bioconductor.org/biocLite.R")
#biocLite("ISLR")
# http://www.johnverostek.com/wp-content/uploads/2014/02/ISLR-First-Printing_CH10-Clustering.pdf
library("ISLR")

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


