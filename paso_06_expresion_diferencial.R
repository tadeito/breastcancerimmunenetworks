#########################################
### Analisis de expresion diferencial ###
#########################################

# Analisis de exp. dif. para comparar los patrones entre muestras normales y los distintos subtipos de cancer de mama
library(limma)

BiocManager::install("KEGGprofile")

# Definir la ubicacion de los archivos
list.files("data")

files <- c("data/Healthy_PBCMC_ARSYM.txt","data/LumA_PBCMC_ARSYM.txt","data/LumB_PBCMC_ARSYM.txt","data/Her2_PBCMC_ARSYM.txt","data/Basal_PBCMC_ARSYM.txt")

matrices <- lapply(files, read.table, sep="\t",header=TRUE, row.names = 1)
#str(matrices) # Verificar que se hayan cargado correctamente
#lapply(matrices,dim) # verificar que las dimensiones sean las esperadas

ordenar <- function(x){
x <-  x[ order(rownames(x)),] 
}

matrices <- lapply(matrices,ordenar) # Asegurarse de que los rownames correspondan en todas las matrices
# all.equal(rownames(matrices[[1]]),rownames(matrices[[3]]))
# [1] TRUE

matrices <- lapply(matrices,as.matrix)
#lapply(matrices,class)
names(matrices) <- c("Normal","LumA","LumB","Her2","Basal") # Añadir etiquetas a la lista de matrices de expresión.

exp_mat <- cbind(matrices[[1]],matrices[[2]],matrices[[3]],matrices[[4]],matrices[[5]])
# class(exp_mat)
# dim(exp_mat)

## Matriz de diseño ##

# La matriz de diseño decribe los fenotipos a los que pertenece cada muestra en la matriz de expresión. Cada columna es un fenotipo y cada renglón es una muestra. 

design <- matrix(rep(0,5*536), nrow=536)
colnames(design) <- c("Normal","LumA","LumB","Her2","Basal")

design[1:101,1]<- 1
design[102:264,2]<- 1
design[265:322,3]<- 1
design[323:394,4]<- 1
design[395:536,5]<- 1

## Matriz de contrastes ##

#La matriz de contrastes describe que comparaciones se van a hacer entre los fenotipos definidos en la matriz de diseño
cont.matrix = makeContrasts("LumA - Normal","LumB - Normal","Her2 - Normal","Basal - Normal", levels=design)

## Ajuste de modelos y cálculo de los coeficientes de expresión diferencial ##
fit = lmFit(exp_mat, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

save(fit2,file = "results/paso_06_expresion_diferencial/lmfit_ebayes_subtipos.Rdata")

# Para recuperar los coeficientes usamos la función topTable. podemos indicar que nos devuelva todos los coeficientes de la comparación con "n= "

topTable(fit2, coef=1, adjust = "fdr",n=15267)

# sacar las tablas de lfc para todos los subtipos

tablas_expdif <-  list(LumA=topTable(fit2, coef=1, adjust = "fdr",n=15267),LumB=topTable(fit2, coef=2, adjust = "fdr",n=15267),Her2=topTable(fit2, coef=3, adjust = "fdr",n=15267),Basal=topTable(fit2, coef=4, adjust = "fdr",n=15267))

# Hacer un atributo de color para cada tabla

i=1

for(i in 1:length(tablas_expdif)){
DE_status <- ifelse( tablas_expdif[[i]]$logFC >1 & tablas_expdif[[i]]$B>6, "up",
        ifelse(tablas_expdif[[i]]$logFC < (-1) & tablas_expdif[[i]]$B>6, "down"
          
          , "no_change")
        )
tablas_expdif[[i]]$DE_status <- DE_status

write.table(tablas_expdif[[i]],file = paste0("results/paso_06_expresion_diferencial/DE_",names(tablas_expdif)[i],"_vs_normal.txt"),sep = "\t",quote = FALSE, col.names = TRUE,row.names = TRUE) # Añadir al header de forma externa la etiqueta de la columna de genes "genesym"
}


########################################################################
### Comparación de los patrones globales de expresión entre subtipos ###
########################################################################


# correlacion de spearman de los patrones de expresión diferencial entre subtipos

#all.equal(rownames(tablas_expdif[[1]]),rownames(tablas_expdif[[2]]))
#[1] "15256 string mismatches"
# Como los datos estan ordenados de forma distinta, ordenamos los datos antes de unirlos.
i=1

for(i in 1:length(tablas_expdif)){
tablas_expdif[[i]] <- tablas_expdif[[i]][order(rownames(tablas_expdif[[i]])),]
}

#all.equal(rownames(tablas_expdif[[1]]),rownames(tablas_expdif[[2]]))
#[1] TRUE
# Se cumple para todas las demas comparaciones

#LumA vs LumB
cor.test(tablas_expdif[[1]]$logFC,tablas_expdif[[2]]$logFC,method = "spearman")
# Spearman's rank correlation rho
# 
# data:  tablas_expdif[[1]]$logFC and tablas_expdif[[2]]$logFC
# S = 1.3465e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9772969 

#LumA vs Her2
cor.test(tablas_expdif[[1]]$logFC,tablas_expdif[[3]]$logFC,method = "spearman")
# Spearman's rank correlation rho
# 
# data:  tablas_expdif[[1]]$logFC and tablas_expdif[[3]]$logFC
# S = 5.5646e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9061739 


#LumA vs Basal
cor.test(tablas_expdif[[1]]$logFC,tablas_expdif[[4]]$logFC,method = "spearman")
# Spearman's rank correlation rho
# 
# data:  tablas_expdif[[1]]$logFC and tablas_expdif[[4]]$logFC
# S = 2.7777e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.5316459 

#LumB vs Her2
cor.test(tablas_expdif[[2]]$logFC,tablas_expdif[[3]]$logFC,method = "spearman")
# Spearman's rank correlation rho
# 
# data:  tablas_expdif[[2]]$logFC and tablas_expdif[[3]]$logFC
# S = 1.6101e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9728525 


#Lumb vs Basal
cor.test(tablas_expdif[[2]]$logFC,tablas_expdif[[4]]$logFC,method = "spearman")
# Spearman's rank correlation rho
# 
# data:  tablas_expdif[[2]]$logFC and tablas_expdif[[4]]$logFC
# S = 1.9482e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.6715129 


#Her2 vs Basal
cor.test(tablas_expdif[[3]]$logFC,tablas_expdif[[4]]$logFC,method = "spearman")
# Spearman's rank correlation rho
# 
# data:  tablas_expdif[[3]]$logFC and tablas_expdif[[4]]$logFC
# S = 1.1571e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8048907 



# Algunas pruebas sobre la distribucion de datos y las correlaciones
which(rownames(matrices[[1]])=="CD4")
hist(matrices[[5]][2197,],breaks = 20)

which(rownames(matrices[[1]])=="CD8A")
which(rownames(matrices[[1]])=="CD8B")

hist(matrices[[1]][2221,],breaks = 100)
plot(x=matrices[[1]][2221,],y=matrices[[1]][2222,])





