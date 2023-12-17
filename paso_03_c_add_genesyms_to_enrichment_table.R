######################################################################
### Recuperar los genes que hacen los enriquecimientos en la tabla ###
######################################################################

# Utilizando la información de la tabla de enriquecimientos, recuperar los genes que lo producen.
# Tomamos la tabla de enriquecimientos inmunes de la carpeta resultados 03
tabla_inmunes <-read.table(file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/modulos_procesos_inflamacion.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

# Tambien requerimos las listas de genes de GO, generadas en el script del paso anterior:
GO_BP

# Y las listas de genes de los módulos que se generaron en el paso 01
memberships <- lapply(clusters,membership)

# Así como la correspondencia de genesymbol a ensembl
GeneIDuniverse



# obtener numero de modulo desde la tabla
#n=1 # contador para entrada en la tabla

GeneSym <- character()
ENTREZID <- character()

for(n in 1:dim(tabla_inmunes)[1]){
# encontrar los genes en el modulo indicado en la tabla y convertirlos a ENSEMBLID
mod_no <- as.integer(strsplit(tabla_inmunes$module[n], split = "_")[[1]][2])
pheno_no <- which(names(memberships)==tabla_inmunes$phenotype[n])
mod_genesymbols <- names(memberships[[pheno_no]][memberships[[pheno_no]]==mod_no])
mapped_names <- GeneIDuniverse$ENTREZID[GeneIDuniverse$SYMBOL%in%mod_genesymbols]

# Obtener la intersección entre los genes del módulo y los del proceso
intersect_mod_process <- mapped_names[mapped_names%in%GO_BP[[which(names(GO_BP)%in%tabla_inmunes$GOid[n])]]]

#obtener de vuelta los genesybol
results <- GeneIDuniverse$SYMBOL[GeneIDuniverse$ENTREZID%in%intersect_mod_process]
results2 <- paste(results,collapse=" ")

GeneSym[n]<-results2
ENTREZID[n]<-paste(intersect_mod_process,collapse=" ")
}

tabla_inmunes$GeneSym <- GeneSym
tabla_inmunes$ENTREZID <- ENTREZID

write.table(tabla_inmunes,file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/modulos_procesos_inflamacion.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

###############################################
## Listas de genes de los modulos de interés ##
###############################################

# Extraemos las listas de genes de los modulos de interes para facilitar la consulta y hacer comparaciones rapidas. 
n=1

Fenotipo <- character()
Modulo <- character()

# obtenemos la lista de módulos por fenotipo
for(n in 1:length(tabla_inmunes$phenotype)){
Modulo[n] <- as.integer(strsplit(tabla_inmunes$module[n], split = "_")[[1]][2])
Fenotipo[n] <- which(names(memberships)==tabla_inmunes$phenotype[n])
}

Listas_Genes <- data.frame(Fenotipo=Fenotipo,Modulo=Modulo,stringsAsFactors = FALSE)

# simplificamos la lista
Listas_Genes <- unique(Listas_Genes)

Genes <- character()

# obtenemos los genes de cada módulo
n=1

for(n in 1:length(Listas_Genes$Modulo)){
pheno_no <- as.integer(Listas_Genes$Fenotipo[n])
mod_no <- as.integer(Listas_Genes$Modulo[n])
mod_genesymbols <- names(memberships[[pheno_no]][memberships[[pheno_no]]==mod_no])
Genes[n] <- paste(mod_genesymbols,collapse = " ")
Listas_Genes$Fenotipo[n] <- names(memberships)[pheno_no]
}

Listas_Genes$Genes <- Genes


write.table(Listas_Genes,file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Listas_genes_modulos_inmunes.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

Listas_Genes
