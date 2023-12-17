##########################################################################
### Module enrichments for genome-wide networks breast cancer subtypes ###
##########################################################################

# Author: Tadeo  2019/05/12 
# this script is intended to make enrichments along multiple modules in MI networks.
# Modules are defined as lists of names taken from igraph community detection outuput. 
# La entrada para este script es una lista donde cada objeto 


# R version 3.5.1 (2018-07-02)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("HTSanalyzeR")
library("HTSanalyzeR") # Biblioteca para hacer los enriquecimientos
library("GO.db") # Referencia, anotaciones de Gene Ontology
library("org.Hs.eg.db")
library("data.table") # Para convertir listas a data.frames

########
## 01 ## Cargar los clusters (resultados del paso 01)
########

load("Results/paso_01_Community_detection/clusters_top100k")

########
## 02 ## Preparar los objetos para trabajar
########

i=1
clusters <- clusters_top100k

GO_BP_1 <- GOGeneSets(species="Hs", ontologies=c("BP")) # obtenemos las listas de genes para las categorias de GO 
GOTERM_list <- as.list(GOTERM) # Usamos la base de datos GOTERM para obtener la informacion de los GOIDs

# Funcion para extraer la informacion asociada a los GO IDs en forma de data.frame
gotable <- function(x){
  data.frame(GOID = GOID(GOTERM_list[[x]]),
             TERM = Term(GOTERM_list[[x]]))#,
  #Synonym(GOTERM_list[[x]])),
  #Secondary(GOTERM_list[[x]])),
  #Definition(GOTERM_list[[x]])),
  #Ontology(GOTERM_list[[x]]))
}

# Esta parte es muy engorrosa y solo sirve para que los nombres de los IDs se presenten en forma de character y puedan ser filtrados de forma adecuada
goIDmap <-lapply(1:length(GOTERM_list),gotable )
GOIDmap <- as.data.frame(rbindlist(goIDmap))
GOIDmap <- data.frame(ID=as.character(GOIDmap[,1]),Term = as.character(GOIDmap[,2]),stringsAsFactors = FALSE)


# definimos el universo de nombres como los genes totales de la matriz de expresión de TCGA

#cargamos la lista de nombres
lista_todos_genes <- read.table("~/Tadeo/Doctorado_Inflamacion_BC_subtipos/Data/lista_todos_genes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

GeneIDuniverse <-select(org.Hs.eg.db, keys = lista_todos_genes$V1,columns = c("ENTREZID"),keytype = "SYMBOL")
# Eliminamos los mapeos multiples
GeneIDuniverse <- GeneIDuniverse[!duplicated(GeneIDuniverse$SYMBOL) & !duplicated(GeneIDuniverse$ENTREZID),]

universe <- as.character(GeneIDuniverse[!is.na(GeneIDuniverse$ENTREZID),2]) 

# Para asegurarnos que los enriquecimientos no tengan un sesgo debido a genes que no se encuentran en la prueba, filtramos las listas de genes en procesos dejando aquellos genes que se encuentran en la matriz de expresión.

filtro <- function(g,u){
r <-g[g%in%u]
return(r)
}

GO_BP <- lapply(GO_BP_1, filtro, universe)


# Definimos cual de las redes queremos analizar mas adelante
# names(clusters)
# [1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"
# selected_network = 1 # Valor de prueba para escoger los clusters de una de las redes

for(selected_network in 1:length(clusters)) { # Abrimos loop por redes (cortes)
  


# Creamos una carpeta para guardar los resultados de cada red:
system2("mkdir",args = paste0("Results/paso_02_Module_enrichments/cluster_enrichments_",names(clusters)[selected_network]))

Modulos_infomap <-cbind(cluster=membership(clusters[[selected_network]]),symbol= names(membership(clusters[[selected_network]]))) # Obtenemos una tabla donde cada renglon corresponde al nombre de un gen y el numero de modulo al que pertenece
  
Modulos_infomap <- data.frame(cluster= Modulos_infomap[,1],symbol=as.character(Modulos_infomap[,2]),stringsAsFactors = FALSE)
#Modulos_infomap$symbol <- as.character(Modulos_infomap[,2]) # por que a R le vale madre y le gusta que todo sean factores, insistimos en que queremos caracteres


modulos_lista <- unique(Modulos_infomap$cluster) # Obtenemos la lista de módulos unicos y rdenados
# Para repetir el analisis con otra red es necesario especificar el numero dentro de la lista






########
## 03 ## Correr el enriquecimiento
########

# Enriquecer todas las comunidades de la red que tengan 20 genes o mas y guardar el resultado si existe enriquecimiento significativo despues de corregir por FDR
#m=1 # Test value
#t=1 # Test value

for(m in 1:length(modulos_lista)){ # abrimos loop por modulos
  genelist <- Modulos_infomap[Modulos_infomap[,1] == modulos_lista[m],2] # Extraemos la lista de genes asociada a un modulo particular a partir de la tabla
      modID <- genelist[1] # Guardamos el nombre del primer gen del modulo (esto facilita identificar los modulos)
  
  if(sum(genelist %in% GeneIDuniverse$SYMBOL)==0) next # Esta instrucción es necesaria para evitar errores si ninguno de los genes de la lista se encuentra presente en la base de datos
  
  mappedID <- select(org.Hs.eg.db, keys = genelist ,columns = c("ENTREZID"),keytype = "SYMBOL") # Usamos la ista de genes para obtener los EntrezID 
  IDs <- mappedID[!is.na(mappedID$ENTREZID),2] # Eliminamos aquellos que no mapearon y aparecen como N/A
  
  if(length(IDs)<10) next # Si una comunidad tiene menos de 10 genes, pasamos a la siguiente
  
  # Hora de correr el enriquecimiento!
  ORA_comunidad <- multiHyperGeoTest(collectionOfGeneSets = GO_BP, universe = universe, hits =  IDs, pAdjustMethod = "BH", minGeneSetSize = 5 )
  
  # selecionamos los enriquecimientos significativos despues de la corrección del p-value por FDR
  significativos <- as.data.frame(ORA_comunidad)
  significativos$Adjusted.Pvalue[which(is.na(significativos$Adjusted.Pvalue))]<-1 # este es un filtro para eliminar los p-values que se registran como NA
  significativos <- significativos[significativos$Adjusted.Pvalue <= 0.05,]
  
  if(length(rownames(significativos))==0) next # Si un enriquecimiento NO tiene ninguna categoria significativa, pasamos al siguiente modulo        
  # Añadimos la informacion de los terminos asociados a cada GOID
  terms<-row.names(significativos)
  TERM_name <- data.frame(TERM=character(),stringsAsFactors = FALSE) # Iniciamos un data.frame para almacenar los nombres de los terminos

   
  
  # Almacenamos el Term de cada GOID encontrado
  for(t in 1:length(terms)){
    TERM_name[t,1] <-GOIDmap[GOIDmap$ID==terms[t],2]
      }

 
  # Anexamos la columna de Term a la tabla de enriquecimiento
      enrichment_results <- data.frame(GOid=row.names(significativos),TERM_name,significativos)

      enrichment_results <- enrichment_results[order(enrichment_results$Adjusted.Pvalue, decreasing = FALSE),] # ordenamos empezando por los p-values mas significativos (de menor a mayor) 
  
  write.table(enrichment_results,file = paste0("Results/paso_02_Module_enrichments/cluster_enrichments_",names(clusters)[selected_network],"/mod_",modulos_lista[m],"_",modID,".txt"), row.names = FALSE,col.names = TRUE,sep = "\t",quote = FALSE) 

} # Cerramos loop por modulos

} # Cerramos loop por redes (cortes)
# La salida consiste en un archivo por cada comunidad que resulto enriquecida y debe verse así:

# > head(enrichment_results)
#   GO          TERM                                                Universe.Size  Gene.Set.Size  Total.Hits  Expected.Hits  Observed.Hits  Pvalue        Adjusted.Pvalue
# 1 GO:0006357  regulation of transcription by RNA polymerase II    8620           633            20          1.468677       18             5.078495e-19  2.854114e-16




