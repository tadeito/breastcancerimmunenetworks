##################################
### Atributos para cytoscape  ####
##################################

# Generar las tablas de atributos para cytoscape a partir de los datos del clustering con infomap
# los atributos sirven para hacer el mapeo visual de los clusters

# Este proyecto está alojado en:
# > getwd()
# [1] "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos"

# La informacion de los clusters está en el directorio /Results/paso_01_Community_detection/
# La informacion de los enriquecimientos está en el directorio /Results/paso_02_Module_enrichments


library(igraph)
library(randomcoloR) # para generar una paleta de colores para los modulos

# obtenemos los datos de los clusters
load("Results/paso_01_Community_detection/clusters_top100k")
clusters_infomap <- clusters_top100k
#rm(clusters_top100k)

# Obtenemos los modulos enriquecidos en cada una de las redes
directorios <- list.files("Results/paso_02_Module_enrichments/",full.names = TRUE)

ids_modulos_enriquecidos <- list()

for(d in 1:length(directorios)){
enriched_modules <- list.files(directorios[d])
ids_modulos_enriquecidos[[d]] <- as.integer(gsub("_.+txt$","",gsub(".*mod_","",enriched_modules)))
names(ids_modulos_enriquecidos)[d] <- names(clusters_infomap)[d]
}

# Generamos las tablas de atributos y las listas de genes en modulos de diez o mas genes. estas se guardan en el directorio /data/Cytoscape y se usarán en las visualizaciones

# n=1 # Testvalue por red
# bucle para repetir por todas las redes de cancer de mama
for(n in 1:length(clusters_infomap)){ # Abre loop por redes 
modulos <- membership(clusters_infomap[[n]])
modulos <-as.data.frame(cbind(gene_id=names(modulos),modulo=modulos),stringsAsFactors=FALSE)
# columna para clasificar los modulos en enriquecidos (1) y no enriquecidos (0)
modulos$enriquecido <- ifelse(modulos$modulo %in% ids_modulos_enriquecidos[[n]],1,0)

#paleta de colores
nmods_10 <- length(table(modulos$modulo)[table(modulos$modulo)>=10])
palette <- data.frame (modulo = table(modulos$modulo)[table(modulos$modulo)>=10],color =distinctColorPalette(nmods_10))

mod_colors <- as.character(palette[match(modulos$modulo,palette$modulo.Var1),3])
mod_colors[is.na(mod_colors)] <- "#FFF"

modulos$color <- mod_colors
# Lista de genes pertenecientes a los modulos con 10 o mas genes
modulos <- modulos[modulos[,2]%in%names(table(modulos[,2]))[table(modulos[,2]) >= 10],]

mayores_a_10 <- modulos$gene_id

# Escribir resultados en disco
write.table(modulos,file=paste0("Results/paso_04_visualizaciones_Cytoscape/cytoscape_top100k_",names(clusters_infomap)[n],".txt"),row.names = FALSE,col.names = TRUE,sep = "\t", quote = FALSE)
write.table(mayores_a_10,file = paste0("Results/paso_04_visualizaciones_Cytoscape/Mayores_a_10_",names(clusters_infomap)[n],".txt"),row.names = FALSE,col.names = FALSE, quote = FALSE)
} # cierra loop por redes

