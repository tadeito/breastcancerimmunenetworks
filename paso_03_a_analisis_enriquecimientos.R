########################################
### Analisis de los enriquecimientos ###
########################################

library(gplots)
library(RColorBrewer)
library(dplyr)
##########################################################################
## Heatmaps para comparar los enrichments entre modulos de una sola red ##
##########################################################################

directorios <- list.files("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments/",full.names = TRUE)
nnames <- gsub("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments//cluster_enrichments_","",directorios)

# obtenemos los enriquecmientos
#d=1


# Este bucle extrae la matriz de enriquecimientos por cada red y hace un heatmap con todos los enriquecimientos al final de cada vuelta. Los heatmapas se guardan como archivos .pdf

for(d in 1:length(directorios)) { # Abre loop todas las redes
enrichment_files <- list.files(directorios[d],full.names = TRUE)
enrichment_files_names <- list.files(directorios[d])
enrichments_todos <- lapply(enrichment_files, read.delim, header=TRUE)
names(enrichments_todos) <- enrichment_files_names

# Extraemos el universo de procesos enriquecidos

i=2

# Vector para almacenar los resultados
universo <- character()

# extraemos todos los terms de todos los modulos 
for(i in 1:length(enrichments_todos)){
terms <-as.character(enrichments_todos[[i]][,2])
universo <- c(universo,terms)
}

# simplificamos la lista de terminos
universo <- unique(universo)

# Matriz de enrichments vacía para almacenra los resultados
enrichment_matrix <- matrix(nrow = length(universo),ncol = length(enrichments_todos))

# Identificamos cuales terms estan en cada mdulo y hacemos un vector con su respectivo p-val ajustado

for(m in 1:length(enrichments_todos)){
p <- numeric()
for (t in 1:length(universo)){
p[t] <- ifelse(universo[t] %in% as.character(enrichments_todos[[m]][,2]),
       -log10(enrichments_todos[[m]][which(universo[t]==enrichments_todos[[m]][,2]),9]),

              0)
}

enrichment_matrix[,m] <- p
}


row.names(enrichment_matrix) <- universo
colnames(enrichment_matrix)<- names(enrichments_todos)

# Heatmap para los resultados 



my_palette <- colorRampPalette(c("white","blue","orange","red"))(n = 40)


pdf(file = paste0("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Heatmap_procesos_",nnames[d],".pdf"))
heatmap.2( x = enrichment_matrix ,
           main = paste("Enriched modules", nnames[d]),   # heat map title
           density.info = "none",  # turns off density plot inside color legend
           trace = "none",         # turns off trace lines inside the heat map
           margins = c(15,2),     # widens margins around plot
           col = my_palette,       # use on color palette defined earlier
           labRow = "",            # change row and column labels
           #labCol = "",
           ylab = "GO Terms",
           xlab = "Module",
           #Rowv = as.dendrogram(row.cluster), # apply selected clustering method
           #Colv = as.dendrogram(col.cluster), # apply selected clustering method
           keysize = 0.8,
           key.title = "Enrichment",
           key.xlab = " -log10 adj p-val"
           #,        # size of color key
           #ColSideColors = colsidecols ,dendrogram = "column",#Additional Options ## Color labeled columns
           #scale = "column" # cambiar por "column" o "row" segun se quiera la escala de colores normalizada en una direccion particular
           
)

dev.off()

} # cierra loop todas las redes

# Esta secion se podria abordar utilizando un formato tidy como en el caso de la ultima seccion de este script ****


######################################
### Preparar las tablas para latex ###  ### Arreglar esta seción ######
######################################
directorios <- list.files("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments/",full.names = TRUE)
nnames <- gsub("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments//cluster_enrichments_","",directorios)


clusters_inflamacion <- list(
  # hay que mencionar, no el nombre del módulo, si no la posición del archivo que lo contiene en la lista de archivos  
  Basal = c(34,38,44,2,4,7,9,11,17,19,25,26,31),
  Her2 = c(43,45,48,4,5,6,13,14,18,20,26,30,38),
  LumA = c(36,41,42,44,45,47,3,4,5,10,11,16,23,27,31),
  LumB = c(34,35,36,38,5,7,10,11,18,19,27,28),
  Normal= c(2,6,19,22,37,38,40,44,45,48,25,34)
)


# Genera las tablas para latex con las columnas de interes y  para los modulos seleccionados
# d=1
for(d in 1:length(directorios)){

enrichment_files <- list.files(directorios[d],full.names = TRUE)
enrichment_files_names <- list.files(directorios[d])
enrichments_todos <- lapply(enrichment_files, read.delim, header=TRUE)
names(enrichments_todos) <- enrichment_files_names


GO_table <- data.frame(GOid=character(),TERM=character(),Adjusted.Pvalue=numeric())
m <- data.frame(GOid=character(),TERM=character(),Adjusted.Pvalue=numeric())

#i=34
for (i in clusters_inflamacion[[d]]){
e <- enrichments_todos[[i]][,c(1,2,9)]
m[1,] <- data.frame(GOid=" ",TERM=paste(enrichment_files_names[i]),Adjusted.Pvalue=0) 
GO_table <- rbind(GO_table,m,e)

  } # cierra loop modulos inflamacion
GO_table$Adjusted.Pvalue <-format(GO_table$Adjusted.Pvalue,digits = 3)

write.table(GO_table,file = paste0("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/tabla_Latex_",nnames[d],"_solo_SI.txt"),sep="&",col.names = TRUE,row.names = FALSE,quote = FALSE)
} # cierra loop directorios


##################################################################################
### Comparación de enrichments procesos inmunologicos entre diferentes redes  ####
##################################################################################

directorios <- list.files("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments/",full.names = TRUE)
nnames <- gsub("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments//cluster_enrichments_","",directorios)



clusters_inflamacion <- list(
  # hay que mencionar, no el nombre del módulo, si no la posición del archivo que lo contiene en la lista de archivos  
  Basal = c(34,38,44,2,4,7,9,11,17,19,25,26,31),
  Her2 = c(43,45,48,4,5,6,13,14,18,20,26,30,38),
  LumA = c(36,41,42,44,45,47,3,4,5,10,11,16,23,27,31),
  LumB = c(34,35,36,38,5,7,10,11,18,19,27,28),
  Normal= c(2,6,19,22,37,38,40,44,45,48,25,34)
)
d=1


# lista para guardar los enrichments
enrichments_todos <- list()

# loop para leerlos todos
for(d in 1:length(directorios)){
# listamos todos los modulos
  enrichment_files <- list.files(directorios[d],full.names = TRUE)
# filtramos quedandonos con aquellos que tienen procesos relacionados a inflamacion
enrichment_files <- enrichment_files[clusters_inflamacion[[d]]]
enrichment_files_names <- list.files(directorios[d])
enrichment_files_names <- enrichment_files_names[clusters_inflamacion[[d]]]
enrichments_todos[[d]] <- lapply(enrichment_files, read.delim, header=TRUE)
names(enrichments_todos[[d]]) <- enrichment_files_names
}
names(enrichments_todos) <- nnames

# Convertimos la lista en una sola tabla, incluyendo la red y el modulo de procedencia de cada categoria.
tabla_inmunes <- bind_rows(lapply(enrichments_todos , bind_rows, .id = "module" ),.id= "phenotype")

tabla_inmunes$GObp_term <- paste(tabla_inmunes$GOid,tabla_inmunes$TERM)
tabla_inmunes$pheno_module <- paste(tabla_inmunes$phenotype,tabla_inmunes$module)

write.table(tabla_inmunes,file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/modulos_procesos_inflamacion.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

#### En el paso 03_c añadimos los genesymbl que dan origen al enrichment en cada proceso.


# extraer el universo de terms inmunologicos
universo_inmune <- unique(tabla_inmunes$GObp_term)
universo_modulos <- unique(tabla_inmunes$pheno_module)

# Construimos la matriz

# matriz vacía para almacenar los enriquecimientos
matriz_todas_inflamacion <- matrix(nrow = length(universo_inmune),ncol = length(universo_modulos))
colnames(matriz_todas_inflamacion) <- universo_modulos
rownames(matriz_todas_inflamacion) <- universo_inmune

#p=1 #valor por procesos
#m=1 #valor de prueba por modulo

# Llenamos las entradas de la matriz basado en la tabla de enrichments
for(p in 1:length(universo_inmune)){  # abre loop por procesos
for(m in 1:length(universo_modulos)){ # abre loop por modulos
matriz_todas_inflamacion[p,m] <- ifelse(  length(which(universo_inmune[p] == tabla_inmunes$GObp_term & universo_modulos[m] == tabla_inmunes$pheno_module))==0, 
         0,
      -log10(tabla_inmunes$Adjusted.Pvalue[which(universo_inmune[p] == tabla_inmunes$GObp_term & universo_modulos[m] == tabla_inmunes$pheno_module)])
      )
} # cierra loop por modulos
} # Cierra loop por procesos



# Heatmap de los resultados

#Notes from: https://tousu.in/qa/?qa=592906/
#lmat is a matrix describing how the screen is to be broken up. By default, heatmap.2 divides the screen into a four element grid, so lmat is a 2x2 matrix. The number in each element of the matrix describes what order to plot the next four plots in. Heatmap.2 plots its elements in the following order:

# 1.-Heatmap,
# 2.-Row dendrogram,
# 3.-Column dendrogram,
# 4.-Key

lmat2 <- rbind(c(4,0),c(0,3),c(2,1))



my_palette <- colorRampPalette(c("white","blue","orange","red"))(n = 40)

# Vector de colores para la barra de anotación en el heatmap
# Vector de colores para la barra de anotación que aparece arriba del heatmap
col_annot_bar = c( rep("VIOLET", 13), 
                   rep("GREEN", 13),
                   rep("RED", 15),
                   rep("ORANGE", 12),
                   rep("BLUE", 12)
)




pdf(file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Heatmap_modulos_inmune_todas_las_redes.pdf",height = 35,width = 15)
heatmap.2( x = matriz_todas_inflamacion ,
           main = "Modules with immune system enrichments",   # heat map title
           density.info = "none",  # turns off density plot inside color legend
           trace = "none",         # turns off trace lines inside the heat map
           margins = c(13,20),     # widens margins around plot
           col = my_palette,       # use on color palette defined earlier
           #lmat=lmat2,
           lwid = c(2,4),
           lhei = c(1,4),
           #labRow = "                            ",            # change row and column labels
           #labCol = "",
           ylab = "GO Terms",
           xlab = "Module",
           #Rowv = as.dendrogram(row.cluster), # apply selected clustering method
           #Colv = as.dendrogram(col.cluster), # apply selected clustering method
           keysize = 0.1,
           key.title = "Enrichment",
           key.xlab = " -log10 adj p-val",
           colCol = col_annot_bar,
           ColSideColors = col_annot_bar           
           #,        # size of color key
           #ColSideColors = colsidecols ,dendrogram = "column",#Additional Options ## Color labeled columns
           #scale = "column" # cambiar por "column" o "row" segun se quiera la escala de colores normalizada en una direccion particular
           
)

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Basal", "Her2", "LumA","LumB","Normal"), # category labels
       col = c("VIOLET", "GREEN", "RED","ORANGE","BLUE"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)


dev.off()

##############################################################################################
##############################################################################################
##############################################################################################



# Heatmap comparando todos los enriquecimientos de todos los modulos enriquecidos entre todas las redes.


# Obtener la dirección de los archivos
directorios <- list.files("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments/",full.names = TRUE)
nnames <- gsub("/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_02_Module_enrichments//cluster_enrichments_","",directorios)

# obtenemos los enriquecmientos
#d=1

# lista para guardar los enrichments
enrichments_todos <- list()

# loop para leerlos todos
for(d in 1:length(directorios)){
  # listamos todos los modulos
  enrichment_files <- list.files(directorios[d],full.names = TRUE)
  enrichment_files_names <- list.files(directorios[d])
  enrichments_todos[[d]] <- lapply(enrichment_files, read.delim, header=TRUE)
  names(enrichments_todos[[d]]) <- enrichment_files_names
}
names(enrichments_todos) <- nnames

# Convertimos la lista en una sola tabla, incluyendo la red y el modulo de procedencia de cada categoria.
tabla_todos <- bind_rows(lapply(enrichments_todos , bind_rows, .id = "module" ),.id= "phenotype")

tabla_todos$GObp_term <- paste(tabla_todos$GOid,tabla_todos$TERM)
tabla_todos$pheno_module <- paste(tabla_todos$phenotype,tabla_todos$module)

#revisar el resultado
head(tabla_todos)


# extraer el universo de terms inmunologicos para construir la matriz de terminos
universo_terminos <- unique(tabla_todos$GObp_term)
universo_modulos <- unique(tabla_todos$pheno_module)

# Construimos la matriz

# matriz vacía para almacenar los enriquecimientos
matriz_todas <- matrix(nrow = length(universo_terminos),ncol = length(universo_modulos))
colnames(matriz_todas) <- universo_modulos
rownames(matriz_todas) <- universo_terminos

#p=1 #valor por procesos
#m=1 #valor de prueba por modulo

# Llenamos las entradas de la matriz basado en la tabla de enrichments
for(p in 1:length(universo_terminos)){  # abre loop por procesos
  for(m in 1:length(universo_modulos)){ # abre loop por modulos
    matriz_todas[p,m] <- ifelse(  length(which(universo_terminos[p] == tabla_todos$GObp_term & universo_modulos[m] == tabla_todos$pheno_module))==0, 
                                              0,
                                              -log10(tabla_todos$Adjusted.Pvalue[which(universo_terminos[p] == tabla_todos$GObp_term & universo_modulos[m] == tabla_todos$pheno_module)])
    )
  } # cierra loop por modulos
} # Cierra loop por procesos



# Heatmap de los resultados
my_palette <- colorRampPalette(c("white","blue","orange","red"))(n = 40)

#Basal   Her2   LumA   LumB Normal 
#324    248    273    273    345 

# Vector de colores para la barra de anotación que aparece arriba del heatmap
col_annot_bar = c( rep("VIOLET", 44), 
                   rep("GREEN", 49),
                   rep("RED", 48),
                   rep("ORANGE", 39),
                   rep("BLUE", 48)
                   )
# Vector de colores para marcar los nombres de las columnas que aparecen abajo del heatmap

mod_annot_bar = rep("BLACK",228) #rellenamos de "negro" el vector

mod_annot_bar[c(34,38,44,2,4,7,9,11,17,19,25,26,31,87,89,91,48,49,50,57,58,62,64,70,74,82,129,134,135,137,138,140,96,97,98,103,104,109,116,120,124,175,176,177,179,146,148,151,152,159,160,168,169,181,186,199,202,217,218,220,224,225,228,205,213)] <- "RED" # asignamos "rojo" a las posiciones de módulos con enrichments de SI e Inflamación


pdf(file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Heatmap_modulos_todas_las_redes.pdf",height =8 ,width =22  )
heatmap.2( x = matriz_todas ,
           main = "Comparison All Enriched Modules In All Networks",   # heat map title
           density.info = "none",  # turns off density plot inside color legend
           trace = "none",         # turns off trace lines inside the heat map
           margins = c(15,2),     # widens margins around plot
           col = my_palette,       # use on color palette defined earlier
           labRow = "",            # change row and column labels
           #labCol = "",
           ylab = "GO Terms",
           xlab = "Module",
           #Rowv = as.dendrogram(row.cluster), # apply selected clustering method
           #Colv = as.dendrogram(col.cluster), # apply selected clustering method
           keysize = 0.8,
           key.title = "Enrichment",
           key.xlab = " -log10 adj p-val",
           colCol = mod_annot_bar,
           ColSideColors = col_annot_bar
           #,        # size of color key
           #ColSideColors = colsidecols ,dendrogram = "column",#Additional Options ## Color labeled columns
           #scale = "column" # cambiar por "column" o "row" segun se quiera la escala de colores normalizada en una direccion particular
           
)

par(lend = 1)           # square line ends for the color legend
legend("left",      # location of the legend on the heatmap plot
       legend = c("Basal", "Her2", "LumA","LumB","Normal"), # category labels
       col = c("VIOLET", "GREEN", "RED","ORANGE","BLUE"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)




dev.off()


###################################################################
# Heatmap version with colored dendrogram
###################################################################
# Create dendrogram with colored clusters:

library(gplots)
library(dendextend)
library(colorspace)

# distance & hierarchical clustering
class(matriz_todas)

dend1 <- as.dendrogram(hclust(dist(t(matriz_todas),method = "euclidean"))) #Matrix transposed to obtain columns clustering
c_group <- 30 # number of clusters
dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl) # add color to the lines
dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl)   # add color to the labels

# reorder the dendrogram, must incl. `agglo.FUN = mean`
cMeans <- colMeans(matriz_todas, na.rm = T)

# may skip reordering of the dendrogram.
dend1 <- reorder(dend1, colMeans(matriz_todas, na.rm = T), agglo.FUN = mean)

# get the color of the leaves (labels) for `heatmap.2`
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# if plot the dendrogram alone:
# the size of the labels:
dend1 <- set(dend1, "labels_cex", 0.3)

# Horizontal version
par(mar = c(4,1,1,5)) # c("down","left","up","right")
plot_horiz.dendrogram(dend1, side = F) # use side = T to horiz mirror if needed

# Vertical version
par(mar = c(5,5,1,1)) # c("down","left","up","right")
plot(dend1)

#Save just the dendrogram
pdf(file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Dendrograma_modulos_todas_las_redes.pdf",height =10 ,width =22  )
plot(dend1)
dev.off()


#######################################
# Add dendrogram to heatmap
# Heatmap de los resultados
my_palette <- colorRampPalette(c("white","blue","orange","red"))(n = 40)

#Basal   Her2   LumA   LumB Normal 
#324    248    273    273    345 

# Vector de colores para la barra de anotación que aparece arriba del heatmap
col_annot_bar = c( rep("VIOLET", 44), 
                   rep("GREEN", 49),
                   rep("RED", 48),
                   rep("ORANGE", 39),
                   rep("BLUE", 48)
)
# Vector de colores para marcar los nombres de las columnas que aparecen abajo del heatmap

mod_annot_bar = rep("BLACK",228) #rellenamos de "negro" el vector

mod_annot_bar[c(34,38,44,2,4,7,9,11,17,19,25,26,31,87,89,91,48,49,50,57,58,62,64,70,74,82,129,134,135,137,138,140,96,97,98,103,104,109,116,120,124,175,176,177,179,146,148,151,152,159,160,168,169,181,186,199,202,217,218,220,224,225,228,205,213)] <- "RED" # asignamos "rojo" a las posiciones de módulos con enrichments de SI e Inflamación


pdf(file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Heatmap_modulos_todas_las_redes_color_dendro.pdf",height =10 ,width =22  )
heatmap.2( x = matriz_todas ,
           main = "Comparison All Enriched Modules In All Networks",   # heat map title
           density.info = "none",  # turns off density plot inside color legend
           trace = "none",         # turns off trace lines inside the heat map
           margins = c(11,2),     # widens margins around plot
           lhei = c(1,2),
           col = my_palette,       # use on color palette defined earlier
           labRow = "",            # change row and column labels
           #labCol = "",
           ylab = "GO Terms",
           xlab = "Module",
           dendrogram = "col",
           #Rowv = as.dendrogram(row.cluster), # apply selected clustering method
           Colv = dend1, # apply selected clustering method
           keysize = 0.8,
           key.title = "Enrichment",
           key.xlab = " -log10 adj p-val",
           colCol = mod_annot_bar,
           ColSideColors = col_annot_bar,
           sepcolor="grey90",
           colsep = c(0,1,5,8,13,14,17,21,26,30,35,36,37,42,43,48,49,54,55,58,59,71,74,78,81,82,100,105,219,224,228) # position of separation bars between dendrogram clusters
           #,        # size of color key
           #ColSideColors = colsidecols ,dendrogram = "column",#Additional Options ## Color labeled columns
           #scale = "column" # cambiar por "column" o "row" segun se quiera la escala de colores normalizada en una direccion particular
           
)

par(lend = 1)           # square line ends for the color legend
legend("left",      # location of the legend on the heatmap plot
       legend = c("Basal", "Her2", "LumA","LumB","Normal"), # category labels
       col = c("VIOLET", "GREEN", "RED","ORANGE","BLUE"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

dev.off()








