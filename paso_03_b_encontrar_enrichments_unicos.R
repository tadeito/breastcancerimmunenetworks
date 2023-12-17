#############################################################
## Encontrar procesos enriquecidos exclusivos por subtipos ##
#############################################################

library(gplots)
library(RColorBrewer)

# Los resultados en conjunto de los enriquecimientos estan en forma de una matriz de enriquecimientos con el GO term como rownames y el nombre-del-módulo+nombre-de-la-red  como colnames. la matriz se obtiene con el script "analisis_enriquecimientos_subtipos_mama.R" y tiene una version en texto en:  "results/enrichment_matrix_all_subtypes.txt"

# binarizar la matriz para mostrar presencia/ausencia de enriquecimeinto significativo para un proceso

enrichment_matrix <- matriz_todas_inflamacion

# convertir la matriz en ausencia/presencia (T/F)
binary_matrix <- !enrichment_matrix==0

rownames(binary_matrix)
colnames(binary_matrix)

# identificar cuantas veces se repite un proceso a lo largo de todos los modulos
process_reps <- apply(binary_matrix,1,sum)

# obtener los nombres de los procesos que se encuentran solo una vez
unique_enrichments <- names(process_reps[process_reps==1])

# filtrar la matriz por procesos unicos
unique_process_matrix <- binary_matrix[rownames(binary_matrix)%in%unique_enrichments,]

# identificar los modulos enriquecidos con procesos unicos
relevant_modules <- apply(unique_process_matrix,2,sum)
relevant_modules <- names(relevant_modules[!relevant_modules==0])

# filtrar la matriz por modulos relevantes
unique_process_matrix <- enrichment_matrix[rownames(enrichment_matrix)%in%unique_enrichments,colnames(enrichment_matrix)%in%relevant_modules]

# ordenar la matriz por numero de enriquecimientos por modulo
ordered_pm <- unique_process_matrix[,order(apply(!unique_process_matrix==0,2,sum),decreasing = FALSE)]

# Heatmap de los resultados

my_palette <- colorRampPalette(c("white","blue","orange","red"))(n = 40)


pdf(file = "/Users/CSBIG/Tadeo/Doctorado_Inflamacion_BC_subtipos/Results/paso_03_Analysis_of_enrichments/Heatmap_unicos_inmune_todas_las_redes.pdf",height = 25,width = 17)

heatmap.2( x = ordered_pm ,
           main = "Unique processes by module",   # heat map title
           density.info = "none",  # turns off density plot inside color legend
           trace = "none", # turns off trace lines inside the heat map
           margins = c(15,25),     # widens margins around plot
           col = my_palette,       # use on color palette defined earlier
           #labRow = "",            # change row and column labels
           #labCol = "",
           ylab = "GO Terms",
           xlab = "Module",
           Rowv = FALSE, # apply selected clustering method
           Colv = FALSE, # apply selected clustering method
           keysize = 0.8,
           key.title = "Enrichment",
           key.xlab = " -log10 adj p-val",
           #,        # size of color key
           #ColSideColors = colsidecols ,
           dendrogram = "none",#Additional Options ## Color labeled columns
           #scale = "column" # cambiar por "column" o "row" segun se quiera la escala de colores normalizada en una direccion particular
           sepwidth=c(0.01,0.01),
           sepcolor="black",
           colsep=0:ncol(ordered_pm), # where column sepparators must be placed
           #rowsep=1:nrow(ordered_pm)
           rowsep=c(0,2,3,8,9,10,12,16,26,27,28,29,30,34,37,57,59,62,71,72,73,74,86,87,93,94,97,98,103,104,105,125,132,178,182,183,184,186) # where row sepparators must be placed
           
)

dev.off()


  ################################################################################################################

# Comparacion por modulos

