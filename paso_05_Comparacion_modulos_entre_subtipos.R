######################################################
## Similitud de modulos entre las diferentes redes ###
######################################################

# Comparar los modulos entre redes para ver si las composiciones de genes son similares, y determinar si exiten modulos exclusivos de alguno de los subtipos
# Medida de similtud usando el indice de Jaccard

# los modulos a comparar se encuentran en los resultados del paso 01

library(igraph)


load("Results/paso_01_Community_detection/clusters_top100k")

# Definimos la función para comparar los modulos
jaccard <- function(a,b){
  sum(a %in% b)/length(unique(c(a,b)))
}



# Obtenemos la tabla de correspondencias de genes a modulos

s=1 # subtipo
m=1 # modulo

big_modulos <- data.frame(Symbol=character(),module=character(),subtype=character(),stringsAsFactors = FALSE)

for(s in 1:length(clusters_top100k)){

mods <- data.frame(Symbol=names(membership(clusters_top100k[[s]])),
                          module= paste(membership(clusters_top100k[[s]]),
                                        names(clusters_top100k)[s]
                                        
                                        ),
                   subtype = names(clusters_top100k)[s],
                   stringsAsFactors = FALSE
)
big_modulos <- rbind(big_modulos,mods)

}

# Filtramos la tabla para quedarnos solo con aquellos modulos que tengan diez genes o mas

# Encontramos los modulos con 10 genes o mas
mayores_a_10 <- names(table(big_modulos$module)[table(big_modulos$module)>9])

# Filtramos la tabla con los modulos identificados
modulos_filtrados <- big_modulos[big_modulos$module %in% mayores_a_10,]

#Generamos la lista de modulos
nombres_modulos <- unique(modulos_filtrados$module)

#Generamos una matriz vacía para almacenar los resultados:
matriz_similitud_modulos <- matrix(nrow = length(nombres_modulos),ncol = length(nombres_modulos))
rownames(matriz_similitud_modulos) <- nombres_modulos
colnames(matriz_similitud_modulos) <- nombres_modulos

# Hacemos las comparaciones de la composicion de genes entre modulos y las almacenamos como entradas en la matriz

i=1 # renglones
j=1 # columnas

# llenamos las entradas de la matriz
for(i in 1:length(nombres_modulos)){ #loop por renglones
for(j in 1:length(nombres_modulos)){ #loop por columnas
a<-modulos_filtrados$Symbol[modulos_filtrados$module == nombres_modulos[i]]
b<-modulos_filtrados$Symbol[modulos_filtrados$module == nombres_modulos[j]]
matriz_similitud_modulos[i,j] <- jaccard(a,b)
}
}

# Guardamos la matriz como objeto de R en los resultados
save(matriz_similitud_modulos,file = "Results/paso_05_similitud_de_modulos/matriz_similitud_modulos.Rdata")

colnames(matriz_similitud_modulos)

##########################################################################################
### Exportar la similitud como una lista de interacciones para visualizar en cytoscape ###
##########################################################################################

library(igraph)

# Revisamos que la matriz se haya completado correctamente
isSymmetric.matrix(matriz_similitud_modulos)
# [1] TRUE

sum(!diag(matriz_similitud_modulos) == 1)
# [1] 0

# Para convertir a edgelist la matriz de similitud utilizamos igraph, que nos permite omitir los edges con valor de cero y las autointeracciones, así como eliminar las que se encuentras duplicadas:

# Convertimos la matriz a objeto de igraph
g_modulos <- graph_from_adjacency_matrix(adjmatrix = matriz_similitud_modulos,mode = "undirected",diag = FALSE,weighted = TRUE)

# Eliminamos las interacciones duplicadas
simplify(g_modulos)

# Convertimos a lista de interacciones con formato (from, to, weight)

edgelist_similitud_modulos <- as_data_frame(g_modulos,what = "edges")

# necesitamos un atributo para los edges para poder mapear facilmente el color al nivel de similitud 
edgelist_similitud_modulos$similarity_cut <- ifelse(edgelist_similitud_modulos$weight >0.9,1,
       ifelse(edgelist_similitud_modulos$weight > 0.8, 2,
              ifelse(edgelist_similitud_modulos$weight > 0.7, 3,
                     ifelse(edgelist_similitud_modulos$weight > 0.6, 4,
                            ifelse(edgelist_similitud_modulos$weight > 0.5, 5,
                                   ifelse(edgelist_similitud_modulos$weight > 0.4, 6,
                                          ifelse(edgelist_similitud_modulos$weight > 0.3, 7,
                                                 ifelse(edgelist_similitud_modulos$weight > 0.2, 8,
                                                        ifelse(edgelist_similitud_modulos$weight > 0.1, 9,
                     10)))))))))

# Verificamos que la columna de corte se haya añadido de forma correcta
edgelist_similitud_modulos


# Verificamos que se tenga un data.frame. Esto para poder escribir el archivo en disco en formato de texto
class(edgelist_similitud_modulos)

# Escribimos el edgelist en la carpeta de resultados
write.table(edgelist_similitud_modulos,file = "Results/paso_05_similitud_de_modulos/edgelist_similitud_modulos_subtipos.txt",sep = "\t", quote = FALSE,row.names = FALSE)


#############################################################
### Generar la tabla de atributos visuales para cytoscape ###
#############################################################

big_modulos # Aqui se encuentra la informacion de pertenencia de los genes a los modulos

# Calculamos el tamaño de cada modulo y generamos un data.frame con los resultados
size_modulos <- as.data.frame(table(big_modulos$module),stringsAsFactors=FALSE)

class(size_modulos)
# [1] "data.frame"
class(size_modulos$Var1)
# [1] "character"

dim(size_modulos)
# [1] 1076    2

colnames(size_modulos)
# [1] "Var1" "Freq"

# Filtramos aquellos modulos con diez genes o más
size_modulos <- size_modulos[size_modulos$Freq > 9,]
# Cambiamos el nombre de la columna para usarlo en la tabla de cytoscape
colnames(size_modulos)[2] <- "module size" 

# Generamos una vector con el subtipo de origen de cada modulo
subtype_origin <- character()
for(m in 1:length(size_modulos$Var1)){
subtype_origin[m]<- big_modulos$subtype[which(big_modulos$module %in% size_modulos$Var1[m])[1]]
}

# añadimos la columna con el subtipo de origen de cada modulo
size_modulos$subtype_origin <- subtype_origin

# Escribimos la tabla de atributos que usaremso en cytoscape.
write.table(size_modulos,file = "Results/paso_05_similitud_de_modulos/Cytoscape_size_modulos_subtipos.txt",sep = "\t", quote = FALSE,row.names = FALSE)

#######################################################################################
## Analisis de agrupamiento para determinar conjuntos de módulos similares entre subtipos ###
#######################################################################################










