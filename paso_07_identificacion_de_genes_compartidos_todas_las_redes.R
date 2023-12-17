############################################################################
## Paso 07 identificación de nucleo de genes compartidos entre las redes ##
############################################################################

# Utilizamos la tabla con los nombres de los genes y el modulo de procedencia obtenida en el paso 05
# Esta tabla contiene los genes que se encuentran en las comunidades de 10 o mas genes

head(modulos_filtrados)
# Symbol    module subtype
# 7  LRRC37A2  41 Basal   Basal
# 8   LRRC37A  41 Basal   Basal
# 9   POLR2J3 158 Basal   Basal
# 10  POLR2J2 158 Basal   Basal
# 11   UPK3BL 158 Basal   Basal
# 12   CTAGE9 148 Basal   Basal

# Queremos los genes que se encuentren en todas las redes, basados en el resultado de los modulos con alta similitud 

length(unique(modulos_filtrados$Symbol))

table(modulos_filtrados$Symbol) # seleccionamos los genes que se repiten 5 veces

comunes_todos <- names(table(modulos_filtrados$Symbol)[table(modulos_filtrados$Symbol)==5])

# Comparación entre los genes considerados en todas las redes y los compartidos por todas 
length(unique(modulos_filtrados$Symbol))
# [1] 14326
length(comunes_todos)
# [1] 4963

#filtrar por genes presentes en todas las redes
genes_compartidos <- modulos_filtrados[modulos_filtrados$Symbol%in%comunes_todos,]
# ordenar para facilitar la visualizacion
genes_compartidos <- genes_compartidos[order(genes_compartidos$Symbol),]

write.table(comunes_todos,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/shared_genes_all_networks.txt", row.names = FALSE, quote = FALSE,sep = "\t")

write.table(genes_compartidos,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/shared_genes_by_community.txt", row.names = FALSE, quote = FALSE,sep = "\t")

# Ver que ocurre a nivel de los módulos. ¿Cuantos de los genes que los componen son compartidos

compartidos_por_modulo <- as.data.frame(table(modulos_filtrados$module),stringsAsFactors = FALSE)
dim(compartidos_por_modulo)
# [1] 1076    2
compart <- as.data.frame(table(genes_compartidos$module),stringsAsFactors = FALSE)
dim(compart)
# [1] 989   2

# el numero de comunidades es distinto si tomamos a los genes compartidos solamente.

i=5

restantes <- integer()
for(i in 1:1076){
restantes[i] <- ifelse(compartidos_por_modulo$Var1[i]%in%compart$Var1,
compart$Freq[which(compart$Var1==compartidos_por_modulo$Var1[i])],
0)
}

compartidos_por_modulo$compartidos <- restantes

shared_percent <- numeric()
percent <- function(x,y){x/y*100}
for(i in 1:1076){
shared_percent[i] <-  percent(compartidos_por_modulo$compartidos[i],compartidos_por_modulo$Freq[i]) 
  }
compartidos_por_modulo$shared_percent <- shared_percent

# cambio de nombre de las columnas
names(compartidos_por_modulo) <- c("module","size","shared_genes","shared_percentage")

#añadir una columna donde se indique el subtipo para hacer filtros en la tabla

st <- strsplit(compartidos_por_modulo$module,split = " ")
i=1

st_or <- character()
for(i in 1:length(st)){
st_or[i] <-  st[[i]][[2]]
  }

compartidos_por_modulo$subtipo <- st_or

compartidos_por_modulo <- compartidos_por_modulo[order(compartidos_por_modulo$shared_percentage,decreasing = TRUE),]

write.table(compartidos_por_modulo,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/percentages_of_sharing.txt", row.names = FALSE, quote = FALSE,sep = "\t")

pdf(file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/plot_size_vs_percentage.pdf")
plot(x=compartidos_por_modulo$size,y=compartidos_por_modulo$shared_percentage)
dev.off()


#####
# Algunas preguntas sobre las comunidades
####



# ¿En cuantas comunidades se reparten los genes compartidos por todas las redes?

length(unique(genes_compartidos$module[genes_compartidos$subtype=="Normal"]))
#  [1] 318

length(unique(genes_compartidos$module[genes_compartidos$subtype=="LumA"]))
# [1] 167

length(unique(genes_compartidos$module[genes_compartidos$subtype=="LumB"]))
# [1] 157

length(unique(genes_compartidos$module[genes_compartidos$subtype=="Her2"]))
# [1] 186

length(unique(genes_compartidos$module[genes_compartidos$subtype=="Basal"]))
# [1] 161


#########################################################
### Identificar la estructura común a todas las redes ###
#########################################################
library(igraph)

# Ubicacion de las redes
network_files <- c(Basal="Data/cuts_top100k_Basal.txt",Her2="Data/cuts_top100k_Her2.txt",LumA="Data/cuts_top100k_LumA.txt",LumB="Data/cuts_top100k_LumB.txt",Normal="Data/cuts_top100k_Sanos.txt")
# Cargamos todas las redes como data.frames en una sola lista
network_tables <- lapply(network_files, read.table,sep="\t", header=TRUE,stringsAsFactors =FALSE)
# Obtenemos las redes como objetos de igraph
networks <- lapply(network_tables,graph_from_data_frame,directed = FALSE)
# Eliminamos de la red las interacciones duplicadas
s_networks <- lapply(networks, simplify)


# Obtener los nodos y aristas que se encuentran en todas las redes y generar una red con atributos para visualizar de gorma combinada.


# Limpiar las redes, dejando solo los genes en comunidades de 10 o mas miembros

clusters <- clusters_top100k
# Verificar el orden de las redes en las listas
names(clusters)
# [1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"
names(s_networks)
#[1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"

x=1
extraer_genes <- function(x){
modulos_ <- as.numeric(names(table(membership(clusters[[x]]))[table(membership(clusters[[x]]))>=10]))
genes_ <- names(membership(clusters[[x]])[membership(clusters[[x]])%in%modulos_])
return(names(V(s_networks[[x]])[names(V(s_networks[[x]]))%in%genes_]))
}

n <- c(1,2,3,4,5)

genes_coms_ <- lapply(n,extraer_genes)
str(genes_coms_)

# Combinar los nombres de todos los genes de todas las redes 

str(genes_coms_)
#List of 5
#$ : chr [1:9981] "LRRC37A2" "LRRC37A" "POLR2J3" "POLR2J2" ...
#$ : chr [1:10658] "INS-IGF2" "IGF2" "LRRC37A2" "LRRC37A" ...
#$ : chr [1:10229] "NOMO3" "NOMO2" "XIST" "TSIX" ...
#$ : chr [1:10177] "CKMT1B" "CKMT1A" "LRRC37A2" "LRRC37A" ...
#$ : chr [1:9253] "PALM2-AKAP2" "AKAP2" "HSPA1B" "HSPA1A" ...

genes_combinados <- unique( c(genes_coms_[[1]],
                              genes_coms_[[2]],
                              genes_coms_[[3]],
                              genes_coms_[[4]],
                              genes_coms_[[5]]
  ))

length(genes_combinados)
#[1] 14326

# tabla de pertencia combinada de los genes
library(matrixStats)

Pertenencia_genes <- data.frame(symbol=genes_combinados,stringsAsFactors = FALSE)

Pertenencia_genes$Basal <- ifelse(Pertenencia_genes$symbol%in%genes_coms_[[1]],2,1
       )
Pertenencia_genes$Her2 <- ifelse(Pertenencia_genes$symbol%in%genes_coms_[[2]],3,1
)
Pertenencia_genes$LumA <- ifelse(Pertenencia_genes$symbol%in%genes_coms_[[3]],5,1
)
Pertenencia_genes$LumB <- ifelse(Pertenencia_genes$symbol%in%genes_coms_[[4]],7,1
)
Pertenencia_genes$Normal <- ifelse(Pertenencia_genes$symbol%in%genes_coms_[[5]],11,1
)

Pertenencia_genes$combinacion <- rowProds(as.matrix(Pertenencia_genes[,-1]))  

write.table(Pertenencia_genes,"Results/paso_07_comparacion_genes_compartidos_todas_las_redes/Pertenecia_combinada.txt",sep = "\t",quote = FALSE, row.names = FALSE,col.names = TRUE)


table(Pertenencia_genes$combinacion)
#2    3    5    6    7   10   11   14   15   21   22   30   33   35   42   55   66   70   77  105  110  154  165  210 
#211  412  321  144  246  115  808  109  186  182  223  160  437  121  183  319  202  166  249  203  141  155  205 2314 
#231  330  385  462  770 1155 2310 
#229  265  161  307  323  266 4963 


######################################################################
## Analisis de las subredes de las comunidades con 10 o mas genes ###
######################################################################

# hacer la lista de interacciones posibles
# La estrategia a seguir es:
# 1.- obtener las interacciones unicas de cada red
x=1
extraer_subgrafos <- function(x){
  modulos_ <- as.numeric(names(table(membership(clusters[[x]]))[table(membership(clusters[[x]]))>=10]))
  genes_ <- names(membership(clusters[[x]])[membership(clusters[[x]])%in%modulos_])
  nodes_ <- V(s_networks[[x]])[names(V(s_networks[[x]]))%in%genes_]
    return(induced_subgraph(s_networks[[x]],nodes_))
}

class(nodes_)
#[1] "igraph.vs"

genes_

n <- c(1,2,3,4,5)

c_networks <- lapply(n,extraer_subgrafos)
names(c_networks) <- names(s_networks)


# 2.- combinar todas las redes en una sola red
interacciones_combinadas <- data.frame(get.edgelist(c_networks[[1]]),stringsAsFactors = FALSE)
interacciones_combinadas <- rbind(interacciones_combinadas,data.frame(get.edgelist(c_networks[[2]]),stringsAsFactors = FALSE))
interacciones_combinadas <- rbind(interacciones_combinadas,data.frame(get.edgelist(c_networks[[3]]),stringsAsFactors = FALSE))
interacciones_combinadas <- rbind(interacciones_combinadas,data.frame(get.edgelist(c_networks[[4]]),stringsAsFactors = FALSE))
interacciones_combinadas <- rbind(interacciones_combinadas,data.frame(get.edgelist(c_networks[[5]]),stringsAsFactors = FALSE))
dim(interacciones_combinadas)

# 3.- eliminar las interacciones multiples de la red combinada
interacciones_combinadas_g <- graph_from_edgelist(as.matrix(interacciones_combinadas),directed = FALSE)
interacciones_combinadas_g <- simplify(interacciones_combinadas_g)
interacciones_combinadas <- data.frame(get.edgelist(interacciones_combinadas_g),stringsAsFactors = FALSE)


# 4.- identificar cada interaccion como presente/ausente en cada una de las redes.


# tabla de presencia de interacciones para las cinco redes y su combinación
i=1
n=1

dim(interacciones_combinadas)

valores_pertenencia <- c(Basal=2,Her2=3,LumA=5,LumB=7,Normal=11)
class(pertenencia_interacciones)

for(n in 1:5){
  pertenencia_interacciones <- rep(0, length(interacciones_combinadas[,1])) #vector para almacenar los resultados

    for(i in 1:length(interacciones_combinadas[,1])){
      pertenencia_interacciones[i]<- ifelse(interacciones_combinadas[i,1]%in%genes_coms_[[n]]&interacciones_combinadas[i,2]%in%genes_coms_[[n]],
                          ifelse(
  are.connected(s_networks[[n]],interacciones_combinadas[i,1],interacciones_combinadas[i,2]), valores_pertenencia[n],1
),

1)

      }#asignar valores segun si la interacción existe o no en la red
interacciones_combinadas[,n+2]<-pertenencia_interacciones
}

names(interacciones_combinadas) <- c("from","to","Basal","Her2","LumA","LumB","Normal")
interacciones_combinadas$combinacion <- rowProds(as.matrix(interacciones_combinadas[,c(-1,-2)]))
#

#quality control: verificar que todas las interacciones esten al menos dentro de una de las redes (tengan vlor combinado =1)
sum(interacciones_combinadas$combinacion==1) 
# [1] 0

# Distinguir las interacciones inter-modulo de las interacciones intra-modulo

#Donde estan las pertenencias de los genes de las comunidades mayores a 10 
library(igraph)
clusters # contiene la información de las comunidades en cada red

x=1 # Contador por redes
i=1 # Contador por interaccion
for(x in 1:5){
modulos_ <- as.numeric(names(table(membership(clusters[[x]]))[table(membership(clusters[[x]]))>=10]))
genes_ <- names(membership(clusters[[x]])[membership(clusters[[x]])%in%modulos_])
pertenencia_ <- membership(clusters[[x]])
status_ <- character() # Vector para almacenar los resultados
for(i in 1:length(interacciones_combinadas$from)){
status_[i]<- ifelse(
!interacciones_combinadas[i,x+2]==1, # identificar si la interacción existe en la red
  ifelse(pertenencia_[which(names(pertenencia_)==interacciones_combinadas$from[i])] ==
pertenencia_[which(names(pertenencia_)==interacciones_combinadas$to[i])], "intra", "inter") # identificar si ambos genes estan en la misma comunidad
, "absent")
}
interacciones_combinadas[,x+8]<- status_ # anexar el status de la red a la tabla de interacciones
}

names(interacciones_combinadas)[9] <- "status_Basal"
names(interacciones_combinadas)[10] <- "status_Her2"
names(interacciones_combinadas)[11] <- "status_LumA"
names(interacciones_combinadas)[12] <- "status_LumB"
names(interacciones_combinadas)[13] <- "status_Normal"

head(interacciones_combinadas)

write.table(interacciones_combinadas,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/interacciones_combinadas.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)

##################################################################
## Preguntas a la tabla de interacciones #########################
##################################################################
# Estas son las subredes de las comunidades con 10 o mas genes

# La codificación de la pertenencia de interacciones es la siguiente:
# Basal = 2
# Her2 = 3
# LumA = 5
# LumB = 7
# Normal = 11
# Si la interaccion pertenece a mas de una red, su codificación es la multiplicación de sus pertenencias combinadas, por ejemplo, si una interacción se encuentre en Basal y en Her2 pero no en las otras, su codificacion es 2*3 = 6

# ¿Cuantas interacciones son exclusivas de la red de Normales?
sum(interacciones_combinadas$combinacion==11)
# [1] 80828

# ¿Cuantas interacciones son exclusivas de la red de Basal?
sum(interacciones_combinadas$combinacion==2)
#[1] 25234

# ¿Cuantas interacciones son exclusivas de la red de Her2?
sum(interacciones_combinadas$combinacion==3)
#[1] 28548

# ¿Cuantas interacciones son exclusivas de la red de LumA?
sum(interacciones_combinadas$combinacion==5)
#[1] 26322

# ¿Cuantas interacciones son exclusivas de la red de LumB?
sum(interacciones_combinadas$combinacion==7)
#[1] 26041

# ¿Cuantas interacciones son compartidas por todas las redes?
sum(interacciones_combinadas$combinacion==2310)
#[1] 1141

# ¿Cuantas interacciones son compartidas por todos los subtipos tumorales pero no en normales?
sum(interacciones_combinadas$combinacion==210)
#[1] 28140

# Obtener las redes para visualizacion
red_pertenencia_exclusivos_Basal <- interacciones_combinadas[interacciones_combinadas$combinacion==2,]
red_pertenencia_exclusivos_Her2 <- interacciones_combinadas[interacciones_combinadas$combinacion==3,]
red_pertenencia_exclusivos_LumA <- interacciones_combinadas[interacciones_combinadas$combinacion==5,]
red_pertenencia_exclusivos_LumB <- interacciones_combinadas[interacciones_combinadas$combinacion==7,]
red_pertenencia_exclusivos_Normal <- interacciones_combinadas[interacciones_combinadas$combinacion==11,]
red_pertenencia_compartidos_todos <- interacciones_combinadas[interacciones_combinadas$combinacion==2310,]
red_pertenencia_compartidos_tumores<- interacciones_combinadas[interacciones_combinadas$combinacion==210,]

# Escribir los archivos de redes en el disco duro
write.table(red_pertenencia_exclusivos_Basal,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_exclusivos_Basal.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)
write.table(red_pertenencia_exclusivos_Her2,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_exclusivos_Her2.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)
write.table(red_pertenencia_exclusivos_LumA,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_exclusivos_LumA.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)
write.table(red_pertenencia_exclusivos_LumB,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_exclusivos_LumB.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)
write.table(red_pertenencia_exclusivos_Normal,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_exclusivos_Normal.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)
write.table(red_pertenencia_compartidos_todos,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_compartidos_todos.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)
write.table(red_pertenencia_compartidos_tumores,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/red_pertenencia_compartidos_tumores.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)

data.frame(t(table(interacciones_combinadas$combinacion)))

class(interacciones_combinadas$combinacion)

hist(table(interacciones_combinadas$combinacion))
hist(interacciones_combinadas$Basal)


# Preguntas a las comunidades de la red


#¿Cuantas interacciones son intra-comunidad y cuantas inter-comunidad en cada una de las redes?

table(interacciones_combinadas$status_Basal)
# absent  inter  intra 
# 186231   7130  89527 

table(interacciones_combinadas$status_Her2)
# absent  inter  intra 
# 189335  10614  82939 

table(interacciones_combinadas$status_LumA)
# absent  inter  intra 
# 187375   8120  87393

table(interacciones_combinadas$status_LumB)
# absent  inter  intra 
# 186748   7574  88566 

table(interacciones_combinadas$status_Normal)
# absent  inter  intra 
# 198604  33688  50596 


compartidos_todos <- interacciones_combinadas$combinacion==2310
compartidos_tumores <- interacciones_combinadas$combinacion == 210

# Interacciones compartidas por todas las redes en el contexto de la comunidades
table(interacciones_combinadas$status_Basal[compartidos_todos])
# inter intra 
# 22  1119 
table(interacciones_combinadas$status_Her2[compartidos_todos])
# inter intra 
# 25  1116 
table(interacciones_combinadas$status_LumA[compartidos_todos])
# inter intra 
# 29  1112
table(interacciones_combinadas$status_LumB[compartidos_todos])
# inter intra 
# 34  1107 
table(interacciones_combinadas$status_Normal[compartidos_todos])
# inter intra 
# 285   856



# Interacciones compartidas por los tumores pero no en los normales en el contexto de la comunidades
table(interacciones_combinadas$status_Basal[compartidos_tumores])
# inter intra 
# 304 27836 
table(interacciones_combinadas$status_Her2[compartidos_tumores])
# inter intra 
# 325 27815 
table(interacciones_combinadas$status_LumA[compartidos_tumores])
# inter intra 
# 341 27799
table(interacciones_combinadas$status_LumB[compartidos_tumores])
# inter intra 
# 435 27705 
table(interacciones_combinadas$status_Normal[compartidos_tumores])
# absent 
# 28140 


# Interacciones exclusivas de cada red en el contexto de las comunidades
table(interacciones_combinadas$status_Basal[interacciones_combinadas$combinacion == 2])
# inter intra 
# 4888 20346 
table(interacciones_combinadas$status_Her2[interacciones_combinadas$combinacion == 3])
# inter intra 
# 8583 19965 
table(interacciones_combinadas$status_LumA[interacciones_combinadas$combinacion == 5])
# inter intra 
# 5878 20444
table(interacciones_combinadas$status_LumB[interacciones_combinadas$combinacion == 7])
# inter intra 
# 5158 20883 
table(interacciones_combinadas$status_Normal[interacciones_combinadas$combinacion == 11])
# inter intra 
# 32744 48084 

######################################################################################################
## Comparar los genes compartidos/exclusivos de cada red con los genes anotados de respuesta inmune ##
######################################################################################################

genelist_immune_system_process <- read.table("Results/paso_07_comparacion_genes_compartidos_todas_las_redes/genelist_GO_Immune_System_Process.txt", stringsAsFactors = FALSE, col.names = FALSE)
dim(genelist_immune_system_process)
class(genelist_immune_system_process[,1])



x=2

intersectar_genes <- function(x){
genes_i <-  Pertenencia_genes$symbol[Pertenencia_genes$combinacion==x]
  
return(genes_i[genes_i%in%genelist_immune_system_process[,1]])
}

intersecciones <- c(Basal=2,Her2=3,LumA=5,LumB=7,todos_tumores=210,todas_redes=2310)


interseccion_combinacion_immune_system_process <- lapply(intersecciones,intersectar_genes)

dim(genelist_immune_system_process)

str(interseccion_combinacion_immune_system_process)


#####################################################################################
## Identificar la estructura de los módulos en comparación con las otras redes ######
#####################################################################################

library(igraph)

# Identificar a que modulo pertenecen las interacciones intra-comunidad en cada red

# identificar cuales módulos son distintos y cuales son parecidos

head(interacciones_combinadas)
names(interacciones_combinadas)
#  [1] "from"          "to"            "Basal"         "Her2"          "LumA"          "LumB"          "Normal"       
# [8] "combinacion"   "status_Basal"  "status_Her2"   "status_LumA"   "status_LumB"   "status_Normal"

dim(interacciones_combinadas)
# [1] 282888     13

i=1
r=1
status_red <- c(9,10,11,12,13)
comm_interaccion <- c(14,15,16,17,18)

for(r in 1:5){
comunidad_intra_interacc <- integer()
for(i in 1:282888){
comunidad_intra_interacc[i] <- ifelse(
interacciones_combinadas[i,status_red[r]]=="intra",
membership(clusters_top100k[[r]])[which(names(membership(clusters_top100k[[r]]))%in%interacciones_combinadas[i,1])],
0)
}
interacciones_combinadas[,comm_interaccion[r]] <- comunidad_intra_interacc
}

names(interacciones_combinadas)[c(14,15,16,17,18)] <- c("comm_Basal","comm_Her2","comm_LumA","comm_LumB","comm_Normal")
# [1] "from"          "to"            "Basal"         "Her2"          "LumA"          "LumB"          "Normal"       
# [8] "combinacion"   "status_Basal"  "status_Her2"   "status_LumA"   "status_LumB"   "status_Normal" "comm_Basal"   
# [15] "comm_Her2"     "comm_LumA"     "comm_LumB"     "comm_Normal"  

# Actualizar la tabla de interacciones combinadas con la información de a que módulo pertenecen en cada red
write.table(interacciones_combinadas,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/interacciones_combinadas.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)

# Obtener la proporcion de interacciones exclusivas de cada modulo
r=1
phenotype <- c("Basal","Her2","LumA","LumB","Normal") 
comm_red <- c(14,15,16,17,18)
status_red <- c(9,10,11,12,13)
phenocode <- c(2,3,5,7,11)

m=1

mod_interactions_summary <- data.frame(phenotype=character(),module=integer(),interaction_count=integer(),exclusive=integer(),all_networks=integer(),all_tumors=integer(),exclusive_precentage=numeric(),stringsAsFactors = FALSE)


for(r in 1:5){
  
  mod_list <- sort(unique(interacciones_combinadas[,comm_red[r]]))[-1]
  
      for(m in 1:length(mod_list)){
mod_sum <- data.frame(phenotype=character(1),module=integer(1),interaction_count=integer(1),exclusive=integer(1),all_networks=integer(1),all_tumors=integer(1),exclusive_precentage=numeric(1),stringsAsFactors = FALSE)

  
  mod_intlist <- interacciones_combinadas$combinacion[interacciones_combinadas[,comm_red[r]]==mod_list[m]]
mod_summary <- table(mod_intlist)

mod_sum$phenotype[1]<-phenotype[r]
mod_sum$module <- mod_list[m] 
mod_sum$interaction_count[1] <- length(mod_intlist)
mod_sum$exclusive[1]<- ifelse(phenocode[r]%in%as.integer(names(mod_summary)) ,mod_summary[as.integer(names(mod_summary))==phenocode[r]] ,0)
mod_sum$all_networks[1]<- ifelse(2310%in%as.integer(names(mod_summary)) ,mod_summary[as.integer(names(mod_summary))==2310],0)
mod_sum$all_tumors[1]<- ifelse(210%in%as.integer(names(mod_summary)) ,mod_summary[as.integer(names(mod_summary))==210],0)
mod_sum$exclusive_precentage <- (mod_sum$exclusive/mod_sum$interaction_count)*100

mod_interactions_summary <- rbind(mod_interactions_summary,mod_sum)
}
}

head(mod_interactions_summary)

# añadir los porcentajes de interacciones compartidas por todas las redes y por todos los tumores
mod_interactions_summary$all_network_percentage <- (mod_interactions_summary$all_networks/mod_interactions_summary$interaction_count)*100
mod_interactions_summary$all_tumors_percentage <- (mod_interactions_summary$all_tumors/mod_interactions_summary$interaction_count)*100

write.table(mod_interactions_summary,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/comparacion_interacciones_exclusivas_intracomunidad.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

# ¿Qué módulos tienen al menos un 50% de composición de interacciones exclusivas?
mod_interactions_summary[mod_interactions_summary$exclusive_precentage>=50,]

# Rastrear que pasa con los modulos enriquecidos en inflamación:
filtro_inflamacion_fenotipo <- c(rep("Basal",13),rep("Her2",13),rep("LumA",15),rep("LumB",12),rep("Normal",12))
filtro_inflamacion_modulo <- c(49,78,96,102,104,110,112,115,134,137,156,157,184,77,88,93,104,107,111,126,127,133,136,143,166,228,40,75,80,83,84,98,111,112,114,122,125,131,143,160,178,69,70,85,93,111,116,119,120,130,134,161,175,10,13,21,26,56,58,61,86,88,97,283,392)
length(filtro_inflamacion_fenotipo)
length(filtro_inflamacion_modulo)

n=1
SI_modules <- data.frame(stringsAsFactors = FALSE)
for(n in 1:65){
SI_modules_ <-mod_interactions_summary[which(mod_interactions_summary$phenotype==filtro_inflamacion_fenotipo[n]&mod_interactions_summary$module==filtro_inflamacion_modulo[n]),]
SI_modules <-rbind(SI_modules,SI_modules_)
}

SI_modules_relevant <- SI_modules[SI_modules$exclusive_precentage>=50,]
# phenotype module interaction_count exclusive all_networks all_tumors exclusive_precentage
#  Basal    134                24        13            0          0             54.16667
#  Basal    157                12         9            0          0             75.00000
#  Basal    184                 9         9            0          0            100.00000
#  Her2    126                24        17            0          1             70.83333
#  Her2    133                20        16            0          0             80.00000
#  Her2    136                17        12            0          0             70.58824
#  Her2    166                10         8            0          0             80.00000
#  LumA    131                20        14            0          0             70.00000
#  LumA    178                11        10            0          0             90.90909
#  LumB     85               103        67            0          4             65.04854
#  LumB    175                11         6            0          0             54.54545
#  Normal     10               996       890           30          0             89.35743
#  Normal     13               871       752           11          0             86.33754
#  Normal     21               277       214           11          0             77.25632
#  Normal     56                51        40            4          0             78.43137
#  Normal     88                23        23            0          0            100.00000
#  Normal    283                10        10            0          0            100.00000
#  Normal    392                 9         8            0          0             88.88889

write.table(SI_modules_relevant,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/modulos_relevantes.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

Si_modules_conserved_tumors <- SI_modules[SI_modules$all_tumors_percentage>=50,]
# [1] phenotype              module                 interaction_count      exclusive              all_networks          
# [6] all_tumors             exclusive_precentage   all_network_percentage all_tumors_percentage 
# <0 rows> (or 0-length row.names)

Si_modules_conserved_all <- SI_modules[SI_modules$all_network_percentage>=50,]

# phenotype module interaction_count exclusive all_networks all_tumors   exclusive_precentage   all_network_percentage
# 705      LumB    120                59         2           32          0             3.389831               54.23729
# 841    Normal     86                36         2           25          0             5.555556               69.44444
# 852    Normal     97                32        11           17          0            34.375000               53.12500
# all_tumors_percentage
# 705                     0
# 841                     0
# 852                     0

############################################################################################################################
### Comparación de los módulos a nivel de interacciones
############################################################################################################################

# Utilizando la información de a que módulo pertenecen las interacciones, determinar la similitud por indice de jaccard de los módulos entre las distintas redes


head(interacciones_combinadas)

# from        to Basal Her2 LumA LumB Normal combinacion status_Basal status_Her2 status_LumA status_LumB status_Normal
# 1 LRRC37A2   LRRC37A     2    3    5    7      1         210        intra       intra       intra       intra        absent
# 2 LRRC37A2    ARL17A     2    3    5    7      1         210        intra       intra       intra       intra        absent
# 3 LRRC37A2   PLEKHM1     2    1    5    1      1          10        intra      absent       intra      absent        absent
# 4 LRRC37A2    KANSL1     2    3    5    1      1          30        intra       inter       intra      absent        absent
# 5 LRRC37A2 CRHR1-IT1     2    1    5    7      1          70        intra      absent       intra       intra        absent
# 6 LRRC37A2      DND1     2    1    1    7      1          14        inter      absent      absent       intra        absent
# comm_Basal comm_Her2 comm_LumA comm_LumB comm_Normal
# 1         41       195        66       163           0
# 2         41       195        66       163           0
# 3         41         0        66         0           0
# 4         41         0        66         0           0
# 5         41         0        66       163           0
# 6          0         0         0       163           0

# En la tabla, las comunidades estan identificadas con cero si los genes unidos son de distintos modulos y con el numero de comunidad si estan unidos en el mismo modulo. 

coms_basal <- sort(unique(interacciones_combinadas$comm_Basal))[-1]
coms_her2 <- sort(unique(interacciones_combinadas$comm_Her2))[-1]
coms_luma <- sort(unique(interacciones_combinadas$comm_LumA))[-1]
coms_lumb <- sort(unique(interacciones_combinadas$comm_LumB))[-1]
coms_normal <- sort(unique(interacciones_combinadas$comm_Normal))[-1]

r=1
s=3

##### Forma sucia y facil de pensar #########

###################
## Basal vs Her2 ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_basal)){
for(s in 1:length(coms_her2)){
  intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[s],]))
total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[s],]))))
sim_index <- intersection/total

if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones

Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
Jaccard_row[1,1] <- paste0("Basal_",coms_basal[r])
  Jaccard_row[1,2] <- paste0("Her2_",coms_her2[s])
  Jaccard_row[1,3] <- sim_index

Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
}}

Jaccard_BasalvsHer2 <- Jaccard_estructura


###################
## Basal vs LumA ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_basal)){
  for(s in 1:length(coms_luma)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[s],]))))
    sim_index <- intersection/total

    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones    
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("Basal_",coms_basal[r])
    Jaccard_row[1,2] <- paste0("LumA_",coms_luma[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_BasalvsLumA <- Jaccard_estructura

###################
## Basal vs LumB ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_basal)){
  for(s in 1:length(coms_lumb)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[s],]))))
    sim_index <- intersection/total

    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
        
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("Basal_",coms_basal[r])
    Jaccard_row[1,2] <- paste0("LumB_",coms_lumb[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_BasalvsLumB <- Jaccard_estructura

###################
## Basal vs Normal ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_basal)){
  for(s in 1:length(coms_normal)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Basal==coms_basal[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("Basal_",coms_basal[r])
    Jaccard_row[1,2] <- paste0("Normal_",coms_normal[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_BasalvsNormal <- Jaccard_estructura

###################
## Her2 vs LumA ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_her2)){
  for(s in 1:length(coms_luma)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("Her2_",coms_her2[r])
    Jaccard_row[1,2] <- paste0("LumA_",coms_luma[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_Her2vsLumA <- Jaccard_estructura

###################
## Her2 vs LumB ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_her2)){
  for(s in 1:length(coms_lumb)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("Her2_",coms_her2[r])
    Jaccard_row[1,2] <- paste0("LumB_",coms_lumb[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_Her2vsLumB <- Jaccard_estructura

###################
## Her2 vs Normal ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_her2)){
  for(s in 1:length(coms_normal)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_Her2==coms_her2[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("Her2_",coms_her2[r])
    Jaccard_row[1,2] <- paste0("Normal_",coms_normal[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_Her2vsNormal <- Jaccard_estructura


###################
## LumA vs LumB ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_luma)){
  for(s in 1:length(coms_lumb)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("LumA_",coms_luma[r])
    Jaccard_row[1,2] <- paste0("LumB_",coms_lumb[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_LumAvsLumB <- Jaccard_estructura

###################
## LumA vs Normal ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_luma)){
  for(s in 1:length(coms_normal)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_LumA==coms_luma[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("LumA_",coms_luma[r])
    Jaccard_row[1,2] <- paste0("Normal_",coms_normal[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_LumAvsNormal <- Jaccard_estructura

###################
## LumB vs Normal ##
###################

# Iniciamos un data.frame vacío para almacenar los datos
Jaccard_estructura <- data.frame(from=character(),to=character(),jaccard_similitude=numeric(),stringsAsFactors = FALSE)

for(r in 1:length(coms_lumb)){
  for(s in 1:length(coms_normal)){
    intersection <- sum(rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[r],]) %in% rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))
    total <- length(unique(c(rownames(interacciones_combinadas[interacciones_combinadas$comm_LumB==coms_lumb[r],]), rownames(interacciones_combinadas[interacciones_combinadas$comm_Normal==coms_normal[s],]))))
    sim_index <- intersection/total
    
    if(sim_index==0)  next # Descartamos los pares de modulos que no comparten interacciones
    
    Jaccard_row <- data.frame(from=character(1),to=character(1),jaccard_similitude=numeric(1),stringsAsFactors = FALSE)
    Jaccard_row[1,1] <- paste0("LumB_",coms_lumb[r])
    Jaccard_row[1,2] <- paste0("Normal_",coms_normal[s])
    Jaccard_row[1,3] <- sim_index
    
    Jaccard_estructura <- rbind(Jaccard_estructura,Jaccard_row) # unimos el resultado a la tabla
    print(c(r,s)) # mostramos en pantalla los valores de r y s para visualizar el avance
  }}

Jaccard_LumBvsNormal <- Jaccard_estructura

# unir los resultados
Similitud_por_interacciones <- rbind(Jaccard_BasalvsHer2,Jaccard_BasalvsLumA,Jaccard_BasalvsLumB,Jaccard_BasalvsNormal,Jaccard_Her2vsLumA,Jaccard_Her2vsLumB,Jaccard_Her2vsNormal,Jaccard_LumAvsLumB,Jaccard_LumAvsNormal,Jaccard_LumBvsNormal)

write.table(Similitud_por_interacciones,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/Similitud_por_interacciones.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)


# Obtener la tabla de atributos para cytoscape

mods_2 <- sort(unique(c(Similitud_por_interacciones$from,Similitud_por_interacciones$to)))

n=1
fenos <- strsplit(mods_2,"_")

phenos <- character(length(fenos))
for(n in 1:length(fenos)){
phenos[n]<-fenos[[n]][1]
}

cytoscape_attr <- cbind(mods_2,phenos)

write.table(cytoscape_attr,file = "Results/paso_07_comparacion_genes_compartidos_todas_las_redes/Cytoscape_attr_Similitud_por_interacciones.txt",sep = "\t",quote = FALSE,col.names = TRUE, row.names = FALSE)


