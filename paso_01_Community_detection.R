###########################
### Community Detection ###
###########################

# Create an object that contains module structure for each network.
# Save R object to be used in steps to follow.

library(igraph) 

# Ubicacion de las redes
network_files <- c(Basal="Data/cuts_top100k_Basal.txt",Her2="Data/cuts_top100k_Her2.txt",LumA="Data/cuts_top100k_LumA.txt",LumB="Data/cuts_top100k_LumB.txt",Normal="Data/cuts_top100k_Sanos.txt")

# Cargamos todas las redes como data.frames en una sola lista
network_tables <- lapply(network_files, read.table,sep="\t", header=TRUE)
#length(network_tables)
#names(network_tables)

# Obtenemos las redes como objetos de igraph
networks <- lapply(network_tables,graph_from_data_frame,directed = FALSE)
#length(networks)
#names(networks)

# Eliminamos de la red las interacciones duplicadas
s_networks <- lapply(networks, simplify)

# Inferimos la estructura de comunidades utilizando el algoritmo INFOMAP
clusters_top100k <- lapply(s_networks, cluster_infomap,nb.trials=1000)
#length(clusters_top100k)
#names(clusters_top100k)

save(clusters_top100k,file = "Results/paso_01_Community_detection/clusters_top100k")

# Comparación de la composición de genes que integran los modulos de 10 o mas genes

#Obtener las listas de genes de los modulos para cada red

com_numbers <- list()
genelist_coms <- list()

com_numbers$Normal <- as.integer(names(sizes(clusters_top100k$Normal)[sizes(clusters_top100k$Normal)>9]))
genelist_coms$Normal <- names(membership(clusters_top100k$Normal)[membership(clusters_top100k$Normal)%in%com_numbers$Normal])

com_numbers$LumA <- as.integer(names(sizes(clusters_top100k$LumA)[sizes(clusters_top100k$LumA)>9]))
genelist_coms$LumA <- names(membership(clusters_top100k$LumA)[membership(clusters_top100k$LumA)%in%com_numbers$LumA])

com_numbers$LumB <- as.integer(names(sizes(clusters_top100k$LumB)[sizes(clusters_top100k$LumB)>9]))
genelist_coms$LumB <- names(membership(clusters_top100k$LumB)[membership(clusters_top100k$LumB)%in%com_numbers$LumB])

com_numbers$Her2 <- as.integer(names(sizes(clusters_top100k$Her2)[sizes(clusters_top100k$Her2)>9]))
genelist_coms$Her2 <- names(membership(clusters_top100k$Her2)[membership(clusters_top100k$Her2)%in%com_numbers$Her2])

com_numbers$Basal <- as.integer(names(sizes(clusters_top100k$Basal)[sizes(clusters_top100k$Basal)>9]))
genelist_coms$Basal <- names(membership(clusters_top100k$Basal)[membership(clusters_top100k$Basal)%in%com_numbers$Basal])


str(genelist_coms)

#Comparacion de las listas de genes
sum(genelist_coms$Normal %in% genelist_coms$LumA)
#[1] 6643
sum(genelist_coms$Normal %in% genelist_coms$LumB)
#[1] 6653
sum(genelist_coms$Normal %in% genelist_coms$Her2)
#[1] 6874
sum(genelist_coms$Normal %in% genelist_coms$Basal)
#[1] 6579
sum(genelist_coms$LumA %in% genelist_coms$LumB)
#[1] 8517
sum(genelist_coms$LumA %in% genelist_coms$Her2)
#[1] 8562
sum(genelist_coms$LumA %in% genelist_coms$Basal)
#[1] 8447
sum(genelist_coms$LumB %in% genelist_coms$Her2)
#[1] 8647
sum(genelist_coms$LumB %in% genelist_coms$Basal)
#[1] 8520
sum(genelist_coms$Her2 %in% genelist_coms$Basal)
#[1] 8538


pdf(file = "Community_sizes.pdf")
plot(table(sizes(clusters_top100k$Normal)),main= "Normal",xlab = "Module size",ylab = "Frequency")
plot(table(sizes(clusters_top100k$LumA)),main= "LumA",xlab = "Module size",ylab = "Frequency")
plot(table(sizes(clusters_top100k$LumB)),main= "LumB",xlab = "Module size",ylab = "Frequency")
plot(table(sizes(clusters_top100k$Her2)),main= "Her2",xlab = "Module size",ylab = "Frequency")
plot(table(sizes(clusters_top100k$Basal)),main= "Basal",xlab = "Module size",ylab = "Frequency")
dev.off()


