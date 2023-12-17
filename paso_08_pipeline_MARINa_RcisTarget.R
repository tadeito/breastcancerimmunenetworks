######################
## RcisTarget-VIPER ##
######################

# This pipeline is intended to join the network deconvolution of ARACNe, community detection of INFOMAP, TF prediction of RcisTarget and TF priorization of VIPER.
## This script was put together by: Tadeo ##

# To start, you need an MI network where a cutoff point is already decided (by MI or p-value cutoff point).


###############
## Libraries ##
###############

library(igraph) # To handle network objects
library(RcisTarget) # Find transcriptional regulators through motif annotations  
library(viper) # Find TMRs 

##########################
## Networks preparation ###
##########################

# List network file locations, load them and convert them into igraph objects.

# Top 100k networks
#netfiles <- c(Normal="/labs/csbig/e/tadeo/corte_de_redes/data/mamaTCGA/completas/top100k_Sanos.txt",
#                   LumA="/labs/csbig/e/tadeo/corte_de_redes/data/mamaTCGA/completas/top100k_LumA.txt",
#                   LumB="/labs/csbig/e/tadeo/corte_de_redes/data/mamaTCGA/completas/top100k_LumB.txt",
#                   Her2="/labs/csbig/e/tadeo/corte_de_redes/data/mamaTCGA/completas/top100k_Her2.txt",
#                   Basal="/labs/csbig/e/tadeo/corte_de_redes/data/mamaTCGA/completas/top100k_Basal.txt")

# Networks cut by p-value
netfiles <- c(Normal="/labs/csbig/e/tadeo/corte_de_redes/sorted_Sanos_1e10-7.sif",
              LumA="/labs/csbig/e/tadeo/corte_de_redes/sorted_LumA_1e10-7.sif",
               LumB="/labs/csbig/e/tadeo/corte_de_redes/sorted_LumB_1e10-7.sif",
               Her2="/labs/csbig/e/tadeo/corte_de_redes/sorted_Her2_1e10-7.sif",
              Basal="/labs/csbig/e/tadeo/corte_de_redes/sorted_Basal_1e10-7.sif")



# Obtain a list with networks  as DF objects inside.
rawnets <- lapply(netfiles, read.table , stringsAsFactors = FALSE)

# To get a smaller cut:
#rawnets <- lapply(rawnets,head,20000) # remember interactions are duplicated, so 10k needs 20k as input before simplification

# Define a function to convert raw nets to undirected, simplified igraph network objects
net_from_DF <- function(x){ simplify(graph_from_data_frame(x[,c(1,3)],directed = FALSE))
  
}

# Obtain list with igraph networks
networks <- lapply(rawnets,net_from_DF)



############################################
## Infer community structure in networks ##
############################################


# Take networks list and obtain a list of clusters objects.

# Clusters are inferred with INFOMAP
communities <-lapply(networks, cluster_infomap)

# Save results to file
save(communities,file = "results/communities_1e-7.RData")

# En el trabajo de Genes asociados a Inflamación por subtipos usamos las comunidades de la red de 100,000 interacciones.
names(clusters_top100k)


communities <- clusters_top100k


# before starting the next step we need to generate gene lists from each cluster and for each network.

# Define a function to extract gene names from communities and make a list of genelists
comm_genelists <- function(x) {
  valid_comms <-as.integer(
  names(
    table(
      membership(x))[table(membership(x))>9] # Find all communities with 10 or more genes
  )
)
  genelists <- list()
  
for(i in 1:length(valid_comms)){
  genelists[[i]] <- names(membership(x)[membership(x)==valid_comms[i]])  
}
  names(genelists) <- as.character(valid_comms)
print(genelists)
  }

genelists_all_networks <- lapply(communities,comm_genelists)
save(genelists_all_networks,file = "results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/genelists_all_networks_top100k_10_o_mas_genes")

########################################################
## Infer regulators for each cluster for each network ##
########################################################

# This section is based on the workflor for RcisTarget that can be found at: "https://bioconductor.org/packages/devel/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html"

# Take clusters list, filter communities of desired size and infer regulators for each one
# Obtain a list of lists that contains the regulators by module by network

## Preparation
# To run the motif enrichment, we have to provide the lists of genes, the motiff ranking and motiff annotations.

# Motiff rankings
# These are contained in a .soft file and can be obtained from: "https://resources.aertslab.org/cistarget/"
# It is important that the rankings are generated from the same reference as the data, in our case it its GRCh38

#motifRankings <- importRankings("data/RcisTarget/db/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
motifRankings <- importRankings("Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")

# Motiff annotations
# Motiff annotations are included in RcisTarget and can be called with the following:
data(motifAnnotations_hgnc)

# Run the motif enrichment (complete version with details hidden)

# Run motif enrichments over all networks
motifEnrichments_all_networks <- lapply(genelists_all_networks,cisTarget ,motifRankings,motifAnnot=motifAnnotations_hgnc)
save(motifEnrichments_all_networks,file="results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/motifEnrichments_500bp_all_networks_top100k.RData")

# To free memory space, remove motifRankings database
rm(motifRankings)

#############################################################
## Convert regulators table into weighted directed network ##
#############################################################

# Take regulators table and extract a regulon.
# Regulon must have a weight.
# Final format might be taken by VIPER 
# Output is a list of regulons networks

GetSIF <- function(x){ #starts function
  
# Get motif enrichments table
tftable <- data.frame(x,stringsAsFactors = FALSE)

# Clean table from motifs that are not associated to a TF (annotated as "" or empty) 
tftable <-tftable[!(tftable[,5]==""),]
# For some reason, some motifs show no annotated targets and need to be cleaned up too
tftable <-tftable[!(tftable[,9]==""),]

# tftable[r,5] # TF names
# tftable[r,9] # Target names

# Empty DF to store all regulons
Regulon <- data.frame(TF=character(),target=character(),motif=character(),NES=numeric(),cluster=character(),stringsAsFactors = FALSE)

#r <- 1 # row counter
for(r in 1:length(tftable$motif)){ # open loop by rows

  TFs <- unlist(strsplit(gsub(" ","",gsub("\\s*\\([^\\)]+\\).","",tftable[r,5])),"[;]")) # Sepparate TF names 
  
#tf <- 1 # TF in row counter
for(tf in 1:length(TFs)){ # loop repeated by tfs

Targets <- unlist(strsplit(unlist(tftable[r,9]),"[;]"))

# Local DF to store regulon
miniRegulon <- data.frame(TF=character(length(Targets)),target=character(length(Targets)),motif=character(length(Targets)),NES=numeric(length(Targets)),cluster=character(length(Targets)),stringsAsFactors = FALSE)

miniRegulon[,] <- NA

#tg <-1 # target in row counter
for(tg in 1:length(Targets)){
  miniRegulon$TF[tg] <- TFs[tf]
  miniRegulon$target[tg] <- Targets[tg]
  miniRegulon$motif[tg] <- tftable[r,2]
  miniRegulon$NES[tg] <- tftable[r,3]
  miniRegulon$cluster[tg] <- tftable[r,1]
  
    
} # cierra loop por targets
Regulon <- rbind(Regulon,miniRegulon) # Add miniregulon to regulons network
} # cierra loop por tfs
} # close loop by rows

return(Regulon) # Print onject to be captured as output

} # ends function GetSIF()

# Obtain results and save them into target object 
SIF_todas <- lapply(motifEnrichments_all_networks,GetSIF)

#######################
## Collapse regulons ##
#######################

# Regulons obtained in the preceeding step potentially contain multiple interactions between a TF and a target gene through many motifs.
# we collapse multiple interactions combining each attribute in to one


simplify_regulon <- function(x){

    g <- graph_from_data_frame(x,directed = TRUE) # get directed igraph object
  g_simplified <- simplify( g, edge.attr.comb =list(motif=toString,NES="sum",cluster="random")) # collapse multiple edges
    # Multiple "motif" attributes are concatenated in a string to keep the reference.
    # The NES score was used as a measure of the likelihood of the interaction in the context of the regulon. We choose to sum NES values when a TF is linked to a gene by multiple motifs. 
  # the "cluster" attribute is chosen at random because a gene can belong only to one cluster.
  SIF_simplificado <- cbind(as.data.frame(get.edgelist(g_simplified)),as.data.frame(edge_attr(g_simplified))) # re-convert to data.frame
  return(SIF_simplificado) # This is the output of the function  
}

SIFs_colapsados <- lapply(SIF_todas,simplify_regulon)

############################################
## Filter out TFs not expressed in tissue ##
############################################

# Obtain the list of all genes in the expression matrix
#expressed_genes <- read.table("data/expressed_genes_brca_TCGA_symbol.txt", header = FALSE,stringsAsFactors=FALSE)

# Alternativamente, si tenemos la matriz de expresión cargada podemos obtener los genes expresados a partir de los rownames de la matriz:

expressed_genes <- data.frame(symbol=rownames(eset),stringsAsFactors = FALSE)
  
# Define function to filter the data 
filter_expressed <- function(x) {
  x[as.character(x[,1]) %in% expressed_genes[,1] & as.character(x[,2]) %in% expressed_genes[,1],]  
}

# Batch filter
SIFs_colapsados_clean <- lapply(SIFs_colapsados,filter_expressed)

# We also want to express weights inthe network with values between 0 and 1. for this we re-scale the NES column dividing by the largest (accumulated) NES value
rescale_NES <- function(x){
x[,4] <- x[,4]/max(x[,4])
return(x)
}

SIFs_colapsados_clean <- lapply(SIFs_colapsados,rescale_NES)

# Save object to file

save(SIFs_colapsados_clean,file = "results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/SIF_top100k_500bp_todas.RData")

# Write individual regulon table files (For use in CYTOSCAPE)
write.table(SIFs_colapsados_clean$Normal,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Normal.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(SIFs_colapsados_clean$LumA,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_LumA.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(SIFs_colapsados_clean$LumB,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_LumB.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(SIFs_colapsados_clean$Her2,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Her2.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(SIFs_colapsados_clean$Basal,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Basal.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)



# Regulons doesn't seem to be working in the .3col format. I'm going to try and convert them to .adj format
adj_3col <- list()

for(i in 1:4){
adj <- SIFs_colapsados_clean[[i]][,c(1,2,4)]
adj <-adj[order(adj[,1]),]
adj <- adj[adj[,1]%in%expressed_genes[,1],]
adj_3col[[i]] <- adj
}

names(adj_3col)<- c("Basal","Her2","LumA","LumB")
# Write individual regulons in 3col format to use as MARINa input

write.table(adj_3col[[3]],file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_LumA.3col",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(adj_3col[[4]],file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_LumB.3col",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(adj_3col[[2]],file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Her2.3col",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(adj_3col[[1]],file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Basal.3col",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)


####################
## VIPER (MARINa) ##  
####################
 
# Rank TFs according to their relevance in the differential expression patterns

# Rank TFs according to their relevance in the differential expression patterns

# VIPER uses four inputs: a regulons network, an expression dataset, a molecular signature and a null model.
# Regulons network was obtained in the preceding stages
# Expression sets need to be loaded
# Molecular signatures and null models are generated with the expression sets

# Recover regulons from RcisTarget:
load("Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/SIF_top100k_500bp_todas.RData")
# Define expression set locations

Eset_locations <- c(Normal="Data/Normal_PBCMC_ARSYM.txt",
                    LumA="Data/LumA_PBCMC_ARSYM.txt",
                    LumB="Data/LumB_PBCMC_ARSYM.txt",
                    Her2="Data/Her2_PBCMC_ARSYM.txt",
                    Basal="Data/Basal_PBCMC_ARSYM.txt")

# Load esets
Esets_all <- lapply(Eset_locations,read.table,header=TRUE,row.names=1, stringsAsFactors=FALSE)
# Convert Esets to "matrix" to make them readable in VIPER

Esets_all <- lapply(Esets_all,as.matrix)

# Verify that all gene names (rownames) are identical and in the same order in all datasets:

order_rows_to_match <- function(x){
  x <- x[order(rownames(x)),]
  return(x)
}

Esets_all <- lapply(Esets_all,order_rows_to_match)

# > all.equal(row.names(Esets_all$Normal),row.names(Esets_all$LumA))
# [1] TRUE
# > all.equal(row.names(Esets_all$Normal),row.names(Esets_all$LumB))
# [1] TRUE
# > all.equal(row.names(Esets_all$Normal),row.names(Esets_all$Her2))
# [1] TRUE
# > all.equal(row.names(Esets_all$Normal),row.names(Esets_all$Basal))
# [1] TRUE



# Generate expression signatures
signatures <- lapply(Esets_all[-1], function(x) rowTtest(x, y = Esets_all$Normal))

# z-score values for the GES to maintain consistency with the null model

signature_znormalize <- function(x) {
  signature <- (qnorm(x$p.value/2, lower.tail=F) * sign(x$statistic))[, 1]
  return(signature)
}

zsignatures <- lapply(signatures,signature_znormalize)


# NULL model by sample permutations
#los modelos nulos nos sirven para comparar si lo que encontramos en nuestros casos es azaroso 
# o hay un fenomeno subyacente en la siguiente funcion de viper se hacen 2000 permutaciones per=1000 a partir de nuestros datos repos=T que se ropngan los datos T

nmods <- function(x){    
  nullmodel <- ttestNull(x,Esets_all$Normal , per=2000, repos=T)
  return(nullmodel)
}

nullmodels <- lapply(Esets_all[-1],nmods)

# Generate the regulons

files_3col <- c("foo",
                "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_LumA.3col",
                "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_LumB.3col",
                "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Her2.3col",
                "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/regulon_top100k_500bp_Basal.3col"
)


# This version of the function joins Normal+Tumor esets to generate a dset
# This is important because in the inference of the regulons, a mixture distribution inference is performed and requires both data sources ## !!! Preguntar a Enrique !!! ##
# Este paso resuelve el problema de que no haya convergencia.
# Usamos los regulones a partir de los SIF pesados por que esta funcion ordena las interacciones basado en los valores de interacción TF - regulador

#****Continuar apartir de aqui ***
# Resolver problema para cargar los regulones:
#Loading the dataset...
#Generating the regulon objects...
#Error in postprobs[, w] : subscript out of bounds
# Este problema se resualve filrando adecuadamente los genes del regulón, eliminando aquellos que no se encuentran en la matriz de expresión.


Regulons <- list()
for(n in 2:5){
  eset <- cbind(Esets_all[[n]],Esets_all$Normal) 
  regulon <- aracne2regulon(afile = files_3col[n],eset = eset ,format="3col")
  Regulons[[n-1]] <- regulon
}

names(Regulons) <- c("LumA","LumB","Her2","Basal")

########  M A R I N a ######## 
#  MARINa uses the molecular signature, regulon and null model to rank TFs by their efect on gene expression changes between two conditions.


TMRs <- list() # To store results

for(n in 1:4){  
  marina <- msviper(zsignatures[[n]], Regulons[[n]] , nullmodels[[n]])
  TMRs[[n]] <- marina
}

names(TMRs) <- c("LumA","LumB","Her2","Basal")

summaries_TMRs <- lapply(TMRs,summary,mrs=30)

names(summaries_TMRs) <- c("LumA","LumB","Her2","Basal")


# Ploting masters regulators
for(n in 1:4){
pdf(paste0("Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/MARINaTop10_",names(TMRs)[n],"_top100k_RcisTarget500pb.pdf"),width=6,height=7)
plot(TMRs[[n]],mrs=10, cex=.7)
dev.off()
}

save(TMRs,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/MARINa_500bp_top100k.RData")


##########################################################
## Preparar los regulones para visualizar en Cytoscape ##
##########################################################


library(igraph)

class(clusters)
names(clusters) # Aqui se encuentran los modulos a que pertence cada gen en las redes
#[1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"


names(SIFs_colapsados_clean) # Este es el regulon que vamos a filtrar
#[1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"


Selecccionar_modulos <- list(
  Basal=c("49","78","96","102","104", "110", "112", "115" ,"134" ,"137" ,"156"  ,"157" ,"184" ),
  Her2=c( "77" ,"88" ,"93" ,"104" ,"107" ,"111" ,"126" ,"127","133","136","143","166","228" ),
  LumA=c("40","75","80","83","84","98","111","112","114","122","125","131","143","160","178" ),
  LumB=c("69","70","85","93","111","116","119","120","130","134","161","175" ),
  Normal=c("10","13","21","26","56","58","61","86","88","97","283","392")
)

n=1
g=1


Regulon_inmunes <- list()
for(n in 1:5){ # Abre loop por redes
Modulo <- character() # Identificamos el modulo de la red al que pertenece cada uno de los targets del regulon
for(g in 1:length(SIFs_colapsados_clean[[n]]$V2)){ # abre loop por interacciones
Modulo[g] <-  membership(clusters[[n]])[which(names(membership(clusters[[n]])) %in% SIFs_colapsados_clean[[n]]$V2[g])]

} # Cierra loop interacciones
Regulon_inmunes[[n]] <- SIFs_colapsados_clean[[n]][Modulo %in% Selecccionar_modulos[[n]],]
print(n)
} # Cierra loop por redes

names(Regulon_inmunes) <- names(SIFs_colapsados_clean)

# Escribir los regulones en disco
getwd()

for(n in 1:5){
write.table(Regulon_inmunes[[n]],file = paste("Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/Regulon_Clusters_Inflamacion",names(Regulon_inmunes)[n],".tsv"),sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
}


###########################################################
# Identificar los TFs compartidos por los subtipos #
########################################################

names(SIFs_colapsados_clean)
# [1] "Basal"  "Her2"   "LumA"   "LumB"   "Normal"

class(SIFs_colapsados_clean[[1]]$V1)
# [1] "factor" 
str(SIFs_colapsados_clean[[1]])



length(unique(c(
  as.character(SIFs_colapsados_clean[[1]]$V1),
  as.character(SIFs_colapsados_clean[[2]]$V1),
  as.character(SIFs_colapsados_clean[[3]]$V1),
  as.character(SIFs_colapsados_clean[[4]]$V1),
  as.character(SIFs_colapsados_clean[[5]]$V1)
  )))
# [1] 1849 reguladores asociados.


combined_TFs <- sort(unique(c(
  as.character(SIFs_colapsados_clean[[1]]$V1),
  as.character(SIFs_colapsados_clean[[2]]$V1),
  as.character(SIFs_colapsados_clean[[3]]$V1),
  as.character(SIFs_colapsados_clean[[4]]$V1),
  as.character(SIFs_colapsados_clean[[5]]$V1)
)))
# Tabla para almacenar los resultados
# indico de antemano cuantos renglones tiene para irlos llenando sin que marque error de dimensiones.
TF_comparisons_and_pertenence <- data.frame(TF=character(1849), #1
                                            in_Basal=numeric(1849), #2
                                            in_Her2=numeric(1849), #3
                                            in_LumA=numeric(1849), #4
                                            in_LumB=numeric(1849), #5
                                            in_Normal=numeric(1849), #6
                                            combination=numeric(1849), #7
                                            target_count_Basal=numeric(1849), #8
                                            target_count_Her2=numeric(1849), #9
                                            target_count_LumA=numeric(1849), #10
                                            target_count_LumB=numeric(1849), #11
                                            target_count_Normal=numeric(1849), #12
                                            reached_modules_Basal=numeric(1849), #13
                                            reached_modules_Her2=numeric(1849), #14
                                            reached_modules_LumA=numeric(1849), #15
                                            reached_modules_LumB=numeric(1849), #16
                                            reached_modules_Normal=numeric(1849), #17
                                            stringsAsFactors = FALSE
                                            )


r=1
for(r in 1:1849){
TF_comparisons_and_pertenence$TF[r] <- combined_TFs[r]
TF_comparisons_and_pertenence$in_Basal[r] <- ifelse(combined_TFs[r] %in% as.character(SIFs_colapsados_clean[[1]]$V1),2,1)
TF_comparisons_and_pertenence$in_Her2[r]<- ifelse(combined_TFs[r] %in% as.character(SIFs_colapsados_clean[[2]]$V1),3,1)
TF_comparisons_and_pertenence$in_LumA[r]<- ifelse(combined_TFs[r] %in% as.character(SIFs_colapsados_clean[[3]]$V1),5,1)
TF_comparisons_and_pertenence$in_LumB[r]<- ifelse(combined_TFs[r] %in% as.character(SIFs_colapsados_clean[[4]]$V1),7,1)
TF_comparisons_and_pertenence$in_Normal[r]<- ifelse(combined_TFs[r] %in% as.character(SIFs_colapsados_clean[[5]]$V1),11,1)
TF_comparisons_and_pertenence$combination[r] <- prod(unlist(TF_comparisons_and_pertenence[r,c(2,3,4,5,6)]))
TF_comparisons_and_pertenence$target_count_Basal[r] <- sum(as.character(SIFs_colapsados_clean[[1]]$V1) %in% combined_TFs[r])
TF_comparisons_and_pertenence$target_count_Her2[r]<- sum(as.character(SIFs_colapsados_clean[[2]]$V1) %in% combined_TFs[r])
TF_comparisons_and_pertenence$target_count_LumA[r]<- sum(as.character(SIFs_colapsados_clean[[3]]$V1) %in% combined_TFs[r])
TF_comparisons_and_pertenence$target_count_LumB[r]<- sum(as.character(SIFs_colapsados_clean[[4]]$V1) %in% combined_TFs[r])
TF_comparisons_and_pertenence$target_count_Normal[r]<- sum(as.character(SIFs_colapsados_clean[[5]]$V1) %in% combined_TFs[r])
TF_comparisons_and_pertenence$reached_modules_Basal[r] <- length(unique(as.character(SIFs_colapsados_clean[[1]]$cluster)[as.character(SIFs_colapsados_clean[[1]]$V1) %in% combined_TFs[r]])) 
TF_comparisons_and_pertenence$reached_modules_Her2[r]<- length(unique(as.character(SIFs_colapsados_clean[[2]]$cluster)[as.character(SIFs_colapsados_clean[[2]]$V1) %in% combined_TFs[r]]))
TF_comparisons_and_pertenence$reached_modules_LumA[r]<- length(unique(as.character(SIFs_colapsados_clean[[3]]$cluster)[as.character(SIFs_colapsados_clean[[3]]$V1) %in% combined_TFs[r]]))
TF_comparisons_and_pertenence$reached_modules_LumB[r]<- length(unique(as.character(SIFs_colapsados_clean[[4]]$cluster)[as.character(SIFs_colapsados_clean[[4]]$V1) %in% combined_TFs[r]]))
TF_comparisons_and_pertenence$reached_modules_Normal[r]<- length(unique(as.character(SIFs_colapsados_clean[[5]]$cluster)[as.character(SIFs_colapsados_clean[[5]]$V1) %in% combined_TFs[r]]))
}

table(TF_comparisons_and_pertenence$combination)
#2    3    5    6   11   14   15   21   22   30   33   42   55   66   70   77  105  110  154  165  210  231  330  385 
#1    3    3    1    3    1    5    1   10    4   18    1   10   31    1    9    3   10    4   25    7   21   63    9 
#462  770 1155 2310 
#41   20   71 1473 

pdf("Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/TF_target_count_vs_modules_reached.pdf")
plot(x=TF_comparisons_and_pertenence$target_count_Basal,y=TF_comparisons_and_pertenence$reached_modules_Basal)
plot(x=TF_comparisons_and_pertenence$target_count_Her2,y=TF_comparisons_and_pertenence$reached_modules_Her2)
plot(x=TF_comparisons_and_pertenence$target_count_LumA,y=TF_comparisons_and_pertenence$reached_modules_LumA)
plot(x=TF_comparisons_and_pertenence$target_count_LumB,y=TF_comparisons_and_pertenence$reached_modules_LumB)
plot(x=TF_comparisons_and_pertenence$target_count_Normal,y=TF_comparisons_and_pertenence$reached_modules_Normal)
dev.off()

write.table(TF_comparisons_and_pertenence,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/TF_comparisons_and_pertenence.txt",sep="\t",quote = FALSE, row.names = FALSE, col.names = TRUE)
