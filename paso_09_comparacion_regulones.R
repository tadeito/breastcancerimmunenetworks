# Comparacion estructura de los regulones
# Identificar que interacciones se encuentran compartidas entre las redes. Determinar numericamente el parecido entre los regulones de los diferentes subtipos de cancer de mama


str(SIFs_colapsados_clean) # Contiene los regulones de todos los subtipos 

names(SIFs_colapsados_clean)

# Obtener el universo de interacciones de los regulones.

# Unimos todos los regulones en uno solo
int_universe <- rbind(SIFs_colapsados_clean[[1]][,c(1,2)],SIFs_colapsados_clean[[2]][,c(1,2)],SIFs_colapsados_clean[[3]][,c(1,2)],SIFs_colapsados_clean[[4]][,c(1,2)],SIFs_colapsados_clean[[5]][,c(1,2)])

# Simplificamos eliminando las interacciones duplicadas
int_universe <- unique(int_universe)

# Generamos una tabla para almacenar las comparaciones
Comparacion_regulones <- data.frame(Basal=numeric(719153),Her2=numeric(719153),LumA=numeric(719153),LumB=numeric(719153),Normal=numeric(719153))

# La estrategia para saber si una interacción se encuentra en uno de los regulones consiste en unir, para cada interacciòn, elnombre de el regulador y su blanco como una cadena de caracteres. luego se identifica si esa cadena se encuentra en el universo de interacciones. El procesamiento interno de R de cada función hace el cálculo mucho mas rapido que con for loops. 
univ_ref <- paste0(int_universe$V1,int_universe$V2) # Referencia de interacciones del universo

Basal_ref <- paste0(SIFs_colapsados_clean[[1]]$V1,SIFs_colapsados_clean[[1]]$V2) # referencia de interacciones en cada regulon
Her2_ref <- paste0(SIFs_colapsados_clean[[2]]$V1,SIFs_colapsados_clean[[2]]$V2)
LumA_ref <- paste0(SIFs_colapsados_clean[[3]]$V1,SIFs_colapsados_clean[[3]]$V2)
LumB_ref <- paste0(SIFs_colapsados_clean[[4]]$V1,SIFs_colapsados_clean[[4]]$V2)
Normal_ref <- paste0(SIFs_colapsados_clean[[5]]$V1,SIFs_colapsados_clean[[5]]$V2)



univ_ref %in% Basal_ref

# Hacemos las comparaciones y escribimos el resultado en la tabla.
Comparacion_regulones$Basal <- ifelse(univ_ref %in% Basal_ref ,2,1)
Comparacion_regulones$Her2 <- ifelse(univ_ref %in% Her2_ref ,3,1)
Comparacion_regulones$LumA <- ifelse(univ_ref %in% LumA_ref ,5,1)
Comparacion_regulones$LumB <- ifelse(univ_ref %in% LumB_ref ,7,1)
Comparacion_regulones$Normal <- ifelse(univ_ref %in% Normal_ref ,11,1)

Comparacion_regulones$Combinacion <- apply(Comparacion_regulones,MARGIN = 1,FUN = prod)

table(Comparacion_regulones$Combinacion)

#         2      3      5      6      7     10     11     14     15     21     22     30     33     35     42     55     66     70   77    105    110    154    165
#[1,]  66296  96283  74588  11997  67956  13087 169652  15217  11621  13124   6633  10873   7391  12056  11760   6365   2543  12550   6076  10627   2342   2832   2007

#        210    231    330    385    462    770   1155   2310
#[1,]  49278   2105   3105   2041   3045   3145   2329  20229

# completamos la tabla con las columnas de regulador - blanco

Comparacion_regulones <- cbind(int_universe,Comparacion_regulones)
names(Comparacion_regulones)[1] <- "Regulador"
names(Comparacion_regulones)[2] <- "Blanco"

# Escribir la tabla en disco


write.table(Comparacion_regulones,file = "Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/Comparacion_Regulones.txt",sep = "\t",quote = FALSE, row.names = FALSE)


##########################################################################################
# ¿Que reguladores se encuentran asociados a las interacciones compartidas por todas las redes?
# ¿Que reguladores definen las interacciones no compartidas por las redes?


Reguladores_no_compartidas <- table( c(
as.character(unique(Comparacion_regulones$Regulador[Comparacion_regulones$Combinacion==2])),
             as.character(unique(Comparacion_regulones$Regulador[Comparacion_regulones$Combinacion==3])),
                          as.character(unique(Comparacion_regulones$Regulador[Comparacion_regulones$Combinacion==5])),
                                       as.character(unique(Comparacion_regulones$Regulador[Comparacion_regulones$Combinacion==7])),
                                                    as.character(unique(Comparacion_regulones$Regulador[Comparacion_regulones$Combinacion==11]))
))

# No compartidos entre regulones
names(Reguladores_no_compartidas[Reguladores_no_compartidas==1])
#[1] "ANHX"   "AVEN"   "CANX"   "DMRTB1" "GPAM"   "IL24"   "NXPH3"  "PKM"    "SEBOX"  "TCEAL2" "ZNF763"

# Compartidos en todos los regulones, pero con blancos distintos
names(Reguladores_no_compartidas[Reguladores_no_compartidas==5])

# 1386 reguladores distintos


Comparacion_regulones[Comparacion_regulones$Regulador=="ANHX",]
#Regulador         Blanco Basal Her2 LumA LumB Normal Combinacion
# 1205102      ANHX TMEM256-PLSCR3     1    1    1    1     11          11
# 1205103      ANHX           MSH6     1    1    1    1     11          11
# 1205104      ANHX         SLC6A2     1    1    1    1     11          11
# 1205105      ANHX         ANGPT2     1    1    1    1     11          11
# 1205106      ANHX          CRHBP     1    1    1    1     11          11
# 1205107      ANHX          EGLN2     1    1    1    1     11          11
# 1205108      ANHX          FSCN1     1    1    1    1     11          11
# 1205109      ANHX        GALNT11     1    1    1    1     11          11
# 1205110      ANHX         SDHAF4     1    1    1    1     11          11
# 1205111      ANHX          TCEB3     1    1    1    1     11          11
# 1205112      ANHX           BCL9     1    1    1    1     11          11
# 1205113      ANHX        SLC30A1     1    1    1    1     11          11
Comparacion_regulones[Comparacion_regulones$Regulador=="AVEN",]
# Regulador Blanco Basal Her2 LumA LumB Normal Combinacion
# 729434      AVEN FBXO32     1    1    5    1      1           5
# 729435      AVEN   LAG3     1    1    5    1      1           5
# 729436      AVEN  H2AFJ     1    1    5    1      1           5
# 729437      AVEN   RERG     1    1    5    1      1           5
Comparacion_regulones[Comparacion_regulones$Regulador=="CANX",]
# Regulador  Blanco Basal Her2 LumA LumB Normal Combinacion
# 493149       CANX   GPSM1     1    3    1    1     11          33
# 1202569      CANX    MYRF     1    1    1    1     11          11
# 1202570      CANX    MSRA     1    1    1    1     11          11
# 1202571      CANX  FAM46A     1    1    1    1     11          11
# 1202573      CANX    PRR7     1    1    1    1     11          11
# 1202574      CANX   PAQR6     1    1    1    1     11          11
# 1202575      CANX  MINPP1     1    1    1    1     11          11
# 1202576      CANX PIKFYVE     1    1    1    1     11          11
# 1202577      CANX  PRSS23     1    1    1    1     11          11
# 1202578      CANX     ME3     1    1    1    1     11          11
# 1202579      CANX   APBB3     1    1    1    1     11          11
# 1202580      CANX SLC6A11     1    1    1    1     11          11
# 1202581      CANX FAM102A     1    1    1    1     11          11
# 1202582      CANX SLC15A1     1    1    1    1     11          11
Comparacion_regulones[Comparacion_regulones$Regulador=="DMRTB1",]
# Regulador Blanco Basal Her2 LumA LumB Normal Combinacion
# 1160745    DMRTB1  PDE7A     1    1    1    1     11          11
# 1160746    DMRTB1  SCN1B     1    1    1    1     11          11
# 1160747    DMRTB1   TPM3     1    1    1    1     11          11
# 1160748    DMRTB1  DOCK8     1    1    1    1     11          11
# 1160749    DMRTB1  CLCN3     1    1    1    1     11          11
# 1160750    DMRTB1 PFKFB2     1    1    1    1     11          11
# 1160751    DMRTB1   BCO2     1    1    1    1     11          11
# 1160752    DMRTB1   DEF8     1    1    1    1     11          11
# 1160753    DMRTB1   MPP3     1    1    1    1     11          11
# 1160754    DMRTB1 ZSWIM4     1    1    1    1     11          11
# 1160755    DMRTB1 CEP131     1    1    1    1     11          11
Comparacion_regulones[Comparacion_regulones$Regulador=="GPAM",]
# Regulador  Blanco Basal Her2 LumA LumB Normal Combinacion
# 232912      GPAM    ETV5     2    1    1    1      1           2
# 232913      GPAM    ACP5     2    1    1    1      1           2
# 232914      GPAM     FYB     2    1    1    1      1           2
# 232915      GPAM  LILRA2     2    1    1    1      1           2
# 232916      GPAM   LPAR5     2    1    1    1      1           2
# 232917      GPAM   NAPSB     2    1    1    1      1           2
# 232918      GPAM   DACT3     2    1    1    1      1           2
# 232919      GPAM    GAPT     2    1    1    1      1           2
# 232920      GPAM  P2RY13     2    1    1    1      1           2
# 232921      GPAM  CX3CR1     2    1    1    1      1           2
# 232922      GPAM RAPGEF5     2    1    1    1      1           2
# 232923      GPAM    CPN2     2    1    1    1      1           2
# 232924      GPAM    PROC     2    1    1    1      1           2
Comparacion_regulones[Comparacion_regulones$Regulador=="IL24",]
# Regulador Blanco Basal Her2 LumA LumB Normal Combinacion
# 1170956      IL24 POU2F3     1    1    1    1     11          11
# 1170957      IL24   RARG     1    1    1    1     11          11
# 1170958      IL24   PAX6     1    1    1    1     11          11
# 1170959      IL24  PRKCB     1    1    1    1     11          11
# 1170960      IL24  PTPRJ     1    1    1    1     11          11
# 1170961      IL24   SBF2     1    1    1    1     11          11
# 1170962      IL24   SS18     1    1    1    1     11          11
Comparacion_regulones[Comparacion_regulones$Regulador=="NXPH3",]
# Regulador Blanco Basal Her2 LumA LumB Normal Combinacion
# 478183     NXPH3   THRA     1    3    1    1      1           3
# 478184     NXPH3   PAK6     1    3    1    1      1           3
# 478185     NXPH3  ERBB2     1    3    1    1      1           3
# 478186     NXPH3 GPR179     1    3    1    1      1           3
# 478187     NXPH3   GRB7     1    3    1    1      1           3
Comparacion_regulones[Comparacion_regulones$Regulador=="PKM",]
# Regulador Blanco Basal Her2 LumA LumB Normal Combinacion
# 490974       PKM   COMP     1    3    1    1      1           3
# 490975       PKM   DKK2     1    3    1    1      1           3
# 490976       PKM  FGF14     1    3    1    1      1           3
# 490977       PKM    LBR     1    3    1    1      1           3
# 490978       PKM  NALCN     1    3    1    1      1           3
# 490979       PKM SLC6A6     1    3    1    1      1           3
# 490980       PKM  DCAF6     1    3    1    1      1           3
# 490981       PKM   MPC2     1    3    1    1      1           3
# 490982       PKM  PTPRD     1    3    1    1      1           3
# 490983       PKM RNF152     1    3    1    1      1           3
Comparacion_regulones[Comparacion_regulones$Regulador=="SEBOX",]
# Regulador  Blanco Basal Her2 LumA LumB Normal Combinacion
# 493126     SEBOX    MSX2     1    3    1    1      1           3
# 493127     SEBOX GUCY1B2     1    3    1    1      1           3
Comparacion_regulones[Comparacion_regulones$Regulador=="TCEAL2",]
# Regulador  Blanco Basal Her2 LumA LumB Normal Combinacion
# 729319    TCEAL2    GFI1     1    1    5    1      1           5
# 729320    TCEAL2  CLEC4E     1    1    5    1      1           5
# 729321    TCEAL2 IL18RAP     1    1    5    1      1           5
# 729322    TCEAL2    XCL1     1    1    5    1      1           5
Comparacion_regulones[Comparacion_regulones$Regulador=="ZNF763",]
# Regulador   Blanco Basal Her2 LumA LumB Normal Combinacion
# 705897    ZNF763   AKAP11     1    1    5    1      1           5
# 705898    ZNF763    CKAP2     1    1    5    1      1           5
# 705899    ZNF763   FNDC3A     1    1    5    1      1           5
# 705900    ZNF763   GTF2F2     1    1    5    1      1           5
# 705901    ZNF763    LACC1     1    1    5    1      1           5
# 705902    ZNF763    SUGT1     1    1    5    1      1           5
# 705903    ZNF763  TSC22D1     1    1    5    1      1           5
# 705904    ZNF763   NHLRC3     1    1    5    1      1           5
# 705905    ZNF763  PROSER1     1    1    5    1      1           5
# 705906    ZNF763    SPG20     1    1    5    1      1           5
# 705907    ZNF763   CAB39L     1    1    5    1      1           5
# 705908    ZNF763    INTS6     1    1    5    1      1           5
# 705909    ZNF763    NAA16     1    1    5    1      1           5
# 705910    ZNF763 RABGAP1L     1    1    5    1      1           5
# 705911    ZNF763   PRRC2C     1    1    5    1      1           5
# 705912    ZNF763     STX6     1    1    5    1      1           5


####################################################################
## Distribucion de tamaños de regulones por cada subtipo

tam_regulon_Basal <- table(Comparacion_regulones$Regulador[Comparacion_regulones$Basal==2])

tam_regulon_Her2 <- table(Comparacion_regulones$Regulador[Comparacion_regulones$Her2==3])

tam_regulon_LumA <- table(Comparacion_regulones$Regulador[Comparacion_regulones$LumA==5])

tam_regulon_LumB <- table(Comparacion_regulones$Regulador[Comparacion_regulones$LumB==7])

tam_regulon_Normal <- table(Comparacion_regulones$Regulador[Comparacion_regulones$Normal==11])



hist(tam_regulon_Basal,breaks = 100)
hist(tam_regulon_Her2,breaks = 100)
hist(tam_regulon_LumA,breaks = 100)
hist(tam_regulon_LumB,breaks = 100)
hist(tam_regulon_Normal,breaks = 100)

##############################################################################
## Analisis de redundancia de los reguladores


# Redundacia, definida como la proporcion de blancos que solo son regulados por el reguldor.  R= blancos del regulador - blancos no exclusivos / blsncos del regulador. 
# La métrica se encuentra entee 0 y 1, siendo 1 completamente redundante.
# Comparar la redundacia de cada regulador con su numero de blancos.

library(igraph)


str(Regulon_inmunes)
names(Regulon_inmunes)

str(Regulons)


data.frame(Regulator=as.character(SIFs_colapsados_clean[[1]]$V1),Target=as.character(SIFs_colapsados_clean[[1]]$V2),stringsAsFactors = FALSE)
regulon_Basal <- graph_from_data_frame(data.frame(Regulator=as.character(SIFs_colapsados_clean[[1]]$V1),Target=as.character(SIFs_colapsados_clean[[1]]$V2),stringsAsFactors = FALSE),
                      directed = TRUE  
  )

degree(regulon_Basal,mode="all")
degree(regulon_Basal,mode="out")

indegrees <-degree(regulon_Basal,mode="in")

r=1
t=1
n=1

#Lista para almacenar los resultados
TF_redundancies <- list()

for(n in 1:5){

  print(names(SIFs_colapsados_clean)[n]) # indicador de avance

  # Obtenemos los grados de entrada de todos los target de cada red
  data.frame(Regulator=as.character(SIFs_colapsados_clean[[n]]$V1),Target=as.character(SIFs_colapsados_clean[[n]]$V2),stringsAsFactors = FALSE)
  regulon_ <- graph_from_data_frame(data.frame(Regulator=as.character(SIFs_colapsados_clean[[n]]$V1),Target=as.character(SIFs_colapsados_clean[[n]]$V2),stringsAsFactors = FALSE),
                                         directed = TRUE  
  )
  
  indegrees <-degree(regulon_,mode="in")
  
  
    # obtenemos la lista de reguladores
tflist  <-  unique(as.character(SIFs_colapsados_clean[[n]]$V1))
# data.frame para almacenar los resultados
tf_redundancies <- data.frame(TF=tflist,targets=numeric(length(tflist)),exclusive_targets=numeric(length(tflist)),redundancy=numeric(length(tflist)))

for(r in 1:length(tflist)){ # loop reguladores
# obtenemos la lista de blancos de cada regulador
  targetlist <- as.character(SIFs_colapsados_clean[[n]]$V2)[as.character(SIFs_colapsados_clean[[n]]$V1)==tflist[r]]
exclusives = 0

for(t in 1:length(targetlist)){ #loop targets
exclusivestatus <-  ifelse(indegrees[which(names(indegrees) == targetlist[t])] == 1,1,0)
exclusives <- exclusives + exclusivestatus[1]   
print(exclusives)
 } # cierra loop targets
tf_redundancies$targets[r] <- length(targetlist)
tf_redundancies$exclusive_targets[r] <- exclusives
tf_redundancies$redundancy[r] <- (tf_redundancies$targets[r] - tf_redundancies$exclusive_targets[r])/tf_redundancies$targets[r]
print(paste(r/length(tflist)*100,"%"))
} #cierra loop reguladores

TF_redundancies[[n]] <- tf_redundancies 
names(TF_redundancies)[n]<-names(SIFs_colapsados_clean)[n]
}


pdf("Results/paso_08_Reguladores_Transcripcionales_MARINa_iRegulon/Redundancia_blancos.pdf")
plot(x=TF_redundancies[[1]]$redundancy,y=log10(TF_redundancies[[1]]$targets),main = names(TF_redundancies)[1])
plot(x=TF_redundancies[[2]]$redundancy,y=log10(TF_redundancies[[2]]$targets),main = names(TF_redundancies)[2])
plot(x=TF_redundancies[[3]]$redundancy,y=log10(TF_redundancies[[3]]$targets),main = names(TF_redundancies)[3])
plot(x=TF_redundancies[[4]]$redundancy,y=log10(TF_redundancies[[4]]$targets),main = names(TF_redundancies)[4])
plot(x=TF_redundancies[[5]]$redundancy,y=log10(TF_redundancies[[5]]$targets),main = names(TF_redundancies)[5])
dev.off()

TF_redundancies[[1]][TF_redundancies[[1]]$redundancy < 0.95,]
# TF targets exclusive_targets redundancy
# 456  LIN28B      65                 4  0.9384615
# 488  ZBTB34      11                 4  0.6363636
# 593  ZNF883      16                 2  0.8750000
# 868    LHX9      19                 1  0.9473684
# 907   CSTF2      11                 1  0.9090909
# 1243   ALX1      16                 1  0.9375000

TF_redundancies[[2]][TF_redundancies[[2]]$redundancy < 0.95,]
# TF targets exclusive_targets redundancy
# 450    ZNF257      33                 2  0.9393939
# 521   ZNF702P      13                 1  0.9230769
# 576  TMSB4XP8      17                 1  0.9411765
# 1537   ZNF788      10                 1  0.9000000

TF_redundancies[[3]][TF_redundancies[[3]]$redundancy < 0.95,]
# TF targets exclusive_targets redundancy
# 110    ZFP62     131                 7  0.9465649
# 213    ZNF17      34                 2  0.9411765
# 308   ZBTB34      13                 1  0.9230769
# 607   ZNF563      59                 3  0.9491525
# 620   BCL11B      16                 1  0.9375000
# 745   ZNF343      75                 4  0.9466667
# 820    CNOT3      17                 1  0.9411765
# 899   CREBZF      17                 1  0.9411765
# 936     ZIM3      18                 1  0.9444444
# 1092   ZNF12      14                 1  0.9285714
# 1100 ZSCAN31      11                 1  0.9090909

TF_redundancies[[4]][TF_redundancies[[4]]$redundancy < 0.95,]
# TF targets exclusive_targets redundancy
# 461   ZNF416      19                 1  0.9473684
# 762     HSF5      14                 1  0.9285714
# 954   ZNF546      12                 2  0.8333333
# 1047   TERF1      12                 1  0.9166667
# 1068  CSNK2B      19                 1  0.9473684
# 1117 ZSCAN5A       9                 1  0.8888889

TF_redundancies[[5]][TF_redundancies[[5]]$redundancy < 0.95,]
# TF targets exclusive_targets redundancy
# 127   ZNF66     133                27  0.7969925
# 146  NKX6-1     115                13  0.8869565
# 1478 ZNF806      18                 1  0.9444444



