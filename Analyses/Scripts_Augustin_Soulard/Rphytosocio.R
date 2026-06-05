# Chargement des packages
if(!require("foreign")){install.packages("foreign")} ; library("foreign")
if(!require("reshape2")){install.packages("reshape2")} ; library("reshape2")
if(!require("tidyverse")){install.packages("tidyverse")} ; library("tidyverse")
if(!require("xlsx")){install.packages("xlsx")} ; library("xlsx")

# Charger le fichier d'export des releves en recherchant un CSV
# DATAPHYTO = read.csv(choose.files())

# Ouverture direct du SHP et geopackage avec jointure
# Package d'ouverture d'un geopackage
if(!require("sf")){install.packages("sf")} ; library("sf")


#Préparation de la zone d'export
dir.create("OUTPUT")
# Ouverture des fichiers
TAXAQgis = read_sf("DataProjetFlore.gpkg", layer = "TAXAQgis")
Flore = read.dbf("FLore/Flore.dbf")

# Création du fichier DATAPHYTO
DATAPHYTO = left_join(Flore,TAXAQgis,by=c("NomComplet" = "LB_NOM"))
# Message erreur si NA
if(any(is.na(DATAPHYTO$FAMILLE) == TRUE)){
  cat("!!!ATTENTION JOINTURE DES CD NOM INCOMPLETE !!!")
}


#Choix du protocole utilisé
# sigm pour les rélevés sigmatistes
# syne pour les relevés synusiaux
# zhum pour les relevés du protocole réglementaire zones humides



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
protocole = "zhum" # ATTENTION !!!!!!!!!!! VALEUR A ADAPTER
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if(protocole == "zhum"){

DATAPHYTO = DATAPHYTO[DATAPHYTO$REL_ZH == "Oui",] 


#Calcule des indices écologiques sans pondération de recouvrement
source("R_script/indiceEcoReleve.R")
}




if(protocole == "sigm"){
  #Decroiser le tableau
  MELT_DATA <- melt(DATAPHYTO, id.vars = c("CD_NOM","NomComplet","RELEVE"), measure.vars = c("ARBO", "ARBU", "HERB", "MUCINALE"))
  
  
  # Attribution d'un poid pour chaque classe phyto
  MELT_DATA$pondvalue = str_replace_all(MELT_DATA$value,c("5"="10","4"="8","3"="6","2"="4","1"="1","\\+"="0.5","r"="0.25","i"="0.125"))
  MELT_DATA$pondvalue = as.numeric(MELT_DATA$pondvalue)
  colnames(MELT_DATA) = c("CD_NOM","NOM_VALIDE","RELEVE","STRATE","value","pondvalue")
  
  #Somme des scores par espèce
  AGG_DATA = aggregate(pondvalue~CD_NOM+NOM_VALIDE+RELEVE+STRATE+value,data=MELT_DATA,sum)
  
  # Joindre les informations de baseflor
  TAXREF_baseflor <- read.csv("TAXREF_baseflor.csv", sep=";")
  
  JOIN_DATA  = left_join(AGG_DATA,TAXREF_baseflor,by ="CD_NOM")
  
  JOIN_DATA = select(JOIN_DATA,RELEVE,STRATE,CD_NOM,NOM_VALIDE,Nom.scientifique, value,pondvalue,code_CATMINAT)
  
  
  # ANALYSE PAR STRATE !
  source("R_script/Pluristrate.R")
  
  # ANALYSE TOUTES STRATES CONFONDUES
  source("R_script/Monostrate.R")
  
  ################################################################################
  ###################CREATION DU FICHIER RMARKDOWN################################
  ################################################################################
  if(!require("knitr")){install.packages("knitr")} ; library("knitr")
  if(!require("flextable")){install.packages("flextable")} ; library("flextable")
  
  # Fonction qui enregistre les sorties de R pour 
  con <- file("Analyse_phyto.Rmd", open = "wt", encoding = "UTF-8")
  sink(con,split=T)
  cat("---
title: \"Analyse phytosociologique\"
date: \"`r Sys.Date()`\"
output:
  word_document: 
    reference_docx: word_styles_references.dotm
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
``` \n \n")
  LIST_RELEVE = levels(as.factor(DATARELEVE_JOIN$RELEVE))
  for(i in 1:length(LIST_RELEVE)){
    cat("### Relevé : ",LIST_RELEVE[i],"\n \n")
    JOIN_DATA[JOIN_DATA$RELEVE==LIST_RELEVE[i],]
    cat("Tableau des espèces : 
```{r ",paste0("ESP",LIST_RELEVE[i]),",echo=TRUE}
flextable(JOIN_DATA[JOIN_DATA$RELEVE==LIST_RELEVE[",i,"],])
``` 
\n \n")
    cat("Tableau des relevés : 
```{r", paste0("RELEVE",LIST_RELEVE[i]),",echo=TRUE}
knitr::kable(DATARELEVE_JOIN[DATARELEVE_JOIN$RELEVE==LIST_RELEVE[",i,"],])
``` 
\n \n")
    cat(paste0("![](",DATAPHYTO[DATAPHYTO$RELEVE == LIST_RELEVE[i],]$RP_PHOTO1[1],") \n \n"))
    
  }
  sink()
  close(con)
  
  #Test knit
  rmarkdown::render("Analyse_phyto.Rmd",output_format = "word_document",output_file = "Analyse_phyto.docx")
  
}

