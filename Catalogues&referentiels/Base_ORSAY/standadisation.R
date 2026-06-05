if(!require("tidyverse")){install.packages("tidyverse")} ; library("tidyverse")
if(!require("reshape2")){install.packages("reshape2")} ; library("reshape2")
if(!require("readxl")){install.packages("readxl")} ; library("readxl")
all_results <- list()

dossiers <- paste0("Tableaux des relevés/t000", 1:5)
fichiers <- paste0("resultat_", 1:10, ".csv")


for (dossier in dossiers) {
  for (fichier in fichiers) {
    cat(dossier,"Fichier",fichier,"\n")
    # Construire le chemin du fichier
    chemin_fichier <- file.path(dossier, fichier)
    
    # Vérifier si le fichier existe
    if (file.exists(chemin_fichier)) {
      # Lire le fichier CSV
      data <- read.csv2(chemin_fichier)
      
      # Transformer les données en utilisant melt
      data_melted <- melt(data, id.vars = "X", variable.name = "Site", value.name = "Abondance")
      
      # Filtrer les lignes avec des valeurs vides ou '.'
      data_melted <- data_melted[!(data_melted$Abondance == '' | data_melted$Abondance == '.'), ]
      
      # Ajouter le résultat à la liste
      all_results[[length(all_results) + 1]] <- data_melted
    }
  }
}

# Combiner tous les résultats en un seul data frame
final_result <- bind_rows(all_results)
colnames(final_result) = c("espece","releve","indicebb")
final_result$releve = toupper(final_result$releve)

typevege = read.csv2("Types de vegetation.csv",h=T,encoding = "latin1")
colnames(typevege) = c("releve","typevegetation")
typevege$releve = toupper(typevege$releve)

final_result = left_join(final_result,typevege,by="releve")


Loc<- read_excel("Localisation des relevés.xlsx")
Loc$COUNTRY = NULL
colnames(Loc) = c("releve","localisation")
Loc$releve = toupper(Loc$releve)

final_result = left_join(final_result,Loc,by="releve")

write.csv(final_result,"releveORSAY.csv",row.names = F,fileEncoding = 'UTF-8')


#Exemple filtre


sort(unique(final_result$localiation))
sort(unique(final_result$typevegetation))
OSRSAY_sud_pelouse_seche = final_result %>% 
  filter(localiation %in% c("Provence", "Provence calcaire", "Provence cristalline") &
           startsWith(typevegetation, "Pelouse"))


sort(unique(OSRSAY_sud_pelouse_seche$releve))




############Pareil pour head


df = read.csv2("Tableaux des relevés/t0001/header_1.csv",h=T)
colnames(df) = c("X","releve","info")
df$Index <- ave(df$info, df$releve, FUN = seq_along)

df = dcast(df, releve ~ Index, value.var = "info")


all_results <- list()

dossiers <- paste0("Tableaux des relevés/t000", 1:5)
fichiers <- paste0("header_", 1:10, ".csv")


for (dossier in dossiers) {
  for (fichier in fichiers) {
    cat(dossier,"Fichier",fichier,"\n")
    # Construire le chemin du fichier
    chemin_fichier <- file.path(dossier, fichier)
    
    # Vérifier si le fichier existe
    if (file.exists(chemin_fichier)) {
      # Lire le fichier CSV
      data <- read.csv2(chemin_fichier)
      
      # Transformer les données en utilisant melt
      data_melted <- melt(data, id.vars = "X", variable.name = "Site", value.name = "Abondance")
      
      # Filtrer les lignes avec des valeurs vides ou '.'
      data_melted <- data_melted[!(data_melted$Abondance == '' | data_melted$Abondance == '.'), ]
      
      # Ajouter le résultat à la liste
      all_results[[length(all_results) + 1]] <- data_melted
    }
  }
}

# Combiner tous les résultats en un seul data frame
final_result <- bind_rows(all_results)