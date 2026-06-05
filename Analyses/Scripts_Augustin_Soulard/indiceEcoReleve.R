# Chargement de las base de donnees base flore avec CD_NOM
baseflorTAXREFv16 = read.csv("baseflorTAXREFv16.csv", sep=";")


# Jointure de baseflor a DATAOPHYTO
DATAbaseflorJOIN = left_join(DATAPHYTO,baseflorTAXREFv16,by="CD_NOM")

#Verification de la jointure et NA
if(any(is.na(DATAbaseflorJOIN$N._Nomenclatural_BDNFF)==TRUE)){
  cat("!!!ATTENTION JOINTURE INCOMPLETE !!!")
  
}

# Preparation du tableau final des releves
DATA_RELEVE_INDICE = data.frame(
  RELEVE = unique(DATAbaseflorJOIN$RELEVE),
  NUM = as.numeric(substring(unique(DATAbaseflorJOIN$RELEVE), first = 2, last = 3)),
  Humidite_edaph = NA,
  Humidite_ZH_regl = NA, # Non humide >5.5 et 5.9 < Humide
  Lumiere = NA,
  Temperature = NA,
  Humidite_atmosp = NA,
  React_sol_pH = NA,
  Niv_trophiq = NA,
  Salinite = NA,
  Texture = NA,
  Matiere_organiq = NA
)

#Boucle principale d'analyse
for(i in 1:nrow(DATA_RELEVE_INDICE)){
  cat("RELEVE : ",as.character(DATA_RELEVE_INDICE$RELEVE[i])," - ")
  DONNEE_REL_I = DATAbaseflorJOIN[DATAbaseflorJOIN$RELEVE == DATA_RELEVE_INDICE$RELEVE[i],]
  #Calcule des moyens par indices sans pondération
  DATA_RELEVE_INDICE$Humidite_edaph[i] = mean(DONNEE_REL_I$Humidité_édaphique,na.rm = TRUE)
  
  if(DATA_RELEVE_INDICE$Humidite_edaph[i]>=5.9){
    DATA_RELEVE_INDICE$Humidite_ZH_regl[i] = 'Zone humide'
  } else if(DATA_RELEVE_INDICE$Humidite_edaph[i]<=5.5){
    DATA_RELEVE_INDICE$Humidite_ZH_regl[i] = 'Zone non humide'
  } else{ DATA_RELEVE_INDICE$Humidite_ZH_regl[i] = 'Incertitude'}
  
  DATA_RELEVE_INDICE$Lumiere[i] = mean(DONNEE_REL_I$Lumière,na.rm = TRUE)
  DATA_RELEVE_INDICE$Temperature[i] = mean(DONNEE_REL_I$Température,na.rm = TRUE)
  
  DATA_RELEVE_INDICE$Humidite_atmosp[i] = mean(DONNEE_REL_I$Humidité_atmosphérique,na.rm = TRUE)
  
  DATA_RELEVE_INDICE$React_sol_pH[i] = mean(DONNEE_REL_I$Réaction_du_sol_.pH.,na.rm = TRUE)
  
  DATA_RELEVE_INDICE$Niv_trophiq[i] = mean(DONNEE_REL_I$Niveau_trophique,na.rm = TRUE)
  
  DATA_RELEVE_INDICE$Salinite[i] = mean(DONNEE_REL_I$Salinité,na.rm = TRUE)
  
  DATA_RELEVE_INDICE$Texture[i] = mean(DONNEE_REL_I$Texture,na.rm = TRUE)
  
  DATA_RELEVE_INDICE$Matiere_organiq[i] = mean(DONNEE_REL_I$Matière_organique,na.rm = TRUE)
  
    cat("Humidité édaphique : ",DATA_RELEVE_INDICE$Humidite_edaph[i],"\n")
}


# Export du tableau final des relevés
write.csv(DATA_RELEVE_INDICE,"OUTPUT/DATA_RELEVE_INDICE.csv",row.names = FALSE, fileEncoding = "UTF-8")



