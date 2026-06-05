#Chargement des données ZH
Flore_ZH <- read.csv2("data\ZH/FLore_ZH.csv",colClasses = "character")

#Jointure de la flore indicatrice zone humide
DATA_FLORE_ZH = left_join(DATAbaseflorJOIN,Flore_ZH,by="CD_NOM")

#Numérotage des relevés
DATA_FLORE_ZH$NUM = as.numeric(substring(DATA_FLORE_ZH$RELEVE, first = 2, last = 3))
# Ajout de la colonne somme cumulée
DATA_FLORE_ZH = DATA_FLORE_ZH %>% select(NUM,
                                         RELEVE, 
                                         CD_SYNUSIE, 
                                         CD_NOM, 
                                         NomComplet, 
                                         NOM_VERN,
                                         Surface,
                                         Indicatrice.ZH)



# TRI du tableau
DATA_FLORE_ZH = DATA_FLORE_ZH %>% arrange(desc(Surface)) %>% arrange(CD_SYNUSIE) %>% arrange(NUM)


# Calcul de la somme cumulée
DATA_FLORE_ZH = DATA_FLORE_ZH %>% group_by(RELEVE, CD_SYNUSIE) %>% mutate(CUMSUM = cumsum(Surface))



# Preparation du tableau final des releves
DATA_RELEVE_ZH = data.frame(
              RELEVE = unique(DATAbaseflorJOIN$RELEVE),
              NUM = as.numeric(substring(unique(DATAbaseflorJOIN$RELEVE), first = 2, last = 3)),
              ZONE_HUMIDE = NA
           )

#Boucle principale d'analyse
for(i in 1:nrow(DATA_RELEVE_ZH)){
  cat("RELEVE : ",as.character(DATA_RELEVE_ZH$RELEVE[i])," - ")
  DONNEE_REL_I = DATA_FLORE_ZH[DATA_FLORE_ZH$RELEVE == DATA_RELEVE_ZH$RELEVE[i],]
  # Surface >50%
  DONNEE_REL_I = DONNEE_REL_I %>% filter(!(CUMSUM > 50 & lag(CUMSUM) >50 & CD_SYNUSIE == lag(CD_SYNUSIE) & Surface < 20)) %>% filter(Surface>=5)
  if(nrow(DONNEE_REL_I[DONNEE_REL_I$Indicatrice.ZH=="Oui",]) >= nrow(DONNEE_REL_I)/2 ){
    DATA_RELEVE_ZH$ZONE_HUMIDE[i] = "Oui"
  } else(DATA_RELEVE_ZH$ZONE_HUMIDE[i] = "-")
 cat("Zone humide : ",DATA_RELEVE_ZH$ZONE_HUMIDE[i],"\n")
}

# Export du tableau final des relevés
write.csv(DATA_RELEVE_ZH,"OUTPUT/DATA_RELEVE_ZH.csv",row.names = FALSE, fileEncoding = "UTF-8")



# Export du tableau final d'espèce
DATA_FLORE_ZH_EXP = DATA_FLORE_ZH %>% select(RELEVE, 
                                         CD_SYNUSIE, 
                                         CD_NOM, 
                                         NomComplet, 
                                         NOM_VERN,
                                         Surface,
                                         CUMSUM,
                                         Indicatrice.ZH)
colnames(DATA_FLORE_ZH_EXP) = c("Relevé","Strates","CD_NOM","Nom scientifique", "Nom vernaculaire","Recouvrement (%)", "Recouvrement cumulé (%)", "Espèce indicatrice de zone humide")
write.csv(DATA_FLORE_ZH_EXP,"OUTPUT/DATA_FLORE_ZH.csv",row.names = FALSE, fileEncoding = "UTF-8")
