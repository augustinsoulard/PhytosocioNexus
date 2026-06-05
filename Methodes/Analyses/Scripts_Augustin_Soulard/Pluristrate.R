# On retirer les espèces sans code CATMINAT
JOIN_DATA_CALL = JOIN_DATA[!is.na(JOIN_DATA$code_CATMINAT),]

DATARELEVE = aggregate(JOIN_DATA,list(JOIN_DATA$RELEVE,JOIN_DATA$STRATE),mean)
DATARELEVE = DATARELEVE[,c(1,2)]
colnames(DATARELEVE) = c('RELEVE','STRATE')
DATARELEVE$CATMINAT = ""
#Boucle qui tourne pour chaque releve
for(i in 1:nrow(DATARELEVE)){
  cat("RELEVE : ",DATARELEVE$RELEVE[i],"\n")
  JOIN_DATA_CALL = JOIN_DATA[!is.na(JOIN_DATA$code_CATMINAT),]
  JOIN_DATA_CALL = JOIN_DATA_CALL[JOIN_DATA_CALL$RELEVE == DATARELEVE$RELEVE[i] 
                                   & JOIN_DATA_CALL$STRATE == DATARELEVE$STRATE[i],]
  CATMINAT_SEP = separate(JOIN_DATA_CALL, code_CATMINAT, into = c("code1", "code2", "code3", "code4", "code5","code6","code7"), sep = "/|\\.")
  JOIN_DATA_CALL = left_join(CATMINAT_SEP,JOIN_DATA_CALL,by="CD_NOM",keep = FALSE)[,-c(14:16)]
    
  # On analyse les code CATMINAT dans l'ordre d'importance
    for(j in 1:6){
      cat("___VERIFICATION CODE : ",j,"\n")
      col = paste0("code",j)
      if(all(is.na(JOIN_DATA_CALL[col]))){break}
      AGG = aggregate(pondvalue.x~get(paste0("code",j)),data=JOIN_DATA_CALL,sum)
      sort_value = sort(AGG$pondvalue.x,decreasing = T)
      if(sort_value[1]<sort_value[2]*1.2 & length(sort_value)>1){
        cat("______PROXIMITE OU EGALITE DES VALEURS AU NIVEAU",j, "\n")
        break}
      
      #Attribution du code CATMINAT selon l'avancement
      CATMIN = as.character(AGG[AGG$pondvalue.x==max(AGG$pondvalue.x),][1])
      if(j == 1){
        sep1 = "/" ; sep2 = ""
      } else if(j == 2){
        sep1 = "." ; sep2 = ""
      } else if(j == 3){
        sep1 = "" ; sep2 = ""
      } else{
        sep1 = "" ; sep2 = "."
      }
      DATARELEVE$CATMINAT[i] = paste0(DATARELEVE$CATMINAT[i],sep2,CATMIN,sep1)
      JOIN_DATA_CALL = JOIN_DATA_CALL[JOIN_DATA_CALL[col]== CATMIN,]

    }
}


# Fonction pour enlever ".0" à la fin des
enlever_point_zero <- function(x) {
  if (str_sub(x, nchar(x)-1, nchar(x)) == ".0") {
    x <- str_sub(x, 1, nchar(x)-2)
  }
  return(x)
}

DATARELEVE$CATMINAT = sapply(DATARELEVE$CATMINAT, enlever_point_zero)

# Une fois un code attribu? ? chaque relev?, on y joint les informations de baseveg

baseflor <- read.csv("baseflor.csv",sep=";",fileEncoding = "latin1")
DATARELEVE_JOIN = left_join(DATARELEVE,baseflor,by=c("CATMINAT"="code_CATMINAT"))
DATARELEVE_JOIN = aggregate(DATARELEVE_JOIN,list(DATARELEVE_JOIN$RELEVE,DATARELEVE_JOIN$STRATE),unique)

# Ajout des nomenclature de phytosociologie a JOIN_DATA
JOIN_DATA$IDKEY = paste0(JOIN_DATA$STRATE,JOIN_DATA$CD_NOM,JOIN_DATA$RELEVE)
JOIN_DATA_EXP = inner_join(JOIN_DATA,baseflor,by=c("code_CATMINAT"="code_CATMINAT"))
JOIN_DATA_EXP = JOIN_DATA_EXP[!duplicated(JOIN_DATA_EXP$IDKEY),]
JOIN_DATA_EXP = select(JOIN_DATA_EXP,RELEVE,STRATE,CD_NOM,NOM_VALIDE,Nom.scientifique, value,pondvalue,code_CATMINAT,CARACTERISATION_ECOLOGIQUE_.HABITAT_OPTIMAL.,INDICATION_PHYTOSOCIOLOGIQUE_CARACTERISTIQUE)
# On retire les synonymes
# DATARELEVE_JOIN = DATARELEVE_JOIN[!DATARELEVE_JOIN$NIVEAU %in% c("syn =","syn = pp","syn compl inval pp","syn compl pp","syn pp"),]

# On retirer les colonnes inutiles
DATARELEVE_JOIN = DATARELEVE_JOIN %>% select(RELEVE,STRATE,CATMINAT,CARACTERISATION_ECOLOGIQUE_.HABITAT_OPTIMAL.,INDICATION_PHYTOSOCIOLOGIQUE_CARACTERISTIQUE)

# Enregistrement des résultats
dir.create("OUTPUT")
DATARELEVE_JOIN = DATARELEVE_JOIN %>% arrange(RELEVE,STRATE)
write.csv2(DATARELEVE_JOIN,"OUTPUT/BILAN_RELEVE_STRATE.csv",fileEncoding = "UTF-8",row.names = F)

# Export de JOIN_DATA
JOIN_DATA_EXP = JOIN_DATA_EXP %>% arrange(RELEVE ,STRATE,code_CATMINAT)
write.csv2(JOIN_DATA_EXP,"OUTPUT/BILAN_ESP_STRATE.csv",fileEncoding = "UTF-8",row.names = F)
