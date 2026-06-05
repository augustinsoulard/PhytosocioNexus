# On retirer les espèces sans code CATMINAT
JOIN_DATA_CALL = JOIN_DATA[!is.na(JOIN_DATA$code_CATMINAT),]

DATARELEVE = data.frame(RELEVE = levels(factor(JOIN_DATA_CALL$RELEVE)),CATMINAT = "")
#Boucle qui tourne pour chaque releve
for(i in 1:nrow(DATARELEVE)){
  cat("RELEVE : ",DATARELEVE$RELEVE[i],"\n")
  JOIN_DATA_CALL = JOIN_DATA[!is.na(JOIN_DATA$code_CATMINAT),]
  JOIN_DATA_CALL = JOIN_DATA_CALL[JOIN_DATA_CALL$RELEVE == DATARELEVE$RELEVE[i],]
  CATMINAT_SEP = separate(JOIN_DATA_CALL, code_CATMINAT, into = c("code1", "code2", "code3", "code4", "code5","code6","code7"), sep = "/|\\.")
  JOIN_DATA_CALL = left_join(CATMINAT_SEP,JOIN_DATA_CALL,by="CD_NOM",keep = FALSE)[,-c(13:16)]
  # On analyse les code CATMINAT dans l'ordre d'importance
  for(j in 1:5){
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

# Une fois un code attribu? ? chaque relev?, on y joint les informations de baseveg

baseflor <- read.csv("baseflor.csv",sep=";",fileEncoding = "latin1")
DATARELEVE_JOIN = inner_join(DATARELEVE,baseflor,by=c("CATMINAT"="code_CATMINAT"))
DATARELEVE_JOIN = aggregate(DATARELEVE_JOIN,list(DATARELEVE_JOIN$RELEVE),unique)

# Ajout des nomenclature de phytosociologie a JOIN_DATA
JOIN_DATA$IDKEY = paste0(JOIN_DATA$STRATE,JOIN_DATA$CD_NOM,JOIN_DATA$RELEVE)
JOIN_DATA_EXP = inner_join(JOIN_DATA,baseflor,by=c("code_CATMINAT"="code_CATMINAT"))
JOIN_DATA_EXP = aggregate(JOIN_DATA_EXP,list(JOIN_DATA_EXP$IDKEY),unique)
JOIN_DATA_EXP = select(JOIN_DATA_EXP,RELEVE,CD_NOM,NOM_VALIDE,Nom.scientifique, value,pondvalue,code_CATMINAT,CARACTERISATION_ECOLOGIQUE_.HABITAT_OPTIMAL.,INDICATION_PHYTOSOCIOLOGIQUE_CARACTERISTIQUE)
# On retire les synonymes
# DATARELEVE_JOIN = DATARELEVE_JOIN[!DATARELEVE_JOIN$NIVEAU %in% c("syn =","syn = pp","syn compl inval pp","syn compl pp","syn pp"),]

# On retirer les colonnes inutiles
DATARELEVE_JOIN = DATARELEVE_JOIN %>% select(RELEVE,CATMINAT,CARACTERISATION_ECOLOGIQUE_.HABITAT_OPTIMAL.,INDICATION_PHYTOSOCIOLOGIQUE_CARACTERISTIQUE)

# Enregistrement des résultats
dir.create("OUTPUT")
write.csv(DATARELEVE_JOIN,"OUTPUT/BILAN_RELEVE.csv",fileEncoding = "UTF-8",row.names = F)

# Export de JOIN_DATA
JOIN_DATA_EXP = JOIN_DATA_EXP %>% arrange(RELEVE ,code_CATMINAT)
write.csv(JOIN_DATA_EXP,"OUTPUT/BILAN_ESP.csv",fileEncoding = "UTF-8",row.names = F)