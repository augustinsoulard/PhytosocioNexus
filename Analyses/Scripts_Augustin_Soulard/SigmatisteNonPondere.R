##################-----PARAMETRE-----#########################
nb_releve = 3

#############################################################






if(!require("sf")){install.packages("sf")} ; library("sf")
if(!require("foreign")){install.packages("foreign")} ; library("foreign") # Pour read.dbf
if(!require("RVAideMemoire")){install.packages("RVAideMemoire")} ; library("RVAideMemoire")
if(!require("vegan")){install.packages("vegan")} ; library("vegan")
if(!require("reshape2")){install.packages("reshape2")} ; library("reshape2")
source("function/taxabase.R")
source("function/habref_function.R")

#Chargement des données depuis Biodivercity
flore = dbGetQuery(con, "SELECT * FROM donnees.flore
WHERE projet IN ('PNCal') AND releve ~ '^R(1[0-7]|[1-9])$';
;")
rp = dbGetQuery(con, "SELECT * FROM donnees.releves_phytosociologiques;")

#Chargement des données de références

PVF2 = taxon_hab(28)
PVF2$CD_NOM = updatetaxa(PVF2$CD_NOM)

PVF2 = PVF2[!match(PVF2$CD_NOM,FloreRELEVE$CD_NOM,nomatch=0)==0,]


PVF2 = PVF2 %>% filter(!is.na(CD_NOM))
PVF2 = left_join(PVF2,TAXREFv17tojoin,by='CD_NOM')

#Comptage des occurences des syntaxons en correspondance aux espèces
head(PVF2 %>% count(LB_HAB_FR,sort = T),20)


PVF2contingence = PVF2 %>% select(LB_HAB_FR,LB_NOM) %>% mutate(P = 1) %>%
  pivot_wider(names_from = LB_HAB_FR, values_from = P,values_fill = list(P = 0)) # Code parfois disfonctionnel
# Code de remplacement potentiel
# PVF2contingence <- PVF2 %>%
#   select(LB_HAB_FR, LB_NOM, CD_NOM) %>%
#   mutate(P = 1) %>%
#   dcast(LB_NOM ~ LB_HAB_FR, value.var = "P", fun.aggregate = sum, fill = 0)

### pondérer les valeurs de recouvrement ####
flore$pondvalue = str_replace_all(flore$Releve_Recouvrement,c("5"="10","4"="8","3"="6","2"="4","1"="1","\\+"="0.5","r"="0.25","i"="0.125"))
flore$pondvalue = as.numeric(flore$pondvalue)
#Préparation du tableau de contingence
tabContingence = flore %>% pivot_wider(names_from = releve, values_from = Releve_Recouvrement, values_fill = list(P = 0))
tabContingence = full_join(tabContingence,PVF2contingence,by=c("lb_nom"="LB_NOM"))
#Suppression des lb_nom avec NA
tabContingence = tabContingence %>% filter(!is.na(lb_nom))

#Nommer les lignes par les noms d'espèces
row.names(tabContingence) = tabContingence$lb_nom

## Tableau de contingence des relevés avec pondération####
### Sans strate ####
tabContingence <- flore %>%
  group_by(Nom, releve) %>%
  summarise(pondvalue = sum(pondvalue, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = releve,
    values_from = pondvalue,
    values_fill = 0
  )


#Création du tableau DEDOU
DEDOU = tabContingence[,2:ncol(tabContingence)]*-1+1
row.names(DEDOU) = tabContingence$lb_nom
colnames(DEDOU) = paste0('zz',colnames(DEDOU))
DEDOU = cbind(tabContingence,DEDOU)

#Suppression des colonnes avec des NA
rows_with_na <- rowSums(is.na(DEDOU)) > 0
DEDOU <- DEDOU[!rows_with_na,]

#Réalisation de l'AFC sur DEDOU
AFC<-cca(t(DEDOU[,3:ncol(DEDOU)]))
summaryAFC = summary(AFC)
summaryAFC

MVA.synt(AFC)
stressplot(AFC)

# Extraire les scores des espèces et des sites
species_scores <- as.data.frame(scores(AFC, display = "species"))
sites_scores <- as.data.frame(scores(AFC, display = "sites"))

# Ajouter les noms des espèces et des sites aux dataframes
species_scores$species <- rownames(species_scores)
sites_scores$sites <- rownames(sites_scores)

# Filtrer les relevés (sites) qui ne commencent pas par 'n'
sites_scores_filtered <- sites_scores %>%
  filter(!grepl("^zz", rownames(sites_scores))|grepl("^R[1-3]", rownames(sites_scores)))

# Visualiser avec ggplot2 et ggrepel
if(!require("ggrepel")){install.packages("ggrepel")} ; library("ggrepel")

ggplot() +
  geom_point(data = species_scores, aes(x = CA1, y = CA2), color = "blue") +
  geom_point(data = sites_scores_filtered, aes(x = CA1, y = CA2), color = "red") +
  #geom_text_repel(data = species_scores, aes(x = CA1, y = CA2, label = species), color = "blue", max.overlaps = 20) +
  geom_text_repel(data = sites_scores_filtered, aes(x = CA1, y = CA2, label = sites), color = "red", max.overlaps = 20) +
  labs(title = "AFC des données FloreRELEVE", x = "CA1", y = "CA2")
ggsave('AFC.jpg',  width = 3000,height =2000,units='px')

if(!require("xlsx")){install.packages("xlsx")} ; library("xlsx")
write.xlsx(tabContingence,"tabContingence.xlsx")

#Matrice de distance pour vérifier les similitudes
mat_dist = dist(sites_scores)
mat_dist_matrix <- as.matrix(mat_dist)

# Afficher les associations similaires aux relevées
head(sort(mat_dist_matrix[,"R1"]),20)
head(sort(mat_dist_matrix[,2]),20)
head(sort(mat_dist_matrix[,"R3"]),20)
     