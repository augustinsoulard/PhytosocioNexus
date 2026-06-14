# Chargement des packages ####
if(!require("foreign")){install.packages("foreign")} ; library("foreign") # Pour read.dbf
if(!require("RVAideMemoire")){install.packages("RVAideMemoire")} ; library("RVAideMemoire")
if(!require("vegan")){install.packages("vegan")} ; library("vegan")
if(!require("reshape2")){install.packages("reshape2")} ; library("reshape2")
if(!require("ggrepel")){install.packages("ggrepel")} ; library("ggrepel")
if(!require("factoextra")){install.packages("factoextra")} ; library("factoextra")
if(!require("indicspecies")){install.packages("indicspecies")} ; library("indicspecies") #multipatt
if(!require("ggdendro")){install.packages("ggdendro")} ; library("ggdendro")
if(!require("tidyverse")){install.packages("tidyverse")} ; library("tidyverse")


# Chargement des données depuis BiodiversitySQL ####
flore = dbGetQuery(con, "SELECT * FROM donnees.flore
WHERE projet IN ('PNCal') AND releve ~ '^R(1[0-7]|[1-9])$';
;")
rp = dbGetQuery(con, "SELECT * FROM donnees.releves_phytosociologiques;")



## Tableau de contingence des relevés avec pondération####
### pondérer les valeurs de recouvrement ####
flore$pondvalue = str_replace_all(flore$Releve_Recouvrement,c("5"="10","4"="8","3"="6","2"="4","1"="1","\\+"="0.5","r"="0.25","i"="0.125"))
flore$pondvalue = as.numeric(flore$pondvalue)
### Sans strate ####
tabContingence <- flore %>%
  group_by(Nom, releve) %>%
  summarise(pondvalue = sum(pondvalue, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = releve,
    values_from = pondvalue,
    values_fill = 0
  ) %>% column_to_rownames("Nom") %>%
  t() %>%
  as.data.frame()

#Réalisation de l'AFC sur DEDOU
AFC<-cca(tabContingence[,2:ncol(tabContingence)])
summaryAFC = summary(AFC)
summaryAFC

MVA.scoreplot(AFC)
MVA.synt(AFC)
stressplot(AFC)

# Extraire les scores des espèces et des sites
species_scores <- as.data.frame(scores(AFC, display = "species"))
sites_scores <- as.data.frame(scores(AFC, display = "sites"))

# Ajouter les noms des espèces et des sites aux dataframes
species_scores$species <- rownames(species_scores)
sites_scores$sites <- rownames(sites_scores)

# Visualiser avec ggplot2 et ggrepel

ggplot() +
  geom_point(data = species_scores, aes(x = CA1, y = CA2), color = "blue") +
  geom_point(data = sites_scores_filtered, aes(x = CA1, y = CA2), color = "red") +
  #geom_text_repel(data = species_scores, aes(x = CA1, y = CA2, label = species), color = "blue", max.overlaps = 20) +
  geom_text_repel(data = sites_scores_filtered, aes(x = CA1, y = CA2, label = sites), color = "red", max.overlaps = 20) +
  labs(title = "AFC des données FloreRELEVE", x = "CA1", y = "CA2")
ggsave('AFC.jpg',  width = 3000,height =2000,units='px')

#Matrice de distance pour vérifier les similitudes
mat_dist = dist(sites_scores)
mat_dist_matrix <- as.matrix(mat_dist)

# Afficher les associations similaires aux relevées
head(sort(mat_dist_matrix[,"R2"]),20)
head(sort(mat_dist_matrix[,2]),20)
head(sort(mat_dist_matrix[,"R3"]),20)


## Similitude par espèces 
#Matrice de distance pour vérifier les similitudes
mat_dist_sp = dist(species_scores)
mat_dist_sp_matrix <- as.matrix(mat_dist_sp)

# Afficher les espèces similaires
diff_esp = function(esp,nb = 30){
  valeurs <- head(sort(mat_dist_sp_matrix[,esp]), nb)
  
  # Crée un data frame pour ggplot
  df_plot <- data.frame(
    espece = names(valeurs),
    distance = as.numeric(valeurs)
  )
  
  # Ordre des espèces = celui des valeurs (trié)
  df_plot$espece <- factor(df_plot$espece, levels = df_plot$espece)
  
  # Plot ggplot
  ggplot(df_plot, aes(x = espece, y = distance,fill = distance)) +
    geom_col() +
    scale_fill_gradient(
      low = "red", 
      high = "green"
    ) +
    labs(
      title = paste0("Distances à ",esp),
      x = "Espèce",
      y = "Distance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # texte vertical
      plot.title = element_text(hjust = 0.5)
    )
}
diff_esp("Daucus carota subsp. commutatus")


# Analyse CAH ####

# 1. Matrice de distance (euclidienne par défaut)
dist_mat <- dist(tabContingence, method = "euclidean")

# 2. CAH avec méthode de liaison (ex : Ward)
cah <- hclust(dist_mat, method = "ward.D2")

# 3. Visualiser le dendrogramme
plot(cah, main = "Dendrogramme de la CAH", xlab = "", sub = "", cex = 0.8)

# 4. (Optionnel) Couper en k groupes
k <- 7  # nombre de groupes souhaité
groupes <- cutree(cah, k = k)

# Afficher les groupes par releve
sort(groupes)

# renommer les groupes
groupes<- factor(groupes,
                 levels = 1:7,
                 labels = c("Formations à Sonchus asper rudéralisées et nitrophiles","Formations littorals sur sol rocheux","Formations à Lobularia", "Pelouses psammophiles",
                            "Formations typiques à Pistachia lentiscus","Rosmarinion","Formations dégradées de Pistachia lentiscus"))

# Convertir l'objet hclust en objet dendrogramme ggplot-compatible
dend_data <- dendro_data(as.dendrogram(cah), type = "rectangle")

# Les labels (feuilles) sont ici :
labels_df <- dend_data$labels  # contient x et label

# Associer les groupes aux labels
labels_df <- labels_df %>%
  mutate(
    groupe_num = groupes[label],
    groupe_nom = levels(groupes)[groupes[label]]
  )

# Calculer position moyenne par groupe
groupe_positions <- labels_df %>%
  group_by(groupe_num, groupe_nom) %>%
  summarise(x = median(x), .groups = "drop")
p <- fviz_dend(cah, k = k,
               rect = TRUE, rect_fill = TRUE, rect_border = "gray30",
               show_labels = TRUE, main = "Dendrogramme CAH des communautés végétales")

# Ajouter les labels de groupes automatiquement
p + geom_text(
  data = groupe_positions,
  aes(x = x-0.2, y = 0, label = groupe_nom),
  angle = 90, hjust = 0, size = 3
)


# Individue valeur
indval_res <- multipatt(tabContingence, groupes, control = how(nperm=999))
summary(indval_res)


# Résumé sous forme de data frame
indval_df <- as.data.frame(indval_res$sign)

# Ajouter les noms d’espèces
indval_df$espece <- rownames(indval_df)

# Filtrer les espèces avec p < 0.05
indval_sig <- indval_df %>% filter(p.value <= 0.05)


# Espèce associée à un seul groupe
indval_long <- indval_sig %>%
  select(starts_with("s."), stat, p.value, espece) %>%
  pivot_longer(starts_with("s."), names_to = "groupe", values_to = "appartenance") %>%
  filter(appartenance == 1)

# Barplot
ggplot(indval_long, aes(x = reorder(espece, stat), y = stat, fill = groupe)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Espèces indicatrices significatives (IndVal)",
    x = "Espèce",
    y = "Valeur IndVal",
    fill = "Groupe"
  ) +
  theme_minimal()

