# Chargement des librairies ####
if(!require("sf")){install.packages("sf")} ; library("sf")
if(!require("tidyverse")){install.packages("tidyverse")} ; library("tidyverse")
if(!require("foreign")){install.packages("foreign")} ; library("foreign") # Pour read.dbf
if(!require("RVAideMemoire")){install.packages("RVAideMemoire")} ; library("RVAideMemoire")
if(!require("vegan")){install.packages("vegan")} ; library("vegan")
if(!require("reshape2")){install.packages("reshape2")} ; library("reshape2")


# Charement des données ####
point_silene_releve = dbGetQuery(con, "
SELECT *
  FROM bibliotaxa.point_silene t
WHERE EXISTS (
  SELECT 1
  FROM bibliotaxa.point_silene t2
  WHERE 
  t2.id <> t.id -- s'assurer que ce n’est pas lui-même
        AND DATE(t.date_debut) = DATE(t2.date_debut) -- comparer les dates au jour près
        AND t2.observateu = t.observateu
        AND t2.geom = t.geom
);
")


# Créer un ID pour les relevés ####
point_silene_releve = point_silene_releve %>%
  group_by(date_debut, observateu, geom) %>%
  mutate(groupe_id = paste0('S',cur_group_id())) %>%
  ungroup()

# Garder que les relevés avec au moins 3 espèces ####
point_silene_releve <- point_silene_releve %>%
  group_by(groupe_id) %>%
  filter(n() != 2) %>%
  ungroup()

## Tableau de contingence des relevés avec pondération####
### pondérer les valeurs de recouvrement ####
point_silene_releve$pondvalue = 1
### Sans strate ####
tabContingence <- point_silene_releve %>%
  group_by(nom_valide,groupe_id) %>%
  summarise(pondvalue = sum(pondvalue, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = groupe_id,
    values_from = pondvalue,
    values_fill = 0
  ) %>% column_to_rownames("nom_valide") %>%
  t() %>%
  as.data.frame()

#Réalisation de l'AFC sur DEDOU
AFC<-cca(tabContingence)
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

