# Introduction à R

# Installation de R et RStudio ####
# Pour commencer, assurez-vous d'avoir installé R et RStudio.
# Vous pouvez télécharger R depuis https://cran.r-project.org/
# et RStudio depuis https://posit.co/download/rstudio-desktop/.

# Installation de R et RStudio
# 1. Allez sur les sites mentionnés ci-dessus.
# 2. Téléchargez et installez R et RStudio.
# 3. Ouvrez RStudio et créez un nouveau script R.

# Manipulation de données avec R

# Charger les packages nécessaires
library(dplyr)
library(readr)
library(tidyr)

# Importation de données ####
# Utilisons un exemple de jeu de données de relevés phytosociologiques
# ?nous avons un fichier CSV nommé 'data_releve_type.csv'.
# Voici comment importer ces données :

# Exemple de code pour importer des données
data_releve <- read_csv2("data_releve_type.csv")


data_releve = data_releve %>% mutate(abondance_dominance = case_when(
  abondance_dominance == "+" ~ 0.5,
  TRUE ~ as.numeric(as.character(abondance_dominance))
))

# Filtrer les données utiles
# Exemple ne choisir que certains relevés
data_releve_filtre = data_releve %>% filter(releve %in% c("R1","R2","R5"))
# Visualisation de données

# Charger le package ggplot2
library(ggplot2)

# Création de graphiques heatmap pour visualiser les similitudes des relevés

ggplot(data_releve_filtre, aes(x = espece, y = releve, fill = abondance_dominance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +  # Ajustez la palette de couleurs selon vos préférences
  labs(x = "Espèce", y = "Relevé", fill = "Abondance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Utilisons des méthodes d'ordination pour analyser les communautés végétales.
library(vegan)

# Convertir les données en matrice appropriée pour l'analyse NMDS
library(tibble)

# Créer une colonne combinée pour strate et espece
data_releve_combined <- data_releve %>%
  mutate(combination = paste(espece, strate, sep = "_")) %>%
  group_by(releve, combination) %>%
  summarise(abondance_dominance = mean(abondance_dominance, na.rm = TRUE), .groups = 'drop')


# Préparer les données pour la table matricielle avec les données agrégées
data_releve_matrix <- data_releve_combined %>%
  pivot_wider(names_from = combination, values_from = abondance_dominance, values_fill = 0) %>%
  column_to_rownames(var = "releve") %>%  as.matrix()  # Assure une conversion en matrice# Remplacer les valeurs manquantes par 0, au cas où

data_releve_matrix[is.na(data_releve_matrix)] <- 0

# Calculer la matrice des distances de Bray-Curtis
distance_matrix <- vegdist(data_releve_matrix, method = "bray")

# Effectuer la Classification Ascendante Hiérarchique
cah_result <- hclust(distance_matrix, method = "ward.D2")


# Couper le dendrogramme pour obtenir des groupes (par exemple, 3 groupes)
num_groups <- 3
groups <- cutree(cah_result, k = num_groups)

# Convertir les groupes en facteur
groups <- factor(groups)

# Tracer le dendrogramme avec les groupes colorés
plot(cah_result, main = "Classification Ascendante Hiérarchique", xlab = "Relevés", sub = "")
rect.hclust(cah_result, k = num_groups, border = 2:4)

# Exécuter l'analyse NMDS ####
set.seed(123) # pour la reproductibilité
nmds_result <- metaMDS(data_releve_matrix, k = 2, trymax = 100, autotransform = FALSE)
stress_val <- round(nmds_result$stress, 3)

# Récupérer les scores NMDS
nmds_sites <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_sites$label <- rownames(nmds_sites)

# Assigner les groupes aux relevés NMDS
nmds_sites$Groupe <- groups[rownames(nmds_sites)]

# Visualisation avec ggplot2
ggplot(nmds_sites, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Groupe), size = 3) +
  ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 100) +
  labs(title = "Ordination NMDS",
       subtitle = paste("Stress:", stress_val),
       x = "NMDS 1", y = "NMDS 2") +
  theme_minimal() +
  coord_equal()


# Retirer des relevés
data_releve = data_releve %>% filter(!releve %in% c("R21","R8","R24"))
# Relancer maintenant après la création de relevé


# Analyses des données environnementales ####
# Charger les données environnementales
library(readxl)
env_data <- read_excel("envdata_releve.xlsx") # Remplacez par le chemin correct
env_data = env_data[order(env_data$Nom), ]
rownames(env_data) <- env_data$Nom # Assurez-vous que les lignes sont nommées par les relevés

#Retirer les relevés à retirer : 
env_data = env_data %>% filter(!Nom %in% c("R21","R8","R24","R25"))

# Ordonner les tableaux de la même manière

data_releve_matrix = data_releve_matrix[order(rownames(data_releve_matrix)), ]

all(rownames(data_releve_matrix) == rownames(env_data)) # TRUE c'est que les données sont prêtes pour la CCA. AUtrement réordonner les tables

#Filtrer les variables utilisées
env_data = env_data %>% select(Altitude, Pente,Recouvrement_herbacee,Recouvrement_arbustive, Recouvrement_arboree,
                               Hauteur_herbacee,Hauteur_arbustive,Hauteur_arboree)# Supprimer la colonne 'releve' si elle est incluse dans les données environnementales


# Exécuter la CCA
cca_result <- cca(data_releve_matrix ~ ., data = env_data)
####IMPORTANT##### les relevés des tableaux data_releve_matrix et env_data doivent ^tre les mêmes et dans le même ordre de ligne

# Extraire les scores des sites et des espèces
site_scores <- vegan::scores(cca_result, display = "sites")
species_scores <- vegan::scores(cca_result, display = "species")
biplot_scores <- vegan::scores(cca_result, display = "bp") # Variables environnementales

# Créer un dataframe pour les sites
df_sites <- data.frame(
  Site = rownames(site_scores),
  CCA1 = site_scores[, 1],
  CCA2 = site_scores[, 2]
)

# Créer un dataframe pour les espèces
df_species <- data.frame(
  Species = rownames(species_scores),
  CCA1 = species_scores[, 1],
  CCA2 = species_scores[, 2]
)

# Créer un dataframe pour les variables environnementales (biplot)
df_env <- data.frame(
  Variable = rownames(biplot_scores),
  CCA1 = biplot_scores[, 1],
  CCA2 = biplot_scores[, 2]
)

### récupérer la contribution des axes :
# Stocker le résumé de l'analyse
sommaire_cca <- summary(cca_result)

# Récupérer le tableau de la contribution des axes contraints
contribution_axes <- sommaire_cca$concont$importance

# Extraire la proportion de variance expliquée pour les axes 1 et 2
# et la convertir en pourcentage joliment formaté
cca1_percent <- round(sommaire_cca$concont$importance[2, 1] * 100, 1)
cca2_percent <- round(sommaire_cca$concont$importance[2, 2] * 100, 1)

# Créer les nouvelles étiquettes pour les axes
axe_x_label <- paste0("CCA1 (", cca1_percent, "%)")
axe_y_label <- paste0("CCA2 (", cca2_percent, "%)")

# Visualisation avec ggplot2 pour les relevés et variables environnementales
library(ggplot2)
library(ggrepel)

ggplot() +
  # --- Relevés ---
  geom_point(data = df_sites, aes(x = CCA1, y = CCA2),
             color = "black", size = 3) +
  geom_text(data = df_sites, aes(x = CCA1, y = CCA2, label = Site),
            vjust = -0.5, size = 3, color = "black") +
  
  # --- Variables environnementales ---
  geom_segment(data = df_env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text_repel(data = df_env, aes(x = CCA1, y = CCA2, label = Variable),
                  size = 3, color = "red", segment.color = "grey50") +
  
  labs(title = "CCA Biplot", x = axe_x_label, y = axe_y_label) +
  theme_minimal()

# Test de permutation
permutest(cca_result, permutations = 999)

# Summary pour voir la variance expliquée
summary(cca_result)


################################ Indice Value

library(indicspecies)

# Assurez-vous que 'groups' est un facteur
groups <- factor(groups)

# --- CORRECTION ICI ---
# Utilisez votre matrice de communauté (relevés x espèces), PAS la matrice de distance.
# Je suppose qu'elle s'appelle 'data_releve_matrix' d'après votre script précédent.
indval_res <- multipatt(data_releve_matrix, groups, 
                        func = "IndVal.g", 
                        control = how(nperm = 999))

# Extraire le tableau des espèces significatives (p-value <= 0.05 par défaut)
indval_df <- as.data.frame(indval_res$sign)

# Ajouter les noms d'espèces comme une colonne (ils sont dans les noms de lignes)
indval_df$Espèce <- rownames(indval_df)

# Renommer et sélectionner les colonnes pour un affichage propre
indval_df_filtre <- indval_df %>%
  select(Espèce, Groupe_Indicateur = index, Indice_IndVal = stat, `p-value` = p.value)

# Afficher le tableau final
print(indval_df_filtre)


### Appartenance des relevés aux groupes
# Assurez-vous que le package dplyr est chargé
library(dplyr)

# 1. Créer un data frame à partir de votre objet 'groups'
# L'objet 'groups' (créé avec cutree) contient déjà les noms des relevés et leur groupe.
df_groupes <- data.frame(Releve = names(groups), 
                         Groupe = groups)

# 2. Utiliser votre code dplyr (qui est parfait) pour résumer l'information
df_summary <- df_groupes %>%
  group_by(Groupe) %>%
  summarise(Releves_inclus = paste(Releve, collapse = ", "), .groups = 'drop')

# 3. Afficher le tableau final
print("Liste des relevés pour chaque groupe de la CAH :")
print(df_summary)

# Sauvegarder les résultats dans un même excel
library(openxlsx)

# 1. Créer un classeur Excel vide
wb <- createWorkbook()

# 2. Ajouter la première feuille et y écrire le premier tableau
addWorksheet(wb, "Especes_Indicatrices")
writeData(wb, "Especes_Indicatrices", indval_df_filtre)

# 3. Ajouter la deuxième feuille et y écrire le deuxième tableau
addWorksheet(wb, "Releves_Par_Groupe")
writeData(wb, "Releves_Par_Groupe", df_summary)

# 4. Enregistrer le fichier Excel sur votre ordinateur
# Le fichier s'appellera "Resultats_Analyse_Groupes.xlsx"
saveWorkbook(wb, file = "Resultats_Analyse_Groupes.xlsx", overwrite = TRUE)

