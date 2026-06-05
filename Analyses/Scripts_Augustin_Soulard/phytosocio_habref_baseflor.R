if(!require("tidyverse")){install.packages("tidyverse")} ; library("tidyverse")
if(!require("caulisroot")){devtools::install("../caulisroot")} ; library("caulisroot")

source("function/postgres/postgres_manip.R")
source("function/taxabase.R")
source("function/habitats/habref_manip.R")

con = copo() # Connexion à la bdd


# Chargement des données
habref_esp = dbGetQuery(con, '
                        SELECT 
  t2."cd_hab",
  t2."FG_VALIDITE",
  t2."cd_typo",
  t2."lb_code",
  t2."lb_hab_fr",
  t2."lb_hab_fr_complet",
  t2."LB_HAB_EN",
  t2."LB_AUTEUR",
  t2."niveau",
  t2."lb_niveau",
  t2."cd_hab_sup",
  t2."PATH_CD_HAB",
  t2."france",
  t2."lb_description",
  t1."cd_corresp_tax",
  t1."cd_hab_entree",
  t1."cd_nom",
  t1."CD_TYPE_RELATION",
  t1."LB_CONDITION",
  t1."LB_REMARQUES",
  t1."NOM_CITE",
  t1."nom_cite_match",
  t1."VALIDITE",
  t1."DATE_CREA",
  t1."DATE_MODIF" FROM habref.habref_corresp_taxon_70 AS t1
LEFT JOIN habref.habref_70 AS t2
  ON t1."cd_hab_entree" = t2."cd_hab"
WHERE t2."cd_typo" IN (\'4\', \'7\', \'8\', \'22\', \'18\', \'28\', \'107\', \'100\') 
                        and t2."france" = \'true\'
                          AND t1."cd_nom" IS NOT NULL
                          AND t1."cd_nom" <> \'\'
                          AND t2."cd_hab" IS NOT NULL
                          AND t2."cd_hab" <> \'\';'
                        )

habref_esp = rarete_occurence(habref_esp)

# Charger les données de QBiome
charger_gpkg(layers = c("Flore","Releve_Phyto"))

#Ajouter les cd_nom
Flore = findtaxa("Nom",Flore)

#Début de la fonction à créer
#Jointure flore et habiatts
florehab = left_join(st_drop_geometry(Flore),habref_esp,by="cd_nom", 
                     relationship = "many-to-many") %>% filter(!is.na(Releve))




#### Par nombre d'espèce
habitats_releve_esp<- florehab %>%
  distinct(Releve, cd_hab, lb_code, lb_hab_fr, lb_nom) %>%
  group_by(Releve, cd_hab,lb_code, lb_hab_fr,lb_nom) %>%
  summarise(nb = n(), .groups = "drop") %>%
  arrange(Releve, desc(nb)) %>%
  group_by(Releve,lb_nom) %>%
  slice_max(nb, n = 10, with_ties = FALSE) %>%
  ungroup()%>%
  select(Releve,lb_code, lb_hab_fr,lb_nom, nb )%>% filter(!is.na(lb_hab_fr)) 

habitats_releve <- florehab %>%
  distinct(Releve, cd_hab, lb_code, lb_hab_fr, lb_nom) %>%
  group_by(Releve, cd_hab,lb_code,  lb_hab_fr) %>%
  summarise(nb = n(), .groups = "drop") %>%
  arrange(Releve, desc(nb)) %>%
  group_by(Releve) %>%
  slice_max(nb, n = 10, with_ties = FALSE) %>%
  ungroup()%>%
  select(Releve,lb_code,  lb_hab_fr, nb)%>% filter(!is.na(lb_hab_fr))

releves <- unique(habitats_releve$Releve)

for (r in releves) {
  p <- habitats_releve %>%
    filter(Releve == r) %>%
    mutate(
      nom_habitat = paste(lb_code, lb_hab_fr),  # Nouveau label
      nom_habitat = fct_reorder(nom_habitat, nb)
    ) %>%
    ggplot(aes(x = nb, y = nom_habitat)) +
    geom_col(fill = "darkgreen") +
    labs(
      title = paste("Fréquence des habitats – Relevé", r),
      x = "Nb espèces caractéristiques",
      y = "Habitat"
    ) +
    theme_minimal(base_size = 11)
  
  print(p)
  ggsave(paste0("plot/", r, "_habitats.png"), plot = p, width = 9, height = 6)
}


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)

plot_esp_hab_top <- function(releve_id, df, top_n = 15, output_dir = "plot") {
  # Créer le dossier s'il n'existe pas
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Étape 1 : créer nom d’habitat
  df <- df %>%
    filter(Releve == releve_id, !is.na(lb_hab_fr)) %>%
    mutate(nom_habitat = paste(lb_code, lb_hab_fr))
  
  # Étape 2 : habitats les plus riches
  top_habitats <- df %>%
    distinct(nom_habitat, lb_nom) %>%
    count(nom_habitat, name = "nb_especes") %>%
    slice_max(nb_especes, n = top_n, with_ties = FALSE)
  
  # Étape 3 : données pour le graphique
  df_plot <- df %>%
    filter(nom_habitat %in% top_habitats$nom_habitat) %>%
    distinct(nom_habitat, lb_nom) %>%
    mutate(nom_habitat = fct_reorder(nom_habitat, nom_habitat, function(x) sum(x %in% top_habitats$nom_habitat)))
  
  # Étape 4 : graphique
  p <- ggplot(df_plot, aes(x = lb_nom, y = nom_habitat, fill = 1)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(guide = "none") +
    labs(
      title = paste("Espèces × Habitats les plus riches – Relevé", releve_id),
      x = "Espèces",
      y = paste0("Top ", top_n, " habitats")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
  
  # Affichage
  print(p)
  
  # Export PNG
  ggsave(
    filename = file.path(output_dir, paste0(releve_id, "_top", top_n, "_hab_esp.png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# 🔁 Application automatique à tous les relevés
releves <- unique(florehab$Releve)

for (r in releves) {
  plot_esp_hab_top(r, florehab)
}

# Par score #####

habitats_releve_esp <- florehab %>%
  filter(!is.na(lb_hab_fr)) %>%
  distinct(Releve, cd_hab, lb_code, lb_hab_fr, lb_nom, score_total) %>%
  group_by(Releve, cd_hab, lb_code, lb_hab_fr, lb_nom) %>%
  summarise(score_total_sum = sum(score_total, na.rm = TRUE), .groups = "drop") %>%
  arrange(Releve, desc(score_total_sum)) %>%
  group_by(Releve, lb_nom) %>%
  slice_max(score_total_sum, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  select(Releve, lb_code, lb_hab_fr, lb_nom, score_total_sum) %>%
  rename(nb = score_total_sum)  # pour conserver le nom attendu dans ton graphique


habitats_releve <- florehab %>%
  filter(!is.na(lb_hab_fr), !is.na(score_total)) %>%
  group_by(Releve, cd_hab, lb_code, lb_hab_fr) %>%
  summarise(score_total_sum = sum(score_total, na.rm = TRUE), .groups = "drop") %>%
  arrange(Releve, desc(score_total_sum)) %>%
  group_by(Releve) %>%
  slice_max(score_total_sum, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  select(Releve, lb_code, lb_hab_fr, score_total_sum) %>%
  rename(nb = score_total_sum)



releves <- unique(habitats_releve$Releve)

for (r in releves) {
  p <- habitats_releve %>%
    filter(Releve == r) %>%
    mutate(
      nom_habitat = paste(lb_code, lb_hab_fr),  # Nouveau label
      nom_habitat = fct_reorder(nom_habitat, nb)
    ) %>%
    ggplot(aes(x = nb, y = nom_habitat)) +
    geom_col(fill = "darkgreen") +
    labs(
      title = paste("Fréquence des habitats – Relevé", r),
      x = "Score des espèces caractéristiques",
      y = "Habitat"
    ) +
    theme_minimal(base_size = 11)
  
  print(p)
  ggsave(paste0("plot/", r, "_habitats.png"), plot = p, width = 9, height = 6)
}

## Gros calcule pour ggplot espèces
plot_top15_heatmap_habitats <- function(habitats_releve_esp, releve_id, output_dir = "plot") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Étape 1 : filtrer les données pour ce relevé
  df <- habitats_releve_esp %>%
    filter(Releve == releve_id, !is.na(nb)) %>%
    mutate(nom_habitat = paste(lb_code, lb_hab_fr))
  
  # Étape 2 : calcul du score total par habitat
  habitat_score <- df %>%
    group_by(nom_habitat) %>%
    summarise(score_total_habitat = sum(nb, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score_total_habitat)) %>%
    slice_head(n = 15)  # garder les 15 meilleurs
  
  # Étape 3 : filtrer les données uniquement pour ces habitats
  df_top <- df %>%
    filter(nom_habitat %in% habitat_score$nom_habitat) %>%
    left_join(habitat_score, by = "nom_habitat") %>%
    mutate(nom_habitat = fct_reorder(nom_habitat, score_total_habitat))
  
  # Étape 4 : graphique
  p <- ggplot(df_top, aes(x = lb_nom, y = nom_habitat, fill = nb)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "C", name = "Score espèce") +
    labs(
      title = paste("Top 15 habitats (score espèce) – Relevé", releve_id),
      x = "Espèces caractéristiques",
      y = "Habitat (par score total)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
  
  print(p)
  ggsave(
    filename = file.path(output_dir, paste0(releve_id, "_top15_esp_x_hab_heatmap.png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

for (r in unique(habitats_releve_esp$Releve)) {
  plot_top15_heatmap_habitats(habitats_releve_esp, releve_id = r)
}


