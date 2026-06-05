# Chargement des bibliothèques nécessaires
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(shiny)
library(tidyr)
library(dplyr)
library(DT)        # Pour les tableaux interactifs
library(ape)       # Pour la visualisation de l'arbre de classification
library(data.table) # Pour utiliser dcast()
library(vegan)     # Pour NMDS et analyse de similarité
library(indicspecies) # Pour les espèces indicatrices avec multipatt()

# Fonction pour convertir les codes Braun-Blanquet en valeurs numériques (ex: '+' = 1, '1' = 2, etc.)
convert_bb <- function(x) {
  bb_codes <- c("+" = 1, "1" = 2, "2" = 3, "3" = 4, "4" = 5, "5" = 6)
  return(as.numeric(bb_codes[as.character(x)]))
}

# Interface utilisateur principale (barre de navigation)
ui <- navbarPage("Application Phytosociologique",
                 
                 # Onglet principal pour l'analyse
                 tabPanel("Analyse",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("file", "Charger un fichier CSV", accept = ".csv"),
                              uiOutput("select_releves"),
                              numericInput("n_clusters", "Nombre de groupes à identifier:", value = 2, min = 2),
                              actionButton("run", "Lancer l'analyse")
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Données brutes (head)", tableOutput("data_head")),
                                tabPanel("Données pivotées", DTOutput("pivoted_data"),downloadButton("download_pivoted", "Télécharger les données pivotées")),
                                tabPanel("Classification", plotOutput("clustering_plot")),
                                tabPanel("Ordination (NMDS)", plotOutput("nmds_plot")),
                                tabPanel("Espèces caractéristiques",
                                         HTML("<p><strong>Définition :</strong> Cette section utilise la fonction <code>multipatt()</code> du package <code>indicspecies</code> pour identifier les espèces les plus représentatives (indicatrices) des groupes de relevés définis par la classification hiérarchique. Elle calcule pour chaque espèce un score combinant sa fidélité (présence fréquente dans un groupe) et sa spécificité (présence exclusive dans ce groupe). Même sans p-value significative, une forte valeur d'indice peut indiquer une affinité marquée avec un groupe.</p>"),
                                         DTOutput("indval_table"),
                                         DTOutput("groupe_releves_table"),
                                         downloadButton("download_indval", "Télécharger les espèces caractéristiques"),
                                         downloadButton("download_groupes", "Télécharger les relevés par groupe")
                                )
                              )
                            )
                          )
                 ),
                 
                 # Onglet séparé pour visualiser les espèces par relevé
                 tabPanel("Liste des espèces par relevé",
                          DTOutput("especes_par_releve")
                 )
)

# Partie serveur de l'application
server <- function(input, output, session) {
  
  data_input <- reactive({
    req(input$file)
    read.csv(input$file$datapath, stringsAsFactors = FALSE, sep = ";")
  })
  
  output$data_head <- renderTable({
    head(data_input(), 10)
  })
  
  output$select_releves <- renderUI({
    req(data_input())
    releves <- unique(data_input()$releve)
    checkboxGroupInput("selected_releves", "Relevés à inclure:", choices = releves, selected = releves)
  })
  
  data_pivoted <- eventReactive(input$run, {
    df <- data_input()
    colnames(df) <- tolower(colnames(df))
    required_cols <- c("releve", "espece", "strate", "abondance_dominance")
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(paste("Colonnes manquantes dans le fichier :", paste(missing_cols, collapse = ", ")))
    }
    df <- df %>%
      filter(releve %in% input$selected_releves) %>%
      mutate(abondance_dominance = trimws(as.character(abondance_dominance))) %>%
      mutate(abondance = convert_bb(abondance_dominance)) %>%
      filter(!is.na(abondance)) %>%
      mutate(espece_strate = paste0(espece, "_", strate))
    df_dt <- as.data.table(df)
    df_pivot <- dcast(
      df_dt,
      releve ~ espece_strate,
      value.var = "abondance",
      fun.aggregate = sum,
      fill = 0
    )
    df_pivot <- as.data.frame(df_pivot)
    rownames(df_pivot) <- df_pivot$releve
    df_pivot <- df_pivot[ , !(names(df_pivot) %in% c("releve"))]
    df_pivot[] <- lapply(df_pivot, as.numeric)
    attr(df_pivot, "releve") <- df_dt[, unique(releve)]
    return(df_pivot)
  })
  
  output$pivoted_data <- renderDT({
    req(data_pivoted())
    datatable(data_pivoted(), options = list(pageLength = 10))
  })
  
  output$download_pivoted <- downloadHandler(
    filename = function() {
      paste0("donnees_pivotees_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- data_pivoted()
      releve <- attr(df, "releve")
      df_export <- cbind(releve = releve, df)
      write.csv2(df_export, file, row.names = FALSE)
    }
  )
  
  cluster_membership <- reactive({
    mat <- data_pivoted()
    dist_mat <- vegdist(mat, method = "bray")
    clust <- hclust(dist_mat, method = "ward.D2")
    cutree(clust, k = input$n_clusters)
  })
  
  output$clustering_plot <- renderPlot({
    mat <- data_pivoted()
    dist_mat <- vegdist(mat, method = "bray")
    clust <- hclust(dist_mat, method = "ward.D2")
    plot(clust, main = "Classification hiérarchique (Bray-Curtis)", xlab = "", sub = "")
    rect.hclust(clust, k = input$n_clusters, border = 2:6)
  })
  
  output$nmds_plot <- renderPlot({
    mat <- data_pivoted()
    groups <- factor(cluster_membership())
    nmds <- metaMDS(mat, k = 2, trymax = 100, autotransform = FALSE)
    nmds_sites <- as.data.frame(scores(nmds, display = "sites"))
    nmds_sites$label <- rownames(nmds_sites)
    nmds_sites$Groupe <- groups[rownames(nmds_sites)]
    
    stress_val <- round(nmds$stress, 3)
    
    ggplot(nmds_sites, aes(x = NMDS1, y = NMDS2)) +
      geom_point(aes(color = Groupe), size = 3) +
      ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 100) +
      labs(title = "Ordination NMDS",
           subtitle = paste("Stress:", stress_val),
           x = "NMDS 1", y = "NMDS 2") +
      theme_minimal() +
      coord_equal()
  })
  
  
  output$indval_table <- renderDT({
    mat <- data_pivoted()
    groups <- factor(cluster_membership())
    if (length(unique(groups)) < 2) {
      return(datatable(data.frame(Message = "Moins de 2 groupes détectés")))
    }
    indval_res <- multipatt(mat, groups, func = "IndVal.g", duleg = TRUE, control = how(nperm = 999))
    indval_df <- as.data.frame(indval_res$sign)
    if (nrow(indval_df) == 0) {
      return(datatable(data.frame(Message = "Aucune espèce caractéristique détectée (p > 0.05)")))
    }
    indval_df$Espèce <- rownames(indval_df)
    indval_df <- indval_df %>%
      filter(p.value <= 0.05) %>%
      mutate(Groupe = index) %>%
      select(Espèce, Groupe, stat, p.value)
    colnames(indval_df) <- c("Espèce", "Groupe indicateur", "Indice IndVal", "p-value")
    datatable(indval_df, options = list(pageLength = 10, dom = 'Blfrtip'), filter = 'top')
  })
  
  output$groupe_releves_table <- renderDT({
    clusters <- cluster_membership()
    df_groupes <- data.frame(Releve = names(clusters), Groupe = clusters)
    df_summary <- df_groupes %>% 
      group_by(Groupe) %>% 
      summarise(Relevés = paste(Releve, collapse = ", "))
    datatable(df_summary, options = list(pageLength = 5))
  })
  
  output$download_indval <- downloadHandler(
    filename = function() {
      paste0("especes_caracteristiques_", Sys.Date(), ".csv")
    },
    content = function(file) {
      mat <- data_pivoted()
      groups <- factor(cluster_membership())
      indval_res <- multipatt(mat, groups, func = "IndVal.g", duleg = TRUE, control = how(nperm = 999))
      indval_df <- as.data.frame(indval_res$sign)
      indval_df$Espèce <- rownames(indval_df)
      indval_df <- indval_df %>%
        filter(p.value <= 0.05) %>%
        mutate(Groupe = index) %>%
        select(Espèce, Groupe, stat, p.value)
      colnames(indval_df) <- c("Espèce", "Groupe indicateur", "Indice IndVal", "p-value")
      write.csv2(indval_df, file, row.names = FALSE)
    }
  )
  
  output$download_groupes <- downloadHandler(
    filename = function() {
      paste0("releves_par_groupe_", Sys.Date(), ".xlsx") # Extension .xlsx
    },
    content = function(file) {
      clusters <- cluster_membership()
      
      # 1. Création du tableau "Détail" (Non aggrégé : Relevé | Groupe)
      df_detail <- data.frame(Releve = names(clusters), Groupe = clusters)
      
      # 2. Création du tableau "Résumé" (Aggrégé)
      df_summary <- df_detail %>% 
        group_by(Groupe) %>% 
        summarise(Relevés = paste(Releve, collapse = ", "))
      
      # 3. Création du fichier Excel avec openxlsx
      wb <- createWorkbook()
      
      # Ajout de la feuille 1 : Résumé
      addWorksheet(wb, "Résumé par Groupe")
      writeData(wb, "Résumé par Groupe", df_summary)
      
      # Ajout de la feuille 2 : Détail (votre demande spécifique)
      addWorksheet(wb, "Détail Relevés")
      writeData(wb, "Détail Relevés", df_detail)
      
      # Sauvegarde
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$especes_par_releve <- renderDT({
    req(data_input())
    df <- data_input()
    colnames(df) <- tolower(colnames(df))
    required_cols <- c("releve", "espece", "strate", "abondance_dominance")
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      return(DT::datatable(data.frame(Erreur = paste("Colonnes manquantes:", paste(missing_cols, collapse=", ")))))
    }
    datatable(
      df %>% dplyr::select(releve, espece, strate, abondance_dominance),
      options = list(pageLength = 10, dom = 'Blfrtip'),
      filter = 'top'
    )
  })
}

shinyApp(ui = ui, server = server)
