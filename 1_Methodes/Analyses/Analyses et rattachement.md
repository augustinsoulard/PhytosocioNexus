---
tags:
date_maj: 2026-06-16
auteurs:
  - Augustin Soulard
sources:
---
Les associations végétales sont des groupes polythétiques. C'est-à-dire qu'elles ne sont pas strictement définies par un groupe d'espèce et deux individus d'associations ont souvent des compositions floristiques différentes. Il s'agit donc plus de groupes statistiques. Cependant, il existe de fortes propriétés émergentes de ces groupes qui sont statistiquement improbables d'un point de vu du nombre d'espèces et des combinaisons possibles. Pourtant, les mêmes associations végétales sont présentes dans différents endroits.

Pour les sceptiques de cette science qui pourrait être perçue comme abstraite. J'aime rappeler que les espèces sont aussi des groupes polythétiques dont la grande expérience humaine en botanique fait oublier les propriétés. Tout comme la phytosociologie, pour identifier une espèce, la botanique requiert normalement de regarder plusieurs individus de la populations et plusieurs critères sur chaque individus. Tout comme la phytosociologie nécessite de regarder plusieurs individus d'association et de regarder toutes les espèces (ou une grande partie) dans chaque association. En effet, il arrive qu'un critère pour l'identification d'une espèce ne colle pas sur certains individus de cette espèce et que certains individus de la population ne correspondent globalement pas aux critères d'identification de l'espèce. D'où la nécessité d'une approche statistique. Aujourd'hui cette technique n'est utilisée en botanique que pour la description des espèces ou pour certains groupes complexes. Pour les autres taxons, on se contente de clefs d'identification plus ou moins complexes. Pour la phytosociologie, nous ne connaissons pas encore assez les syntaxons pour avoir des clefs fiables excepté dans certaines régions très étudiées. Le rattachement se fait donc encore avec des méthodes statistiques plus complexes, mais souvent plus fiables que la lecture d'une clef.

Ainsi, nous disposons de **plusieurs outils pour effectuer des rattachements phytosociologiques** tels que :

- Les catalogues de végétation
- HABREF avec le PVF2
- Clefs d'identification des syntaxons
- Les méthodes statistiques multivariées avec des classifications (ex CAH) et des ordinations (ex AFC)
- [[index|PhytosocioNexus]]

D'anciennes méthodes de tri de tableaux (diagonalisation) existent, mais ces méthodes sont trop subjectives et doivent être évitées pour faire des rapprochements entre relevés.

## Les catalogues de végétation

Les catalogues permettent d'obtenir une liste des syntaxons sur un territoire donné. Cela permet d'aller regarder individuellement chaque syntaxon et de le comparer avec l'individu d'association à rattacher.

➕ **Avantages**

- Méthode simple et accessible

➖ **Inconvénients**

- Pas toujours disponible
- Pas forcément exhaustif et difficile de repérer une association non décrite

Il existe désormais [[Lafon et al. - 2024 - Catalogue de la végétation de France métropolitaine (CatVeg).pdf|CatVeg]], mais ce dernier ne liste que les syntaxons existants sans donner d'informations sur l'identification de ces derniers.

## HABREF avec le PVF2

HABREF est la base de données de référence pour les habitats naturels en France. On y retrouve les principales typologies dont plusieurs liées à la phytosociologie. La typologie de référence en phytosociologie aujourd'hui est celle du Prodrome des végétations de France déclinée (PVF2). Cette base de données offre différentes informations telles que les correspondances entre les typologies, les espèces associées ou même les statuts des habitats. C'est une référence à ne pas oublier, mais il manque encore beaucoup d'informations dedans.

➕ **Avantages**

- Open data
- Fonctionne dans toute la France
- Lié aux autres typologies d'habitats (EUNIS, Corine Biotope, Cahier d’habitats)

➖ **Inconvénients**

- Non exhaustif sur les espèces
- Nécessité de connaissances en manipulation de données pour un usage fluide

## Les clefs d'identification

Il s'agit de clefs dichotomiques classiques pour les syntaxons. Elles peuvent être faites directement sur le terrain ou à partir de relevés de végétations.

➕ **Avantages**

- Méthode simple et accessible

➖ **Inconvénients**

- Pas toujours disponible
- Pas forcément exhaustif et difficile de repérer une association non décrite

Depuis 2026, il existe la [[Botineau et Dupont - 2026 - Clé d'identification des unités phytosociologiques supérieures de France métropolitaine|Clé d'identification des unités phytosociologiques supérieures de France métropolitaine]]. Celle-ci est presque exhaustive et va jusqu'aux alliances et sous-alliances.

## Les méthodes statistiques

Il s'agit des méthodes de référence en matière de description de syntaxons. Ce sont également les méthodes les plus fiables et précises pour identifier des associations végétales.

Pour analyser des données en phytosociologie on utilise des analyses multivariées. Il s’agit d’un ensemble de méthodes statistiques permettant d'analyser simultanément plusieurs variables. En phytosociologie, elles sont particulièrement utiles pour étudier les relations complexes entre les espèces végétales et leur environnement. Par exemple, l'**analyse factorielle des correspondances (AFC)**, une technique d'analyse multivariée, est spécifiquement conçue pour étudier les tableaux de contingence. Elle permet de visualiser les relations entre deux variables qualitatives, comme la présence ou l'absence d'espèces dans différents relevés phytosociologiques. De même, la classification hiérarchique ascendante (ou CAH) est une méthode de **classification automatique** qui permet de regrouper des objets (ici, des relevés phytosociologiques) en classes de plus en plus larges. En **phytosociologie**, elle est utilisée pour identifier des **groupes de végétation** homogènes en termes de composition floristique. Le principe est de construire un arbre hiérarchique qui visualise les regroupements successifs des relevés.

➕ **Avantages**

- Permet de détecter les associations à décrire
- Méthode la plus précise

➖ **Inconvénients**

- Complexité de mise en œuvre
- Bases de données de référence manquantes ou non open data➡️ souvent nécessité de construire la sienne à partir d'articles en PDF ou de tableurs non standardisés.

### Rattachement via AFC DEDOU et CAH

- **Principe de fonctionnement**


L’objectif de ces analyses (pour un rattachement classique) est de comparer ses propres relevés avec des relevés de référence qui sont potentiellement similaires. Il y a donc une première étape de tri des relevés de terrains par grand type de végétation puis de collecte des relevés de référence.

Nous cherchons à mettre en lien les différents relevés selon les ‘variables’ espèces. Il faut donc un tableau de contingence qui croise les espèces par les relevés. On peut utiliser un simple tableau d’absence-présence, ou se servir des coefficients d’abondance-dominance.

- **Tableau de contingence**

|          | Relevé 1 | Relevé 2 | Relevé 3 |
| -------- | -------- | -------- | -------- |
| Espèce 1 | 0        | 1        | 1        |
| Espèce 2 | 0        | 1        | 0        |
| Espèce 3 | 1        | 1        | 1        |
L'**Analyse Factorielle des Correspondance (AFC)** permet de visualiser spatialement les similitudes entre les relevés. Elle possède autant d’axes qu’il y a d’espèce, mais seul les 2, 3 voire 4 axes les plus représentatifs de la variabilité des données sont interprétées. La méthode DEDOU consiste en le fait de créer des relevés ‘négatifs’. En reprenant notre exemple précédent on obtiendrait ceci :

**Tableau de la méthode DEDOU**
    
| Relevé 1 | Relevé 2 | Relevé 3 | Relevé 1’ | Relevé 2’ | Relevé 3’ |     |
| -------- | -------- | -------- | --------- | --------- | --------- | --- |
| Espèce 1 | 0        | 1        | 1         | 1         | 0         | 0   |
| Espèce 2 | 0        | 1        | 0         | 1         | 0         | 1   |
| Espèce 3 | 1        | 1        | 1         | 0         | 0         | 0   |
Cette méthode permet de lisser les effets des valeurs extrêmes pour avoir une AFC plus lisible à la projection graphique. A noter que les relevés ‘négatifs’ ne doivent pas être projetés graphiquement ni interprétés. Ils servent juste lors de la création de l’AFC.
    
La **Classification Ascendante Hiérarchique (CAH)** se base sur une matrice de distance réalisée à partir de tableau de contingence. Cela permet de réaliser des arbres classant les relevés de la même manière qu’un arbre phylogénétique classe des espèces. On a donc des relevés plus proches dans la hiérarchisation que d’autres et il peut s’agir de relevés réalisés dans des associations végétales similaires.

- **Matériel**

   - R et R studio
   - Relevé phytosociologique de référence et à identifier en format compatible R avec les relevés en colonne et les espèces en ligne. Il s'agit de la partie la plus délicate car aucune base de données propre et open data n'existe pour intégrer facilement les relevés de référence. Il faut donc la faire soit même à partir de références bibliographiques ou en demandant à d'autres phytosociologues.
   
**Comment choisir les relevés de références ?**

Il s'agit de tous les syntaxons qui sont potentiellement rattachables aux vôtres dans l'état actuel de vos connaissances. Dans l'idéal, il faudrait présenter à R plusieurs relevés pour chaque syntaxons potentiels. Vous trouverez quelques relevés de référence aux liens qui suivent :

- [Phytofrance | SFP (phytosocio.org)](http://phytosocio.org/#/phytofrance)
- [Tableaux synoptiques | SFP (phytosocio.org)](http://phytosocio.org/#/tableaux-synoptiques)

En l’absence de relevés en format de données accessible, les relevés de référence doivent être collectés directement dans les articles scientifiques et réécrit dans un format compatible avec R.

- **Exemple de script R**

```r
  # Chargement des packages    
 if(!require("vegan")){install.packages("vegan")} ; library("vegan")
   
   # Création du tableau DEDOU (relevés négatifs)
   DEDOU = tabContingence*-1+1
   
   # Création de l'Analyse factoriel des correspondances
   AFC<-cca(DEDOU)
   
   #Lecture des résultats
  summary(AFC)
   plot(AFC)
   
   # Extraire les scores des espèces et des sites
   species_scores <- as.data.frame(scores(AFC, display = "species"))
   sites_scores <- as.data.frame(scores(AFC, display = "sites"))
   
   #Création d'une matrice de distance
   mat_dist = as.matrix(dist(sites_scores))
   
   # Réalisation de la CHA
   arbre <- hclust(mat_dist, method = "ward.D2") # Méthode de Ward par défaut
   
  # Visualisation du dendrogramme
  plot(arbre)
  
```
  

Le script complet pour la phytosociologie est disponible [**ici**](https://github.com/augustinsoulard/Caulis/tree/main/script/phytosociologie) (sa réalisation est actuellement en cours et en teste donc c’est le fouillis je préviens !)

Pour plus de scripts R rendez-vous sur mon [**github**](https://github.com/augustinsoulard) !

>[!`caution`] Il est recommandé, quand cela est possible, de combiner plusieurs méthodes de rattachement.