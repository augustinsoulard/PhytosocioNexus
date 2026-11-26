1. **On récupère les espèces du syntaxon** : les espèces caractéristiques et
   constantes citées dans la fiche PhytosocioNexus
   Quand la fiche n'en cite pas elle-même, on reprend celles du syntaxon supérieur.
2. **On regarde où chaque espèce est présente** : pour chaque département de France
   métropolitaine (Corse comprise), d'après nos sources de chorologie (FloreNum, eFlore,
   Préservons la nature).
3. **On calcule un pourcentage par département** : part des espèces du syntaxon qui sont
   présentes dans ce département. Exemple : si le syntaxon compte 10 espèces et que 9 sont
   dans le département, le pourcentage est de 90 %.
4. **On classe le département selon le rang du syntaxon** (lu sur le suffixe du nom :
   `-etum` association, `-ion` alliance, `-etalia` ordre, `-etea` classe…) :
   - associations et sous-associations : **100 %** = fortement probable, **> 80 %** = moyennement
     probable, **> 60 %** = faiblement probable ;
   - alliances et sous-alliances : mêmes seuils **abaissés de 10 points** (90 / 70 / 50) ;
   - en dessous du seuil bas, le département n'apparaît pas (laissé en blanc).
5. **La carte** (une image PNG par syntaxon) colorie les départements selon ces trois classes.

### Légende des couleurs

| Couleur                | Signification                          |
| ---------------------- | -------------------------------------- |
| blanc                  | département non retenu (sous le seuil) |
| `#aade8e` (vert clair) | faiblement probable                    |
| `#349200` (vert)       | moyennement probable                   |
| `#275c0a` (vert foncé) | fortement probable                     |

La carte est une **potentialité** : elle dit où le syntaxon *pourrait* se trouver d'après la
répartition départementale de ses espèces, pas où il a été observé directement.