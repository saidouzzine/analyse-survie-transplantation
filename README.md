# Analyse de survie après transplantation - Étude statistique sous R

## Présentation du projet

Ce projet propose une analyse statistique des trajectoires cliniques après transplantation chez des patients atteints de leucémie myéloïde chronique.

L'objectif est d'étudier la survie globale, la survenue d'une rechute et la mortalité après transplantation, en tenant compte de la censure, des événements concurrents et des transitions successives entre états cliniques.

L'étude repose sur une cohorte de **1 835 patients adultes**, sélectionnés à partir d'une base initiale de 1 977 patients. Six analyses complémentaires sont réalisées sous **R**, depuis l'estimation non paramétrique de la survie jusqu'aux modèles multi-états et à la régression de Fine-Gray

## Positionnement méthodologique

Ce projet s'inscrit dans une démarche d'application et d'approfondissement des méthodes d'analyse de survie enseignées en Master 2 de statistique et de biostatistique.

Il mobilise des outils classiques et avancés, notamment les estimateurs de Kaplan-Meier et d'Aalen-Johansen, le modèle de Cox, les méthodes de risques concurrents et la régression de Fine-Gray.

L'objectif n'est pas de développer de nouveaux estimateurs, mais de mettre en œuvre, comparer et interpréter rigoureusement plusieurs approches statistiques sur une même cohorte, en explicitant leurs hypothèses, leurs paramètres d'intérêt et leurs limites.

## Objectifs scientifiques

- Estimer la probabilité de survie globale après transplantation.
- Comparer la survie selon le score de risque de Gratwohl.
- Étudier l'association du score et de l'âge avec le risque instantané de décès.
- Estimer les probabilités de rechute et de décès sans rechute en présence de risques concurrents.
- Décrire l'évolution des probabilités d'être vivant sans rechute, vivant après rechute ou décédé.
- Étudier l'association du score avec l'incidence cumulée de rechute et son éventuelle variation au cours du temps.

## Données et population étudiée

La base initiale comprend **1 977 patients**. Après exclusion de 142 patients mineurs, la population d'analyse comporte **1 835 adultes**.

| Premier événement | Effectif |
|---|---:|
| Censure | 773 |
| Rechute | 421 |
| Décès sans rechute préalable | 641 |
| **Total** | **1 835** |

Au total, **809 décès** sont observés, dont **168 après une rechute**. Les durées de suivi sont mesurées depuis la transplantation et présentées en années (conversion des jours par division par 365,25).


## Méthodologie statistique

Le projet est organisé autour de six scripts R :

| Script | Analyse | Méthode |
|---|---|---|
| `01_import.R` | Préparation de la cohorte | Contrôles, sélection des adultes et définition des événements |
| `02_kaplan_meier.R` | Survie globale | Kaplan-Meier et test du log-rank |
| `03_cox.R` | Facteurs associés à la mortalité | Modèle de Cox et résidus de Schoenfeld |
| `04_competing_risks.R` | Rechute et décès sans rechute | Fonctions d'incidence cumulée |
| `05_multistate.R` | Trajectoires post-transplantation | Modèle *illness-death* et estimateur d'Aalen-Johansen |
| `06_fine_gray.R` | Incidence cumulée de rechute selon le score | Modèle de Fine-Gray, test de Gray et interactions temporelles |

Ces méthodes répondent à des questions statistiques distinctes et sont utilisées de manière complémentaire.

## Principaux résultats

### Survie globale - Kaplan-Meier

La survie globale estimée est de **63,9 % à un an**, **52,3 % à cinq ans** et **48,6 % à huit ans**.

À cinq ans, elle est de **75,2 %** pour le groupe *Low risk*, **50,5 %** pour *Medium risk* et **17,4 %** pour *High risk*. Le test global du log-rank met en évidence une différence entre les groupes (p = 2,41 × 10⁻²⁵).

### Modèle de Cox

Après ajustement sur l'âge et par rapport au groupe *Low risk*, les rapports de risques instantanés de décès sont estimés à :

- **HR = 2,19** pour *Medium risk* (IC à 95 % : [1,72 ; 2,80]).
- **HR = 4,47** pour *High risk* (IC à 95 % : [3,31 ; 6,04]).

Le diagnostic global fondé sur les résidus de Schoenfeld ne met pas en évidence de violation statistiquement détectable de l'hypothèse de proportionnalité (p = 0,20).

### Risques concurrents

À cinq ans, les incidences cumulées estimées sont de **25,80 % pour la rechute** et de **36,67 % pour le décès sans rechute**. Ces probabilités tiennent compte de la compétition entre les événements lorsqu'on étudie le premier événement après transplantation.

### Modèle multi-états

À cinq ans, les probabilités d'occupation estimées par Aalen–Johansen sont :

| État clinique | Probabilité à cinq ans |
|---|---:|
| Vivant sans rechute | 37,53 % |
| Vivant après rechute | 14,65 % |
| Décédé | 47,82 % |

Le modèle *illness–death* distingue les décès survenus sans rechute préalable de ceux observés après une rechute.

### Régression de Fine-Gray

Le modèle standard estime un rapport de risques de sous-distribution (**sHR**) de **1,32** pour *Medium risk* et de **1,83** pour *High risk*, par rapport à *Low risk*.

L'ajout d'interactions temporelles met en évidence une variation statistiquement significative des effets du score (test global de Wald : p = 0,000192). Les sHR constants doivent donc être interprétés avec prudence.

## Figures et visualisations

Le rapport scientifique comprend sept figures :

1. Courbe de survie globale estimée par Kaplan-Meier, avec intervalle de confiance à 95 %.
2. Courbes de survie globale selon le score de Gratwohl.
3. Diagnostic de proportionnalité des risques du modèle de Cox par les résidus de Schoenfeld.
4. Incidences cumulées de la rechute et du décès sans rechute.
5. Schéma du modèle *illness-death* à trois états et transitions autorisées.
6. Évolution des probabilités d'occupation des trois états estimées par Aalen-Johansen.
7. Incidences cumulées de rechute et de décès sans rechute selon le groupe de risque.

Les figures complètent les estimations numériques et facilitent la lecture des évolutions temporelles.

## Reproductibilité

### Environnement

Le projet utilise **R** et principalement les packages suivants :

- `survival` : estimation de Kaplan–Meier, modèle de Cox et diagnostics ;
- `cmprsk` : incidences cumulées, test de Gray et régression de Fine-Gray ;
- `mstate` : modèles multi-états et probabilités de transition.

Pour installer ces packages :

```r
install.packages(c("survival", "cmprsk", "mstate"))
```

### Exécution des analyses

1. Installer R et les packages requis (voir la section « Environnement »).
2. Les données sources sont fournies dans le dépôt : le fichier `data/leukemia.csv` est lu directement par le script d'importation. Aucun téléchargement n'est nécessaire.
3. Depuis la racine du projet, exécuter les scripts dans l'ordre :

```
Rscript scripts/01_import.R
Rscript scripts/02_kaplan_meier.R
Rscript scripts/03_cox.R
Rscript scripts/04_competing_risks.R
Rscript scripts/05_multistate.R
Rscript scripts/06_fine_gray.R
```
4. Comparer les sorties numériques, les diagnostics et les figures avec le rapport scientifique (`rapport_analyse_survie.pdf`).

**Remarque :** le script `01_import.R` prépare la cohorte (sélection des adultes, définition des événements) et génère le fichier intermédiaire `data/leucemie2.rds`, utilisé ensuite par `05_multistate.R`. Les scripts doivent donc être exécutés dans l'ordre indiqué.

## Limites de l'étude

Les analyses sont **observationnelles** et ne permettent pas d'établir des relations causales. Elles reposent notamment sur une hypothèse de censure suffisamment non informative et, selon les modèles, sur des hypothèses de proportionnalité ou de structure des transitions.

L'ajustement sur les facteurs pronostiques reste limité. Les estimations tardives peuvent également être moins précises lorsque le nombre de patients encore suivis diminue. Enfin, les diagnostics du modèle de Fine–Gray suggèrent une association entre le score et l'incidence cumulée de rechute qui varie au cours du temps.

## Rapport scientifique

Le rapport complet présente les fondements théoriques, la méthodologie, les résultats détaillés, les figures, les diagnostics statistiques, les limites et les références bibliographiques.

**Analyse de survie après transplantation — Rapport scientifique, octobre 2026.**

[Consulter le rapport complet (PDF)](rapport_analyse_survie.pdf)

## Auteur

**Saïd OUZZINE**  
Mathématiques appliquées, statistique et modélisation des données.

## Références bibliographiques

Les méthodes statistiques utilisées dans ce projet s'appuient sur un ouvrage en accès libre, des supports de cours universitaires français et la documentation officielle des packages R.

### A. Ouvrage numérique en accès libre

**[1]** Bender, A., & Sonabend, R. *Machine Learning for Survival Analysis*.  
[Consulter l'ouvrage en ligne](https://www.mlsabook.com/).

### B. Cours et polycopiés universitaires

**[2]** Saint-Pierre, P. *Introduction à l'analyse de survie (analyse des durées de vie)*. Support de cours universitaire, Université Toulouse III – Paul Sabatier, Toulouse, France.

**[3]** Saint-Pierre, P. *Analyse de survie, modèles multi-états et processus de comptage*. Support de cours universitaire, Université Toulouse III – Paul Sabatier, Toulouse, France.

**[4]** Geffray, S. (2008–2009). *Modèles multi-états, risques concurrents et événements récurrents*. Polycopié du cours de Modélisation en survie, Master 2 Pharmacologie clinique et épidémiologie, Université de Nantes, France.  
[Consulter le polycopié PDF](https://irma.math.unistra.fr/~geffray/cours/cours-nantes/M2coursCompetingRisksRecurrentEvents2008-2009.pdf).

### C. Logiciel et packages statistiques

**[5]** R Core Team. *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienne, Autriche.  
[Site officiel de R](https://www.R-project.org/).

**[6]** Therneau, T. M. *survival: Survival Analysis*. Package R, Comprehensive R Archive Network (CRAN).  
[Documentation officielle — survival](https://doi.org/10.32614/CRAN.package.survival).

**[7]** Gray, B. *cmprsk: Subdistribution Analysis of Competing Risks*. Package R, Comprehensive R Archive Network (CRAN).  
[Documentation officielle — cmprsk](https://doi.org/10.32614/CRAN.package.cmprsk).

**[8]** Putter, H., de Wreede, L. C., Fiocco, M., et al. *mstate: Data Preparation, Estimation and Prediction in Multi-State Models*. Package R, Comprehensive R Archive Network (CRAN).  
[Documentation officielle — mstate](https://doi.org/10.32614/CRAN.package.mstate).

---

Pour une présentation détaillée des fondements théoriques, de la méthodologie et des résultats statistiques, se référer au rapport scientifique associé à ce dépôt.
