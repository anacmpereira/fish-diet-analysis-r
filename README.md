# fish-diet-analysis-r
Statistical analysis in R comparing the diet of native and invasive fish using ecological indices, NMDS, PERMANOVA and GLMM.
# Native vs Invasive Fish Diet Analysis

Statistical analysis in **R** comparing the feeding ecology of the native fish *Phalloceros harpagos* and the invasive guppy *Poecilia reticulata* through ecological indices, multivariate analyses and generalized linear mixed models.

## About the project

This project investigates dietary patterns and trophic overlap between a native and a non-native fish species from freshwater ecosystems.

The workflow combines classical ecological metrics with modern statistical approaches to evaluate niche width, dietary overlap, community composition and prey abundance patterns. The analyses were developed as part of my research in aquatic ecology and biodiversity.

## Objectives

* Compare the diet of native and invasive fish species.
* Estimate trophic niche width.
* Evaluate dietary overlap between species.
* Analyze differences in prey composition.
* Identify the food items that contribute most to dietary variation.
* Produce reproducible statistical analyses and publication-quality figures.

## Species

* *Phalloceros harpagos* (native species)
* *Poecilia reticulata* (invasive species)

## Analyses performed

### Ecological indices

* Levins' Niche Width Index
* Pianka's Niche Overlap Index
* Schoener's Similarity Index
* Frequency of Occurrence

### Multivariate analyses

* Hellinger transformation
* Bray-Curtis dissimilarity
* Non-metric Multidimensional Scaling (NMDS)
* PERMANOVA
* Convex hull visualization

### Statistical modeling

* Generalized Linear Mixed Models (GLMM)
* Model comparison using AICc
* Overdispersion diagnostics
* Residual homogeneity assessment
* Coefficient visualization

## Workflow

<AsyncImage query="R statistical workflow diagram data cleaning analysis visualization" aspectRatio="16:6" width="100%" maxHeight=220/>

1. Import and organize dietary data.
2. Calculate ecological indices.
3. Transform abundance data.
4. Perform NMDS and PERMANOVA.
5. Fit GLMM models.
6. Validate model assumptions.
7. Compare competing models.
8. Generate publication-ready figures.

## Main packages

| Package       | Purpose                    |
| ------------- | -------------------------- |
| `tidyverse`   | Data manipulation          |
| `dplyr`       | Data processing            |
| `vegan`       | Community ecology analyses |
| `ggplot2`     | Data visualization         |
| `lme4`        | GLMM                       |
| `DHARMa`      | Model diagnostics          |
| `performance` | Model evaluation           |
| `car`         | Statistical tests          |
| `AICcmodavg`  | Model selection            |

## Example outputs

The script produces figures such as:

* NMDS ordination with convex hulls.
* GLMM coefficient plots.
* Mean prey abundance charts.
* Shepard plots for NMDS validation.


## Skills demonstrated

* Ecological data analysis
* Statistical modeling
* Community ecology
* Reproducible research
* Data visualization
* Scientific programming in R
* Model validation and diagnostics

## Author

**Ana Carolina Martins Pereira**

Biologist | MSc in Applied Ecology | PhD Candidate in Ecology, Biodiversity and Conservation

📧 [anamartinsper.bio@gmail.com](mailto:anamartinsper.bio@gmail.com)

🔗 LinkedIn
