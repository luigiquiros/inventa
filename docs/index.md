---
layout: default
---
>  Graphical abstract
![Graphical abstract](/assets/img/graphical_abstract.png)


## Description 

Here, we introduce Inventa, an untargeted mass spectrometry-based prioritization workflow for natural products extract libraries. The bioinformatic workflow is composed of four components that aim at estimating the potential of natural extracts for chemical novelty: 

"The feature component" considers the presence (intensity/area) and specificity of the mass spectrometric features along with their annotation status (annotated/ unannotated) with experimental and/or in silico databases. The postulate of this score is that the presence of a highly specific and non annotated metabolome is an indication of a particular chemistry. 

"The literature component" gives a score based on the number of compounds reported in the literature at the species, genus and family level. The more a plant is studied  (species, genus, family), the less likely it is to find structurally new compounds. 

"The class component" takes advantage of the prediction capacity of chemical classes based only on the fragmentation pattern without a formal putative structure annotation through [**CANOPUS**](). The predicted chemical classes in each sample are compared against those reported in the literature. The proposition of new chemical classes in particular samples being a potential indication  of new compounds. 

"The similarity component" harnesses the spectral diversity of the samples through vectorization of the fragmentation spectra data with [**MEMO**](https://doi.org/10.3389/fbinf.2022.842964) and application of machine learning outliers detectors. This component assumes that samples classified as ‘outlier’ hold a particular pool of metabolites  with a specific particular chemistry.

The combined score (adition of the four components) can be modulated acording to the user preference. The ouput consist of .tsv file with all the information generared along the final rank of the samples.

![Rank conception](/assets/img/priority_rank.png)
>Figure 1. Conceptual overview of Inventa’s priority rank and its individual components. (A) Feature Component (FC): a ratio from 0 to 1 considering the proportion of features with a ‘specificity’ higher than the minimum specificity value set by the user and lack of annotation per sample. The ‘Feature specificity’ shows the proportion of features with a ‘specificity’ higher than the minimum specificity value set by the user per sample. (B) Literature Component (LC): a ratio from 0 to 1 considering a ‘minimum of reported compounds’ set by the user and the number of compounds reported in the databases for the taxon. The closer the value to ‘1’ the less compounds are reported for the sample. Additional columns show the total number of reported compounds at the species and genus levels. (C) Similarity Component (SC): a score considering the spectral dissimilarity of the samples within the set. A machine learning algorithm is used to detect the outliers (value ‘1’) from the rest (value ‘0’) in the set based on a vectorized matrix. (D) Class Component (CC): a score considering the presence of predicted known chemical classes new to the species (value ‘1’). It is a product of the difference found between the chemical classes predicted by SIRIUS and the ones reported in the databases for a specific species. Additional columns show the specific ‘new’ known chemical classes at the species and genus levels. The Priority rank (PR) is the addition of the four components. A modulating factor (wn) allows to give a relative weight to each component according to the user’s preferences. The higher the value the higher is the rank of the sample. 
> 


<!-- toc -->

## Table of contents

### 1.[**Installation**](installation.md) 

### 2.[**Getting started**](getting-started.md) 

### 3.[**Configurations and running**](configuration-options.md)

<!-- tocstop -->

<!-- toc -->

##
##
## Publication interactive figures

### Fig 2.A. [**Celastraceae Library Chemical Space TMAP**](/assets/img/Celastraceae_annotation_vs_lotusdnp_tmap.html)
### Fig 2.B. [**Celastraceae Library Chemical Coverage**](/assets/img/Chemical_class_Celastraceae.html)  

### Fig 3.[**PCoA & UMAP projections for the Celastraceae library**](/assets/img/PCoA_UMAP_2D.html)

<!-- tocstop -->

#### Credits
Images were created by Luis Quiros-Guerrero using [bioRender](https://biorender.com/) (© BioRender 2022)

#### Copyright and license

Code and documentation copyright 2011–2022 the authors. Code released under the [MIT License](https://github.com/luigiquiros/inventa/blob/main/LICENSE).
