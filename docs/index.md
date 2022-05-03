---
layout: default
---
>  Graphical abstract
![Graphical abstract](/assets/img/graphical_abstract.png)


## Description 

This is a workflow to explore the potential of a set of samples to contain not previously reported compounds. It is based on the results from the data treatment usign MZmine, spectral organization through Molecular Networking, in-silico dereplication and prediction. 

It is composed of 4 independend components: 

The **Feature Component (FC)** is a percentage indicating the proportion of specific non annotated features by sample. A user-defined `minimum specificity` value is used to designate the minimum percentage at which a feature is considered specific, by default the value is set to 0.9 (90%). Additionally, a column wiht the Feature specificity (not considering the annotations) is calculated. This allows a direct visualization of the proportion of specific features that were annotated. A ratio of non annotated specific features with a hihg quality Molecular formula predictions is calculated if the Zodiac MF prediction from Sirius is included in the input.

The **Literature Component (LC)** is a score based on the number of compounds reported in the literature for the taxon. The output includes one column with the LC score and at least two additional columns of metadata containing the number of reported compounds at the species and genus levels. As the LC values approaches to 1, the lower number of compounds are reported in the literature.

The **Similarity Component (SC)** is a score based on the spectral similarity of the sample within the set. Multiple outlier detection machine learning algorithms are implemented to spot the dissimilar samples. A weight of ‘1’ is given to the sample considered anomalies in at least one detection method.

The **Class Component (CC)** is a score based on the presence of possible known chemical classes not previously reported in the taxon (new for the taxon). The CC will be considered an integer ‘1’ if there are new chemical classes at the species level. Both, not previously reported chemical classes in the species and in the genus are displayed next to the CC.

The combined score (adition of the four components) can be modulated acording to the user preference. The ouput consist of .tsv file with all the information generared along the final rank of the samples.

![Rank conception](/assets/img/priority_rank.png)
> Figure 1. Conceptual overview of the Priority rank and its individual components. FC: Feature Component. LC: Literature Component. CC: Class Component. SC: Similarity Component. A modulating factor (wn) allows to give a relative weight to each component according to the user preferences.
> 


<!-- toc -->

## Table of contents

### 1.[**Installation**](installation.md) 

### 2.[**Getting started**](getting-started.md) 

### 3.[**Configurations and running**](configuration-options.md)

<!-- tocstop -->

<!-- toc -->

## Publication interactive figures

### Fig 1.[**Conceptual overview of Inventa**](/assets/img/Detailed_priority_rank.png) 

### Fig 2.A. [**Celastraceae Library Chemical Space TMAP**](/assets/img/Celastraceae_annotation_vs_lotusdnp_tmap.html)
### Fig 2.B. [**Celastraceae Library Chemical Coverage**](/assets/img/Chemical_class_Celastraceae.html)  

### Fig 3.[**PCoA & UMAP projections for the Celastraceae library**](/assets/img/PCoA_UMAP_2D.html)

<!-- tocstop -->

#### Credits
Images were created by Luis Quiros-Guerrero using [bioRender](https://biorender.com/) (© BioRender 2022)

#### Copyright and license

Code and documentation copyright 2011–2022 the authors. Code released under the [MIT License](https://github.com/luigiquiros/inventa/blob/main/LICENSE).
