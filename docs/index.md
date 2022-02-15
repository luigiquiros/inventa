---
title:  "Home"
layout: default
#permalink: /docs/index/
---

# INVENTA: Prioritization of natural extracts for chemical originality discovery

![Graphical abstract](/docs/assets/img/graphical_abstract.png)

## Description 

This is a workflow to explore the potential of a set of samples to contain not previously reported compounds. It is based on the results from data treatment usign MZmine, spectral organization through Molecular Networking, in-silico dereplication and prediction. 
It is composed of 4 independend components: 

The **Feature Component (FC)** is a ratio that considers compounds with a specificity higher than, for instance, 90% per sample (specified by the user), and without putative annotation. Results: one column with the FC ratio and one column with the sample specificity (ratio of peaks higher thant the specificied % with or without annotations).

The **Literature Component (LC)** is a score based on the number of compounds reported in the literature for the taxon. The output includes one column with the LC score and at least two additional columns of metadata containing the number of reported compounds at the species and genus levels.

The **Similarity Component (SC)** is a score based on the spectral similarity of the sample within the set. Multiple outlier detection machine learning algorithms are implemented to spot the dissimilar samples. A weight of ‘1’ is given to the sample considered anomalies in at least one detection method.

The **Class Component (CC)** is a score based on the presence of possible new chemical classes in the taxon, not previously reported before. The CC will be considered an integer ‘1’ if there are new chemical classes at the species level, and an additional ‘1’ if those new chemical classes are not present in the genus either.

The combined score (adition of the four components) can be modulated acording to the user preference. The ouput consist of .tsv file with all the information generared along the final rank of the samples.

![Rank conception](/docs/assets/img/priority_rank.png)