Readme
Comparative data divergence database
Compiled in 2010 by Josef C Uyeda; uyedaj@science.oregonstate.edu
Description: This file includes node-averaged divergence values in darwin numerators (standardized by dimensionality) for all comparative data found in Uyeda et al. 2011. 

Data Columns:
Node: The number of the node in the phylogeny

Node.label: The name of the node in the phylogeny

Class: Taxonomic class of all descendants of the node

Order: Taxonomic order of all descendants of the node

Family: Taxonomic family of all descendants of the node

Genus: Taxonomic genus of all descendants of the node

Interval.my: Total branch length separating taxa at the given node 

Log10Years: Log10 of the number of years separating taxa at the given node 

d: Average divergence between all pairwise comparisons between taxa spanning the given node in darwin numerators standardized by dimensionality (k). All measurements represent linear body size measurements.

Sd: Standard deviation of the pairwise comparisons spanning the node

Taxa: Number of descendant tips
Comparisons: Number of pairwise comparisons that span the node 

Within.order: Whether or not tips spanning the node are in the same taxonomic order or not. 0-all comparisons spanning the node are between orders. 1- all comparisons spanning the node are within the same order. 0.5- Some comparisons spanning the node are within the same order (the order is not monophyletic)

Within.family: Whether or not tips spanning the node are in the same taxonomic family or not. 0-all comparisons spanning the node are between family. 1- all comparisons spanning the node are within the same family. 0.5- Some comparisons spanning the node are within the same family (the family is not monophyletic)

Within.genus: Whether or not tips spanning the node are in the same taxonomic genus or not. 0-all comparisons spanning the node are between genus. 1- all comparisons spanning the node are within the same genus. 0.5- Some comparisons spanning the node are within the same genus (the genus is not monophyletic)

Tree.source: Citation of the phylogeny source

Trait.source: Citation of the trait data source

Comments: Comments about use of the tree and source of time-calibration


