Readme
Microevolutionary and paleontological divergence database
Compiled in 2010 by Josef C Uyeda; uyedaj@science.oregonstate.edu
Description: This file includes the Hendry et al. 2008 microevolutionary rate database, the Gingerich 2001 paleontological (mostly) rate database, and additional studies reported in the supplementary material of Uyeda et al. 2011. Only body size correlated traits were used in plots in the main text of Uyeda et al. 2011. 

Description of the data columns:

Data: A unique ID number given to each data point

Study: The study number, each organism within a single paper is given a study number

Trait: The trait number within the study

Sample 1: The name of the first population sample

Sample 2: The name of the second population sample

AuthorDate: The author and date of publication

Data.type: Either “Fieldstudy” if both populations are measured from extant populations, or “Fossil” if measurements are taken from fossil population samples.

Citation: Journal citation of the study

System: A brief description of the study system

Order: Taxonomic order of the species studied

Sex: 1- Male, 2- Females

GenPhen: Type of change observed; 1- Genetic, 2- Phenotypic

AlloSyn: Allochronic or synchronic divergence. Allochronic means that data is taken from the same population or lineage at two points in time, synchronic means that data is taken from two populations existing at the same time, and the interval is determined by doubling the time to the most recent common ancestor. Allo- 1, Syn- 2

Authonomous: Autonomous rates are taken between adjacent samples in a time-series. Non-autonomous divergence is measured between non-adjacent samples in time-series. Because non-autonomous rates are highly non-independent within a single time-series, we average non-autonomous 
divergence by binning all possible pairwise population divergence into equally spaced time bins. The number of bins equals the number of population samples in the time-series.

Trait: If available, the name of the trait measured

Dimension: The dimensionality of the measurement. 0- Unknown, 1- linear, 2- squared (e.g. area), 3- cubic (e.g. volume, mass)

BodySizeCorrelated: Whether or not the trait is a body size trait (mostly, any linear morphometric trait). 1- Body size trait, 2- Not correlated strongly to body size
Years: Number of years separating population samples. Synchronic intervals are twice what is reported in Hendry et al. 2008, because we take the total evolutionary time separating populations rather than 
the time to the ancestral population.

Glength: Number of years per generation

Darwins: Evolutionary rate in Darwins

Darwins.num: Darwin numerator

Darwins.num.abs: Absolute value of the Darwin numerator

Darwins (standardized by k): Darwins corrected by the dimensionality of the measurement

Haldanes: Evolutionary rate in Haldanes

Generations: Total number of generations separating population samples

log10.gen: Log10 of the total number of generations separating population samples

h.num: Haldane numeraotr

Darwin.Numerator: Absolute value of the Darwin numerator standardized by dimensionality, k

log10.years: Log10 of the total number of years separating populations

d: Darwin numerator divided by the dimensionality of the measurement, k. This measurement of divergence is used in the figures in the main text of Uyeda et al. 2011, plotted against log10.years.

InSitu: Whether the study system represents a disturbed system in which evolution is occuring in response to known biological disturbance, or represents in situ natural variation, as defined by Hendry et al. 2008. 1- in situ “natural” evolution 2- Disturbance-mediated evolution, NA- not enough information to assign to either class.
