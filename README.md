# hyperbolic_rates

Analyses used for hyperbolic rate paper

Overall data structure

* **BD_reanalyses_Yule_tests**: Contains the analyses for the simulated Yule trees and re-analysis of Henao-Diaz et al. data
* **Markov_tests**: Analyses for the effect of noise on rates on a Markov model
* **hyperanalysis**: The analyses making up the majority of the paper

Not included: the R package developed for analyses: https://github.com/bomeara/hyperr8. 

Within the **hyperanalysis** directory, the workflow is done with [`targets`](https://docs.ropensci.org/targets/) ([Landau 2021](https://doi.org/10.21105/joss.02959)).

* **_targets.R**: The main file for the workflow
* **functions.R**: Functions used in the workflow
* **figures.qmd**: Quarto document that renders to generate the figures
* **run.R**: Script to run the workflow: `Rscript run.R` to make it run.
* **data**: Directory containing the data used in the analyses. Note that we use data from other sources: definitely cite them, and it's best practice to go to them for the original data if possible (unless trying to do a reanalysis of our study, where any errors in our use of the data are part of what to assess).
* **outputs**: The analyses result in a huge data.frame (matrix), with one row per datapoint per model per replicate (1 original and 100 random replicates), 52,070,550 rows and 36 columns in total. This is too large as a single CSV file, so it has been broken into sections based on datasets and stored compressed in the outputs directory.
* **other files**: This includes the figures from the paper and supplement as well as full csv files of summarized results.
