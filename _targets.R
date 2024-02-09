library(targets)
library(tarchetypes)
library(crew)

ncores_to_allocate <- 7


tar_option_set(
        packages = c("ggplot2", "TreeSim", "geiger", "ape", "dplyr", "Rmpfr", "tidyr", "parallel", "pbapply", "nloptr", "dentist", "scales")#, controller = crew_controller_local(workers = ncores_to_allocate)
)

source("functions.R")

list(

  tar_target(all_data, aggregate_all_data()),
  tar_target(nreps_set, 100), 
  tar_target(hyperr8_analysis, hyperr8::hyperr8_run(all_data, nreps=nreps_set)),
  tar_target(hyperr8_analysis_yule_funny, hyperr8::hyperr8_run(get_funny_yule(), nreps=nreps_set)),
  tar_target(hyperr8_analysis_yule_funny_and_regular, merge_yule_funny_and_regular(hyperr8_analysis, hyperr8_analysis_yule_funny)),
  tar_target(hyperr8_unique_compare, get_unique_compared_to_original(hyperr8_analysis)),
  tar_target(hyperr8_unique_compare_summarized, compute_percentiles(hyperr8_unique_compare)),
  tar_target(hyperr8_unique_compare_summarized_csv, write.csv(hyperr8_unique_compare_summarized, file= "hyperr8_unique_compare_summarized.csv")),
  tar_target(hyperr8_datum_percentiles, compute_per_datum_percentiles_faster(hyperr8_analysis)),
  tar_target(hyperr8_datum_summarized_percentiles, summarize_per_datum_percentiles(hyperr8_datum_percentiles)),
  tar_target(r2_results, compute_coefficient_of_determination(hyperr8_analysis)),
  tar_target(r2_results_pretty, prettily_summarize_coefficient_of_determination(r2_results)),
  tar_target(r2_results_saved, write.csv(r2_results_pretty, file="r2_results.csv")),
  tar_target(r2_results_random_params, r2_from_prediction_from_randomizations(hyperr8_analysis)),
  tar_target(r2_results_random_params_summarized, summarize_r2_from_prediction_from_randomizations(r2_results_random_params)),
  tar_target(r2_merged, merge_r2_tables(r2_results_random_params_summarized, r2_results_pretty)),
  tar_target(r2_merged_saved, write.csv(r2_merged, file="r2_merged.csv")),
  tar_target(hyperr8_analysis_saved, write.csv(hyperr8_analysis, file=gzfile("hyperr8_analysis.csv.gz"))),
  tar_target(hyperr8_analysis_yule_funny_saved, write.csv(hyperr8_analysis_yule_funny, file=gzfile("hyperr8_analysis_yule_funny.csv.gz")))
)