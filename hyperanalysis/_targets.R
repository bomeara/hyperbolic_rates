library(targets)
library(tarchetypes)
library(crew)

ncores_to_allocate <- 7


tar_option_set(
        packages = c("ggplot2", "TreeSim", "geiger", "ape", "dplyr", "Rmpfr", "tidyr", "parallel", "pbapply", "nloptr", "dentist", "scales", "hyperr8")#, controller = crew_controller_local(workers = ncores_to_allocate)
)

source("functions.R")

list(

  tar_target(all_data, aggregate_all_data()),
  tar_target(nreps_set, 100), 
  tar_target(hyperr8_analysis, hyperr8::hyperr8_run(all_data, nreps=nreps_set)),
  tar_target(hyperr8_analysis_yule_funny, hyperr8::hyperr8_run(get_funny_yule(), nreps=nreps_set)),
  tar_target(hyperr8_analysis_yule_funny_and_regular, merge_yule_funny_and_regular(hyperr8_analysis, hyperr8_analysis_yule_funny)),
  tar_target(hyperr8_analysis_all, rbind(hyperr8_analysis, hyperr8_analysis_yule_funny)),
  tar_target(supplemental_table, summarize_for_supplemental_table(hyperr8_analysis_all, all_r2_funny_and_regular)),
  tar_target(supplemental_table_saved, write.csv(supplemental_table, file="supplemental_table.csv")),
  tar_target(all_r2_regular, compute_all_possible_r2(hyperr8_analysis)),
  tar_target(all_r2_regular_quick_summary, summarize_all_r2(all_r2_regular)),
  tar_target(all_r2_funny, compute_all_possible_r2(hyperr8_analysis_yule_funny)),
  tar_target(all_r2_funny_quick_summary, summarize_all_r2(all_r2_funny)),
  tar_target(all_r2_funny_and_regular, rbind(all_r2_regular_quick_summary, all_r2_funny_quick_summary)),
  tar_target(animated, create_gganimate_plot(hyperr8_analysis)),
  #tar_target(hyperr8_unique_compare, get_unique_compared_to_original(hyperr8_analysis)),
  #tar_target(hyperr8_unique_compare_summarized, compute_percentiles(hyperr8_unique_compare)),
  #tar_target(hyperr8_unique_compare_summarized_csv, write.csv(hyperr8_unique_compare_summarized, file= "hyperr8_unique_compare_summarized.csv")),
  #tar_target(hyperr8_datum_percentiles, compute_per_datum_percentiles_faster(hyperr8_analysis)),
  #tar_target(hyperr8_datum_summarized_percentiles, summarize_per_datum_percentiles(hyperr8_datum_percentiles)),
#   tar_target(r2_results, compute_coefficient_of_determination(hyperr8_analysis)),
#   tar_target(r2_results_pretty, prettily_summarize_coefficient_of_determination(r2_results)),
#   tar_target(r2_results_saved, write.csv(r2_results_pretty, file="r2_results.csv")),
#   tar_target(r2_results_random_params, r2_from_prediction_from_randomizations(hyperr8_analysis)),
#   tar_target(r2_results_random_params_summarized, summarize_r2_from_prediction_from_randomizations(r2_results_random_params)),
#   tar_target(r2_merged, merge_r2_tables(r2_results_random_params_summarized, r2_results_pretty)),
#   tar_target(r2_merged_saved, write.csv(r2_merged, file="r2_merged.csv")),
   tar_target(hyperr8_analysis_saved, save_file_in_chunks(hyperr8_analysis)),
  tar_target(hyperr8_analysis_yule_funny_saved, save_file_in_chunks(hyperr8_analysis_yule_funny)),
  tar_target(raw_info_for_dentist, hyperr8::optimization_over_all_data(all_data, nstep_dentist=10000)) # yes, rerunning, which is wasteful, but fast
#   tar_target(r2_results_funny, compute_coefficient_of_determination(hyperr8_analysis_yule_funny)),
#   tar_target(r2_results_pretty_funny, prettily_summarize_coefficient_of_determination(r2_results_funny)),
#   tar_target(r2_results_random_params_funny, r2_from_prediction_from_randomizations(hyperr8_analysis_yule_funny)),
#   tar_target(r2_results_random_params_summarized_funny, summarize_r2_from_prediction_from_randomizations(r2_results_random_params_funny)),
#   tar_target(r2_merged_funny, merge_r2_tables(r2_results_random_params_summarized_funny, r2_results_pretty_funny)),
#   tar_target(r2_merged_saved_funny, write.csv(r2_merged_funny, file="r2_merged_funny.csv"))
)