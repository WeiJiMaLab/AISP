% subject_idx 1-9 belong to Experiment 1, and subject_idx 10-14 belong to
% Experiment 2

% Create log likelihood grid for subjects in Experiment 1
% 1 refers to the optimal model, 1:31 are the lambda (precision) value index.
% results are saved in the folder real_data_results/fix_s/raw_tables/opt
compute.generate_prediction_table_cluster(1:9, 1:31, 1)

% Combine segmented tables to a single table
% results are saved in the folder real_data_results/fix_x/opt/tables
compute.combine_tables(1:9, 1)

% Pick the fit parameters
% results are saved in the folder real_data_results/fix_x/opt/fit_pars
compute.fit_pars(1:9, 1)

% compute the evidence for each model for individual subjects, including
% bmc, aic, bic, aicc.
% results are saved in the folder real_data_results/fix_x/opt/evi
compute.evi(1:9, 1)

% compute the p(right|model,fit parameters) for each trial
% results are saved in the folder real_data_results/fix_x/opt/fit_results
compute.fit_predictions(1:9, 1)

% compute average p(right|model,fit parameters) for each quantile in the
% 2d space (target x distractor)
% results are saved in the folder
% real_data_results/fix_x/opt/fit_results_2d2
compute.fit_predictions_2d2(1:9, 1)

% compute average p(right|model,fit parameters) for specific target
% orientation ranges
% results are saved in the folder
% real_data_results/fix_x/opt/fit_results_1d
compute.fit_predictions_1d(1:9, 1)

% compute average p(right|model,fit parameters) for specific distractor
% orientation ranges
% results are saved in the folder
% real_data_results/fix_x/opt/fit_results_1d_dist
compute.fit_predictions_1d_dist(1:9, 1)

% plot psychometric curve (2d) for the real data
plots.psymetric_curve_2d2_data(1:9)

% plot psychometric curve (2d) for model prediction
plots.psymetric_curve_2d2_all(1:9, 1)

% plot psychometic curve (1d target) for real data
plots.psymetric_curve_1d_data(1:9)

% plot psychometic curve (1d target) for model prediction
plots.psymetric_curve_1d_all(1:9, 1)

% plot psychometic curve (1d distractor) for real data
plots.psymetric_curve_1d_dist_data(1:9)


% plot psychometic curve (1d distractor) for model_prediction
plots.psymetric_curve_1d_dist_all(1:9, 1)


