%{
This is the primary script to run all the analyses for the epilepsy
laterality project.
%}

% If this is set to 1, it will re-run the full analysis, which takes
% several hours. If this is set to 0, it will look for .mat files
% containing intermediate results from the machine learning algorithm in
% your results path and load these in order to generate plots (this will
% take a few minutes at most).
do_full_pipeline = 1; 

% End users should not change this
doing_from_github = 1;

%% Get file locs
locations = epilepsy_laterality_locs;
plot_folder = locations.el_plots_folder;
script_folder = locations.el_script_folder;

% add the codebase to the path
addpath(genpath(script_folder))

% Remove the results.html file
if exist(fullfile(plot_folder,'results.html'),'file') ~= 0
    delete(fullfile(plot_folder,'results.html'))
end

if exist(fullfile(plot_folder,'supplemental_results.html'),'file') ~= 0
    delete(fullfile(plot_folder,'supplemental_results.html'))
end

%% Do model analyses
% This takes the intermediate file containing electrode-level features,
% calculates AI, and then does models to classify SOZ. This takes a while,
% in large part because it needs to run the subsampling analyses (takes
% 6 hours or so). It generates .mat files in the plot folder containing the
% model results. If this is set to 0, then the remaining parts of the
% function use the existing .mat model result files to generate figures. If
% set to 1, it replaces these .mat files with new ones.
if do_full_pipeline
    fprintf('\nDoing full pipeline, which will take at least several hours\n')
    model_ext_validation_files
end

%% Now do figures

% Figure 1 is methods figure

% Table 1
if ~doing_from_github
    % This requires loading files that have a lot of clinical data and are
    % too large to share
    make_table_1
end

% Table S1 is feature description table

% Figure S1
correlation_figure

% Figure 2 and Fig S2
combined_univariate_fmri_plots



% Figure 3 and Fig S3, S4, S5
model_plots

% Figure 4 and Fig S6 and S7
outcome_plots

close all

