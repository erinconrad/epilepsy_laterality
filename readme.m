%{
Epilepsy Laterality project

This is the codebase for the project using interictal EEG data to
lateralize temporal lobe epilepsy. Here we describe steps to begin with a
dataset containing electrode contact-level features and perform the
analysis to calculate patient-level asymmetry indices for each feature, and
then perform the machine learning algorithms to predict seizure onset zone
laterality and generate the figures from the epilepsy laterality paper.

To run the analysis, follow these steps:
1) Download the codebase from: 
    https://github.com/penn-cnt/cnt_tle_laterality
2) Download the datasets from: 
    https://upenn.box.com/s/67upxhl9wam135jn99jtjb72mmmz9aip
3) Create a file called "epilepsy_laterality_locs.m" that contains paths to
the codebase, the data folder, and the results folder. It should be
structured as follows:

    function locations = epilepsy_laterality_locs
    
    % Locations needed for epilepsy laterality project
    locations.el_data_folder = ***PATH TO THE DATA FOLDER***;
    locations.el_plots_folder = ***PATH TO THE RESULTS FOLDER (where plots will be output)***;
    locations.el_script_folder = ***PATH TO THE CODEBASE***;;
    
    end

Put this file in your path.

4) If you wish to only run the code to generate figures from intermediate
results files, then put the three .mat files from the Box results folder
into your results path (locations.el_plots_folder). If you wish to re-run 
the analysis from scratch, this is not necessary.
5) If you wish to only run the code to generate figures from intermediate
results files, then edit plots/do_all_results.m so that do_full_pipeline =
0.
6) Then, to run the analysis, navigate to plots/ and run:
    >> do_all_results


If you rerun only the code to generate figures, this took several minutes
(on a 2020 MacBook Air with an Apple M1 chip). 

Erin Conrad
University of Pennsylvania
August 2023
erin.conrad@uphs.upenn.edu

%}