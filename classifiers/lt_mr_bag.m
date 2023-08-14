function [trainedClassifier, validationAccuracy] = lt_mr_bag(trainingData)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% Returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns as those imported into the app.
%
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double containing the accuracy as a
%       percentage. In the app, the Models pane displays this overall
%       accuracy score for each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict

% Auto-generated by MATLAB on 11-Nov-2022 11:36:32


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = {'bp_1_bipolar', 'bp_2_bipolar', 'bp_3_bipolar', 'bp_4_bipolar', 'bp_5_bipolar', 'fc_1_bipolar', 'coh_1_bipolar', 'coh_2_bipolar', 'coh_3_bipolar', 'coh_4_bipolar', 'coh_5_bipolar', 'coh_6_bipolar', 'spikes_1_car', 'rl_1_car', 'bp_1_car', 'bp_2_car', 'bp_3_car', 'bp_4_car', 'bp_5_car', 'fc_1_car', 'coh_1_car', 'coh_2_car', 'coh_3_car', 'coh_4_car', 'coh_5_car', 'coh_6_car'};
predictors = inputTable(:, predictorNames);
response = inputTable.soz_lats;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateTree(...
    'MaxNumSplits', 35, ...
    'NumVariablesToSample', 'all');
classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'Bag', ...
    'NumLearningCycles', 30, ...
    'Learners', template, ...
    'ClassNames', {'bilateral'; 'left'; 'right'});

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'bp_1_bipolar', 'bp_1_car', 'bp_2_bipolar', 'bp_2_car', 'bp_3_bipolar', 'bp_3_car', 'bp_4_bipolar', 'bp_4_car', 'bp_5_bipolar', 'bp_5_car', 'coh_1_bipolar', 'coh_1_car', 'coh_2_bipolar', 'coh_2_car', 'coh_3_bipolar', 'coh_3_car', 'coh_4_bipolar', 'coh_4_car', 'coh_5_bipolar', 'coh_5_car', 'coh_6_bipolar', 'coh_6_car', 'fc_1_bipolar', 'fc_1_car', 'rl_1_car', 'spikes_1_car'};
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2022a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = {'bp_1_bipolar', 'bp_2_bipolar', 'bp_3_bipolar', 'bp_4_bipolar', 'bp_5_bipolar', 'fc_1_bipolar', 'coh_1_bipolar', 'coh_2_bipolar', 'coh_3_bipolar', 'coh_4_bipolar', 'coh_5_bipolar', 'coh_6_bipolar', 'spikes_1_car', 'rl_1_car', 'bp_1_car', 'bp_2_car', 'bp_3_car', 'bp_4_car', 'bp_5_car', 'fc_1_car', 'coh_1_car', 'coh_2_car', 'coh_3_car', 'coh_4_car', 'coh_5_car', 'coh_6_car'};
predictors = inputTable(:, predictorNames);
response = inputTable.soz_lats;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
