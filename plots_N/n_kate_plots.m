function kate_plots

%{
The goal here is to make 2 separate plots: 1 showing the ability to
lateralize epilepsy using just fmri, and 1 showing the ability by combining
fmri and spikes
%}

%% Get file locs
locations = epilepsy_laterality_locs;
data_folder = locations.el_data_folder;
fmri_folder = [data_folder,'fmri_data/'];
plot_folder = locations.el_plots_folder;
script_folder = locations.el_script_folder;
addpath(genpath(script_folder))

%% Just fmri analysis
%% fmri locs
%file_path = [locations.main_folder,'Alfredo_code/fmri_analysis_AL_3_28_23/'];
file_path = fmri_folder;
csv_path = [file_path,'out_csvs/'];

%% Load files
T = readtable([file_path,'df.csv']);
bT = readtable([file_path,'BNA_subregions.xlsx']);

%% Remove controls
controls = strcmp(T.Final_Lat,'Control');
T(controls,:) = [];

%% Define temporal ROIs
temporal_hippo_amygdala_left = [108 110 112 114 116 118  74  78  86 212 214 216];
temporal_hippo_amygdala_right = [109 111 113 115 117 119  75  79  87 213 215 217];

% Add one to the indices of the regions (python to matlab)
temporal_hippo_amygdala_left = temporal_hippo_amygdala_left + 1;
temporal_hippo_amygdala_right = temporal_hippo_amygdala_right + 1;

% Show the names of the regions
lids = bT.LabelID_L;
rids = bT.LabelID_R;
regions = bT.Var6;
alt_regions = bT.LeftAndRightHemisphere;

left_regions = alt_regions(ismember(lids,temporal_hippo_amygdala_left));
right_regions = alt_regions(ismember(rids,temporal_hippo_amygdala_right));
assert(isequal(left_regions,right_regions))

%% Get fcon
npatients = size(T,1);
all_fcon = nan(npatients,246,246);

% Loop over the patients
for i = 1:npatients
    % Get the subject id
    subj_id = T.Subject{i};

    % Load the csv containing the patient-specific fcon
    fcon_T = readtable([csv_path,subj_id,'.csv']);
    fcon = table2array(fcon_T);
    all_fcon(i,:,:) = fcon;

end

%% Get the left strength and right strength
% Take abs value of connectivity and then mean across each temporal lobe
% (averaging all left-left connections and all right-right connections)
left_str = mean(abs(all_fcon(:,temporal_hippo_amygdala_left,temporal_hippo_amygdala_left)),[2 3]);
right_str = mean(abs(all_fcon(:,temporal_hippo_amygdala_right,temporal_hippo_amygdala_right)),[2 3]);

%% Define AI
AI = (left_str-right_str)./(left_str+right_str);
T.AI = AI;

%% Get lateralities
lat = T.Final_Lat;

%% Check lats
% Load manual validation file
mT = readtable('Manual validation.xlsx','Sheet','Comparing lats');


alats = cell(size(mT,1),1);
% loop over rows of this
for r = 1:size(mT,1)
    
    % grab rid
    rid = mT.RIDs(r);

    % turn it into format of T
    if rid < 100
        rid_text = sprintf('sub-RID00%d',rid);
    else
        rid_text = sprintf('sub-RID0%d',rid);
    end

    % find the matching row of Alfredo's table
    arow = strcmp(T.Subject,rid_text);

    if sum(arow) ~= 1, continue; end

    alats{r} = T.Final_Lat{arow};

end
alats(strcmp(alats,'R')) = {'right'};
alats(strcmp(alats,'L')) = {'left'};
alats(strcmp(alats,'B')) = {'bilateral'};

nmT = mT;
nmT = addvars(nmT,alats);

% remove empty
nmT(cellfun(@isempty,alats),:) = [];
assert(isequal(nmT.my_lats,nmT.alats))


% Now go the other way, what are the hup id numbers associated with
% Alfredo's RIDs?
hup_names = cell(size(T,1),1);
hup_lats = cell(size(T,1),1);
for r = 1:size(T,1)
    rid_text = T.Subject{r};
    % get just the number
    rid = str2num(strrep(rid_text,'sub-RID0',''));

    % find the matching RID in my lookup table
    mr = mT.RIDs == rid;

    if sum(mr) ~=1, continue; end
    hup_names{r} = mT.name{mr};
    hup_lats{r} = mT.my_lats{mr};
end

hup_lats(strcmp(hup_lats,'right')) = {'R'};
hup_lats(strcmp(hup_lats,'left')) = {'L'};
hup_lats(strcmp(hup_lats,'bilateral')) = {'B'};

naT = T;
naT = addvars(naT,hup_names);
naT = addvars(naT,hup_lats);

% remove empty
naT(cellfun(@isempty,hup_names),:) = [];
assert(isequal(naT.hup_lats,naT.Final_Lat))


%% change names
lat(strcmp(lat,'L')) = {'left'};
lat(strcmp(lat,'R')) = {'right'};
lat(strcmp(lat,'B')) = {'bilateral'};

%% Prediction analysis
% Do LOO CV
npts = length(AI);

% Prepare combined lats
comb_br = lat;
comb_br(contains(comb_br,'bilateral')|contains(comb_br,'right')) = {'br'};
pos_br = 'left';

comb_bl = lat;
comb_bl(contains(comb_bl,'bilateral')|contains(comb_bl,'left')) = {'bl'};
pos_bl = 'right';

% Initialize scores
all_scores = nan(2,npts);
alt_scores = nan(2,npts);

for ic = 1:2

    % Choose which lats to combine
    if ic == 1
        comb_lat = comb_br;
        pos = pos_br;
    elseif ic == 2
        comb_lat = comb_bl;
        pos = pos_bl;
    end

    comb_lat_bin = zeros(npts,1);
    comb_lat_bin(strcmp(comb_lat,pos)) = 1;

    fT = table(comb_lat,comb_lat_bin,AI,'VariableNames',{'lat','lat_bin','AI'});

    % Loop over patients for LOO
    for ip = 1:npts

        % divide training and testing data
        testing = ip;
        training = [1:ip-1,ip+1:npts];

        ft_test = fT(testing,:);
        ft_train = fT(training,:);

        % train the model on the training data
        mdl = fitglm(ft_train,'lat_bin ~ AI','Distribution','binomial'); % basic LR model
        coef = [mdl.Coefficients.Estimate(1);mdl.Coefficients.Estimate(2)];

        % test the model on the testing data
        score = glmval(coef,ft_test.AI,'logit');
        alt_scores(ic,ip) = 1/(1+exp(-coef(1)-coef(2)*ft_test.AI));
        all_scores(ic,ip) = score;
       
        
    end
end

% Get ROC curves
[XL,YL,~,AUCL] = perfcurve(comb_br,all_scores(1,:),'left');
[XR,YR,~,AUCR] = perfcurve(comb_bl,all_scores(2,:),'right');

% Plot ROC curves
figure
nexttile
ll = plot(XL,YL,'linewidth',2);
hold on
lr = plot(XR,YR,'linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
%
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',AUCR)},'fontsize',15,...
    'location','southeast')
%}
title({'fMRI alone (62 patients, LOO CV)'})
set(gca,'fontsize',15)


%% load mt data
inter_folder = locations.el_data_folder;
mt_data = load(fullfile(inter_folder,'mt_out_epilepsy_laterality.mat'));
mt_data = mt_data.out;

%% Do the lrmt
[eT,features] =  lr_mt(mt_data,3);

%% Add variable for fmri ai
temp_ai = nan(size(eT,1),1);
eT = addvars(eT,temp_ai,'NewVariableNames','fMRI_AI');

%% Find the overlapping patients

for i = 1:size(eT,1)
    % look at mT lookup to get RID
    row_mT = find(strcmp(eT.names{i},mT.name));

    if isempty(row_mT), continue; end
    rid = mT.RIDs(row_mT);

    % put in correct format
    rid_formatted = sprintf('sub-RID0%d',rid);

    % find the row of T
    row_T = find(strcmp(rid_formatted,T.Subject));
    if isempty(row_T), continue; end

    % add ai
    eT.fMRI_AI(i) = T.AI(row_T);

end

%% Find those that are not nan
aT = eT(~isnan(eT.fMRI_AI),:);

%% remove non temporal lobe SOZ
aT(~strcmp(aT.soz_locs,'temporal'),:) = [];

%% define features
spike_feature = 'spikes car sleep';
fmri_feature = 'fMRI_AI';
outcome = 'soz_lats';
soz_bin = nan(size(aT,1),1);
soz_bin(strcmp(aT.soz_lats,'left')) = 1;
soz_bin(strcmp(aT.soz_lats,'bilateral')) = 0;
nT = table(aT.names,aT.(spike_feature),aT.(fmri_feature),aT.(outcome),soz_bin,...
    'VariableNames',{'names','spikes','fmri','soz','soz_bin'});
assert(sum(isnan(nT.soz_bin))==0)

% note there are only left and bilateral SOZs, so just need a single model
%% LOO CV
npts = size(nT,1);
all_probs = nan(npts,1);
for i = 1:npts
    training = [1:i-1,i+1:npts];
    testing = i;

    nt_test = nT(testing,:);
    nt_train = nT(training,:);

    % train the model on the training data
    mdl = fitglm(nt_train,'soz_bin ~ fmri + spikes','Distribution','binomial'); % basic LR model
    coef = [mdl.Coefficients.Estimate(1);mdl.Coefficients.Estimate(2);mdl.Coefficients.Estimate(3)];

    % test the model on the testing data
    score = glmval(coef,[nt_test.spikes,nt_test.fmri],'logit');
    all_probs(i) = score;
end

% Get ROC curves
[XL,YL,~,AUCL] = perfcurve(nT.soz,all_probs,'left');

% Plot ROC curves
nexttile
ll = plot(XL,YL,'linewidth',2);
hold on
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
%
legend([ll],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    },'fontsize',15,...
    'location','southeast')
%}
title({'fMRI + iEEG spikes (6 patients, LOO CV)'})
set(gca,'fontsize',15)


end