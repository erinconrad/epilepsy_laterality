function [T,features,Ts,sleep_info] =  lr_mt(mt_data,which_sleep_stages,rm_bad_spikes)

%{
This function takes a dataset (mt_data) containing electrode-contact-level
features, and calculates AI values for a bunch of features.
%}

%% Parameters
only_hup = 0; % do on all patients including MUSC

%% Default options for missing variables
if ~exist('which_sleep_stages','var')
    which_sleep_stages = [2 3]; % all = 1, wake =2, sleep = 3
end
if ~exist('rm_bad_spikes','var')
    rm_bad_spikes = 0;
end
which_montages = [1 2 3]; % machine = 1,car = 2, bipolar = 3

% Frequency band names
freq_names = {'delta','theta','alpha','beta','gamma','broadband'};

%% Load Manual validation file (in the github project)
T = readtable('Manual validation.xlsx','Sheet','outcome');
sT = readtable('Manual validation.xlsx','Sheet','EDF pipeline');
%sozT = readtable('Manual validation.xlsx','Sheet','SOZ');

%% get variables of interest
all_missing = cellfun(@isempty,mt_data.all_spikes(:,1,1));
names = mt_data.all_names;
npts = length(names);
surgery = mt_data.all_surgery;
resection_lat = mt_data.all_resec_lat;
resection_loc = mt_data.all_resec_loc;
ablation_lat = mt_data.all_ablate_lat;
ablation_loc = mt_data.all_ablate_loc;
engel_yr1 = mt_data.all_engel(:,1);
engel_yr2 = mt_data.all_engel(:,2);
ilae_yr1 = mt_data.all_ilae(:,1);
ilae_yr2 = mt_data.all_ilae(:,2);
soz_lats = mt_data.all_soz_lat;
soz_locs = mt_data.all_soz_loc;
disconnected = mt_data.all_disconnected;
all_n_wake_sleep_connected = mt_data.all_n_wake_sleep_connected;

% Find and exclude patients for whom bulk of record is disconnected
most_disconnected = sum(disconnected == 1,2) >= 0.9* size(disconnected,2);

% Identify patients for whom there is little wake or sleep connected
no_wake = all_n_wake_sleep_connected(:,1) < 5;
no_sleep = all_n_wake_sleep_connected(:,2) < 5;
n_wake = all_n_wake_sleep_connected(:,1);
n_sleep = all_n_wake_sleep_connected(:,2);
n_connected = sum(disconnected == 0,2);

% get number of symmetric labels
ref_labels = mt_data.all_labels(:,1);
n_symmetric = cellfun(@length,ref_labels);

% is all_missing ==1 the same as n_symmetric == 0?
assert(isequal(all_missing==1,n_symmetric==0)) % yes

%% Remove patients with bad spikes if I so choose
sT(cellfun(@isempty,sT.name),:) = [];
snames = sT.name;
assert(isequal(snames,names))
min_allowable = 25;
bi_good = sT.x_Correct_outOf50__bi_;
if rm_bad_spikes 
    bad_spikes = bi_good < min_allowable;
else
    bad_spikes = false(length(names),1);
end

%% Clean SOZ localizations and lateralities
soz_lats(cellfun(@isempty,soz_lats)) = {''};
soz_locs(cellfun(@isempty,soz_locs)) = {''};
soz_lats(strcmp(soz_lats,'diffuse')) = {'bilateral'}; % make diffuse laterality be the same as bilateral
soz_locs(contains(soz_locs,'temporal')) = {'temporal'}; % make any temporal be temporal

%% Get more specific soz localization

%% Consensus ablation or resection lat
surg_lat = cell(npts,1);
for i = 1:npts
    if isempty(resection_lat{i}) && isempty(ablation_lat{i})
    elseif strcmp(resection_lat{i},'na') && strcmp(ablation_lat{i},'na')
    elseif strcmp(resection_lat{i},'na')
        surg_lat{i} = ablation_lat{i};
    elseif strcmp(ablation_lat{i},'na')
        surg_lat{i} = resection_lat{i};
    else
        if ~strcmp(resection_lat{i},ablation_lat{i})
            error('what');
        end
        surg_lat{i} = resection_lat{i};
    end

end


%% Consensus ablation or reseciton loc
surg_loc = cell(npts,1);
for i = 1:npts
    if isempty(resection_loc{i}) && isempty(ablation_loc{i})
    elseif strcmp(resection_loc{i},'na') && strcmp(ablation_loc{i},'NA')
    elseif strcmp(resection_loc{i},'ATL')
        surg_loc{i} = 'temporal';
    elseif contains(ablation_loc{i},'temporal')
        surg_loc{i} = 'temporal';
    else
        surg_loc{i} = 'other';
    end
end

%% Fix the outcomes for the patients I manually confirmed
% I manually confirmed these myself, last updated May 2023
[engel_yr1,engel_yr2,ilae_yr1,ilae_yr2,surg_lat] = replace_with_my_outcomes(names,...
    engel_yr1,ilae_yr1,engel_yr2,ilae_yr2,T,soz_locs,surg_lat,surgery);


%% Parse surgery
surgery(cellfun(@isempty,surgery)) = {''};
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);
outcome(~resection_or_ablation) = {''}; % make non resection or ablation nan


%% Get features
% Initialize feature table
Ts = table(names,engel_yr1,engel_yr2,ilae_yr1,ilae_yr2,surgery,surg_lat,surg_loc,soz_locs,soz_lats,no_wake,no_sleep,n_wake,n_sleep,n_connected,n_symmetric,most_disconnected,bi_good);
features = {};

for which_sleep_stage = which_sleep_stages % all = 1, wake =2, sleep = 3;
    if which_sleep_stage == 1
        sleep_text = 'all';
    elseif which_sleep_stage == 2
        sleep_text = 'wake';
    elseif which_sleep_stage == 3
        sleep_text = 'sleep';
    end

    for which_montage = which_montages % machine = 1,car = 2, bipolar = 3
        
        if which_montage == 1
            montage_text = 'machine';
        elseif which_montage == 2
            montage_text = 'car';
        elseif which_montage == 3
            montage_text = 'bipolar';
        end
    
        
        coh = mt_data.all_coh(:,which_montage,which_sleep_stage);
        pearson = mt_data.all_pearson_squared(:,which_montage,which_sleep_stage);
        bp = mt_data.all_bp(:,which_montage,which_sleep_stage);
        rel_bp = mt_data.all_rel_bp(:,which_montage,which_sleep_stage);
        plv = mt_data.all_plv(:,which_montage,which_sleep_stage);
        rl = mt_data.all_rl(:,which_montage,which_sleep_stage);
        spikes = mt_data.all_spikes(:,which_montage,which_sleep_stage);
        xcor = mt_data.all_xcor(:,which_montage,which_sleep_stage);
        re = mt_data.all_re(:,which_montage,which_sleep_stage);
        se = mt_data.all_se(:,which_montage,which_sleep_stage);
        ll = mt_data.all_ll(:,which_montage,which_sleep_stage);

        coh_iqr = mt_data.all_coh_iqr(:,which_montage,which_sleep_stage);
        pearson_iqr = mt_data.all_pearson_squared_iqr(:,which_montage,which_sleep_stage);
        bp_iqr = mt_data.all_bp_iqr(:,which_montage,which_sleep_stage);
        rel_bp_iqr = mt_data.all_rel_bp_iqr(:,which_montage,which_sleep_stage);
        plv_iqr = mt_data.all_plv_iqr(:,which_montage,which_sleep_stage);
        rl_iqr = mt_data.all_rl_iqr(:,which_montage,which_sleep_stage);
        spikes_iqr = mt_data.all_spikes_iqr(:,which_montage,which_sleep_stage);
        xcor_iqr = mt_data.all_xcor_iqr(:,which_montage,which_sleep_stage);
        re_iqr= mt_data.all_re_iqr(:,which_montage,which_sleep_stage);
        se_iqr = mt_data.all_se_iqr(:,which_montage,which_sleep_stage);
        ll_iqr = mt_data.all_ll_iqr(:,which_montage,which_sleep_stage);

        % Find non-empty pts
        non_empty = find(cellfun(@(x) ~isempty(x), bp));
        first_non_empty = non_empty(1);
        
        % Loop over features
        for which_thing = {'spikes','rl','bp','se','pearson','xcor','coh','plv','ll','re'}
            switch which_thing{1}
                case {'pearson','inter_pearson','near_pearson'}
                    thing = pearson;
                    uni = 0; % not univariate
                    last_dim = 1; % one frequency
                case {'pearson_iqr'}
                    thing = pearson_iqr;
                    uni = 0; % not univariate
                    last_dim = 1; % one frequency
                case {'xcor','inter_xcor'}
                    thing = xcor;
                    uni = 0; % not univariate
                    last_dim = 1; % one frequency
                case {'xcor_iqr'}
                    thing = xcor_iqr;
                    uni = 0; % not univariate
                    last_dim = 1; % one frequency
                case {'coh','near_coh','inter_coh'}
                    thing = coh;
                    uni = 0; % not univariate
                    last_dim = size(coh{first_non_empty},3); % multiple frequencies
                case {'coh_iqr'}
                    thing = coh_iqr;
                    uni = 0; % not univariate
                    last_dim = size(coh{first_non_empty},3); % multiple frequencies
                case {'plv','near_plv','inter_plv'}
                    thing = plv;
                    uni = 0; % not univariate
                    last_dim = size(plv{first_non_empty},3); % multiple frequencies
                case {'plv_iqr'}
                    thing = plv_iqr;
                    uni = 0; % not univariate
                    last_dim = size(plv{first_non_empty},3); % multiple frequencies
                case {'re','inter_re'}
                    thing = re;
                    
                    
                    %thing = cellfun(@(x)
                    %exp(-x),thing,'UniformOutput',false); % do steve fix
                    for i = 1:length(thing)
                        x = thing{i};
                        x(isinf(x)) = nan; % set inf to undefined
                        thing{i} = x;
                    end

                    uni = 0; % not univariate
                    last_dim = size(re{first_non_empty},3); % multiple frequencies
                case {'re_iqr'}
                    thing = re_iqr;
                    uni = 0; % not univariate
                    last_dim = size(re{first_non_empty},3); % multiple frequencies
                case {'bp','inter_bp'}
                    thing = bp;
                    uni = 1; % univariate
                    last_dim = size(bp{first_non_empty},2); % multiple frequencies
                case {'bp_iqr'}
                    thing = bp_iqr;
                    uni = 1; % univariate
                    last_dim = size(bp{first_non_empty},2); % multiple frequencies
                case {'ll'}
                    thing = ll;
                    uni = 1; % univariate
                    last_dim = 1;%size(bp{1},2); % multiple frequencies
                case {'ll_iqr'}
                    thing = ll_iqr;
                    uni = 1; % univariate
                    last_dim = 1;%size(bp{1},2); % multiple frequencies
                case {'rel_bp','inter_rel_bp'}
                    thing = rel_bp;
                    uni = 1; % univariate
                    last_dim = size(rel_bp{first_non_empty},2);% multiple frequencies
                case {'rel_bp_iqr'}
                    thing = rel_bp_iqr;
                    uni = 1; % univariate
                    last_dim = size(rel_bp{first_non_empty},2); % multiple frequencies
                case {'se'}
                    thing = se;
                    uni = 1; % univariate
                    last_dim = 1;%size(bp{1},2); % multiple frequencies
                case {'se_iqr'}
                    thing = se_iqr;
                    uni = 1; % univariate
                    last_dim = 1;%size(bp{1},2); % multiple frequencies
                case 'spikes'
                    thing = spikes;
                    uni = 1; % univariate
                    last_dim = 1; % one frequency
                case 'spikes_iqr'
                    thing = spikes_iqr;
                    uni = 1; % univariate
                    last_dim = 1; % one frequency
                case 'nelecs'
                    thing = cellfun(@(x) ones(length(x),1),spikes,'uniformoutput',false);
                    last_dim = 1;
                    uni = nan;
                case {'rl','inter_rl'}
                    thing = rl;
                    uni = 1;
                    last_dim = 1;
                case {'rl_iqr'}
                    thing = rl_iqr;
                    uni = 1;
                    last_dim = 1;
            end
            
            labels = mt_data.all_labels(:,which_montage);
            
    
            %% Get asymmetry index
            ai_cell = cellfun(@(x,y) clean_ai_calc(x,y,uni,last_dim),labels,thing,...
                'uniformoutput',false);
            ai = cell2mat(ai_cell);
        
            
            %% AI table
            tnames_s = cell(last_dim,1);
            for i = 1:last_dim
                % PUT ACTUAL FREQ BAND IN TEXT
                clean_thing = which_thing{1};
                clean_thing = strrep(clean_thing,'_',' ');
                clean_thing = strrep(clean_thing,'iqr','SD');
                if last_dim == 1
                    tnames_s{i} = [clean_thing,' ',montage_text,' ',sleep_text];
                else
                    tnames_s{i} = [clean_thing,' ',freq_names{i},' ',montage_text,' ',sleep_text];
                end
            end
            features = [features;tnames_s];
        
        
            Ts = addvars(Ts,ai);
            Ts = splitvars(Ts,'ai','newVariableNames',tnames_s);
        end
    
    end

end



%% Pairwise correlations of all features
nfeatures = length(features); % -2 to remove outcome and bilaterality
all_feat = table2array(Ts(:,size(Ts,2)-nfeatures+1:end));

%% Prep  output table
% First, remove those rows missing all columns
T = Ts(~all_missing & ~most_disconnected & ~bad_spikes,:);
sleep_info = mt_data.all_n_wake_sleep_connected;
sleep_info = sleep_info(~all_missing & ~most_disconnected & ~bad_spikes,:);

%% Look for remaining nans
% The following steps demonstrate the reasons for nan values. I believe I
% have accounted for all of them.
if 0
    nT = T(T.no_wake==0,~contains(T.Properties.VariableNames,'rl')); % rl has nans a lot, also get nans if no wake
    nfeatures = features(~contains(features,'rl'));
    nA = table2array(nT(:,nfeatures));
    nan_rows = any(isnan(nA),2);
    nan_names = nT.names(nan_rows);
    
    % Show it
    nT(ismember(nT.names,nan_names),:)
    
    % remaining nan patients include:
    %{
    - HUP141: no bipolar montage
    - HUP160: spikes machine sleep - > no spikes in machine reference sleep
    - HUP172: spikes bipolar wake -> no spikes in bipolar reference wake
    - HUP202: spikes machine wake -> no spikes in machine reference wake 
    
    
    %}
end

%% Double check labels
if 0
good_labels = data.all_labels(not_missing,1);
for i = 1:length(good_labels)
    table(good_labels{i}(contains(good_labels{i},'A')|contains(good_labels{i},'B')|contains(good_labels{i},'C')))
    pause
end
end

%% Restrict to HUP patients
if only_hup
    hup_pt = contains(T.names,'HUP');
    T(~hup_pt,:) = [];
end
    

end