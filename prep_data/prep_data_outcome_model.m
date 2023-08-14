function prep_data_outcome_model

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Initialize result variables

% Not reference-dependent
all_names = cell(npts,1);
all_stereo = nan(npts,1);
all_good_spikes = nan(npts,1);
all_outcome = cell(npts,2,2); % first is engel vs ilae, then 1 vs 2 years
all_surgery = cell(npts,1);
all_resec_lat = cell(npts,1);
all_resec_loc = cell(npts,1);
all_ablate_lat = cell(npts,1);
all_ablate_loc = cell(npts,1);
all_labels = cell(npts,2); % car vs bipolar
all_soz_locs = cell(npts,1);
all_soz_lats = cell(npts,1);


% Structure is car vs bipolar, then all vs wake vs sleep
all_coh = cell(npts,2,3);
all_bp = cell(npts,2,3);
all_spikes = cell(npts,2,3);
all_rl = cell(npts,2,3);
all_pearson = cell(npts,2,3);


%% Loop over patients
for p = 1:npts

    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    spikes = summ.spikes;
    rl = summ.rl;
    labels = summ.labels;
    soz_loc = summ.soz.loc;
    soz_lat = summ.soz.lat;
    name = summ.name;
    fc = summ.avg_fc;
    coh = summ.avg_coh;
    clinical = summ.clinical;
    all_stereo(p) = clinical.stereo;
    good_spikes = summ.good_spikes;
    all_good_spikes(p) = good_spikes;
    bp = summ.bp;
    bipolar_labels = summ.bipolar_labels;
    bipolar_fc = summ.avg_fc_bi;
    bipolar_coh = summ.avg_coh_bi;
    bipolar_bp = summ.bp_bi;
    spikes_ws = summ.spikes_ws;
    rl_ws = summ.rl_ws;
    fc_car_ws = summ.fc_car_ws;
    fc_bi_ws = summ.fc_bi_ws;
    coh_car_ws = summ.coh_car_ws;
    coh_bi_ws = summ.coh_bi_ws;
    bp_bi_ws = summ.bp_bi_ws;
    bp_car_ws = summ.bp_car_ws;


    %% Fill up basic stuff
    all_names{p} = name;
    all_soz_locs{p} = soz_loc;
    all_soz_lats{p} = soz_lat;
    
    %% Make spikes and rl nans if bad spikes
    if ~good_spikes
        spikes = nan(size(spikes));
        rl = nan(size(spikes));
    end    
 
    %% Unwrap coh
    bipolar_coh = wrap_or_unwrap_adjacency_fc_toolbox(bipolar_coh);
    coh = wrap_or_unwrap_adjacency_fc_toolbox(coh);
    coh_car_ws = cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,coh_car_ws,'UniformOutput',false);
    coh_bi_ws = cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,coh_bi_ws,'UniformOutput',false);
    
    %% Remove non-intracranial
    ekg = find_non_intracranial(labels);

    % Labels
    labels(ekg) = [];
    bipolar_labels(ekg) = [];

    % Spikes
    spikes(ekg,:) = [];
    spikes_ws{1}(ekg) = []; 
    spikes_ws{2}(ekg) = [];

    % RL
    rl(ekg,:) = [];
    rl_ws{1}(ekg) = []; 
    rl_ws{2}(ekg) = [];

    % FC
    fc(ekg,:) = []; fc(:,ekg) = [];
    bipolar_fc(ekg,:)=[]; bipolar_fc(:,ekg) = [];
    fc_car_ws{1}(ekg,:) = []; fc_car_ws{1}(:,ekg) = [];
    fc_car_ws{2}(ekg,:) = []; fc_car_ws{2}(:,ekg) = [];
    fc_bi_ws{1}(ekg,:) = []; fc_bi_ws{1}(:,ekg) = [];
    fc_bi_ws{2}(ekg,:) = []; fc_bi_ws{2}(:,ekg) = [];

    % Coherence
    coh(ekg,:,:) = [];
    coh(:,ekg,:) = [];
    bipolar_coh(ekg,:,:) = []; bipolar_coh(:,ekg,:) = [];
    coh_car_ws{1}(ekg,:,:) = []; coh_car_ws{1}(:,ekg,:) = [];
    coh_car_ws{2}(ekg,:,:) = []; coh_car_ws{2}(:,ekg,:) = [];
    coh_bi_ws{1}(ekg,:,:) = []; coh_bi_ws{1}(:,ekg,:) = [];
    coh_bi_ws{2}(ekg,:,:) = []; coh_bi_ws{2}(:,ekg,:) = [];

    % BP
    bp(ekg,:,:) = [];
    bipolar_bp(ekg,:,:) = [];
    bp_car_ws{1}(ekg,:) = []; 
    bp_car_ws{2}(ekg,:) = [];
    bp_bi_ws{1}(ekg,:) = []; 
    bp_bi_ws{2}(ekg,:) = [];
   
    %% Fill up Labels
    all_labels{p,1} = labels;
    all_labels{p,2} = bipolar_labels;
    
    %% Get averages for spikes and rl over times
    avg_spikes = nanmean(spikes,2);
    avg_rl = nanmean(rl,2);
    avg_bp = nanmean(bp,3);
    avg_bp_bipolar = nanmean(bipolar_bp,3);

    %% Fill up other stuff
    % Outcome and surgery
    all_outcome{p,1,1} = clinical.engel{1}; all_outcome{p,1,2} = clinical.engel{2}; % first is engel vs ilae, then 1 vs 2 years
    all_outcome{p,2,1} = clinical.ilae{1}; all_outcome{p,2,2} = clinical.ilae{2};

    all_surgery{p} = clinical.surgery;
    all_resec_lat{p} = clinical.resection_lat;
    all_resec_loc{p} = clinical.resection_loc;
    all_ablate_lat{p} = clinical.ablation_lat;
    all_ablate_loc{p} = clinical.ablation_loc;

    %% Other stuff. Structure is car vs bipolar, then all vs wake vs sleep

    % Coherence
    all_coh{p,1,1} = coh; all_coh{p,1,2} = coh_car_ws{1}; all_coh{p,1,3} = coh_car_ws{2};
    all_coh{p,2,1} = bipolar_coh; all_coh{p,2,2} = coh_bi_ws{1}; all_coh{p,2,3} = coh_bi_ws{2};

    % Bandpower
    all_bp{p,1,1} = avg_bp; all_bp{p,1,2} = bp_car_ws{1}; all_bp{p,1,3} = bp_car_ws{2};
    all_bp{p,2,1} = avg_bp_bipolar; all_bp{p,2,2} = bp_bi_ws{1}; all_bp{p,2,3} = bp_bi_ws{2};

    % Spikes
    all_spikes{p,1,1} = avg_spikes; all_spikes{p,1,2} = spikes_ws{1}; all_spikes{p,1,3} = spikes_ws{2};

    % RL
    all_rl{p,1,1} = avg_rl; all_rl{p,1,2} = rl_ws{1}; all_rl{p,1,3} = rl_ws{2};

    % Pearson correlation
    all_pearson{p,1,1} = fc; all_pearson{p,1,2} = fc_car_ws{1}; all_pearson{p,1,3} = fc_car_ws{2};
    all_pearson{p,2,1} = bipolar_fc; all_pearson{p,2,2} = fc_bi_ws{1}; all_pearson{p,2,3} = fc_bi_ws{2};

    
    %% Error checking
    if 0
        figure
        tiledlayout(1,2)
        nexttile; turn_nans_gray(all_pearson{p,1,1})
        nexttile; turn_nans_gray(all_pearson{p,1,3})
    end

    if 0
        figure
        plot(all_spikes{p,1,2},all_spikes{p,1,3},'o')
    end
    

end

out.all_soz_locs = all_soz_locs;
out.all_soz_lats = all_soz_lats;
out.all_spikes = all_spikes;
out.all_rl = all_rl;
out.all_pearson = all_pearson;
out.all_coh = all_coh;
out.all_labels = all_labels;
out.all_names = all_names;
out.all_stereo = all_stereo;
out.good_spikes = all_good_spikes;
out.outcome = all_outcome;
out.all_surgery = all_surgery;
out.all_resec_lat = all_resec_lat;
out.all_resec_loc = all_resec_loc;
out.all_ablate_lat = all_ablate_lat;
out.all_ablate_loc = all_ablate_loc;
out.all_bp = all_bp;


%% Save
save([out_folder,'main_out.mat'],'out');

end