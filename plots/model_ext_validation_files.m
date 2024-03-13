function model_ext_validation_files

%{
This function re-runs the code to take intermediate dataset containing
electrode contact-level features, calculate patient-level AI values for
each feature, and runs the machine learning algorithms to predict epilepsy
laterality.
%}

%% Parameters (end-users should probably not change)
pca_perc = 95; % the percent variance to explain for pca
rm_non_temporal = 1; % remove patients who are not temporal
rm_wake = 0; % INCLUDE wake segments
which_refs = {'car','bipolar','machine'};

%% Get locations of various files and scripts
locations = epilepsy_laterality_locs;
plot_folder = locations.el_plots_folder;
inter_folder = locations.el_data_folder;
scripts_folder = locations.el_script_folder;
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
addpath(genpath(scripts_folder));

%% Load the file containing intermediate data
mt_data = load(fullfile(inter_folder,'mt_out_epilepsy_laterality.mat'));
mt_data = mt_data.out;

for ir = 1:length(which_refs)

    file_name = sprintf('ext_models_%s.mat',which_refs{ir});

    % Get the AI features
    if rm_wake == 1
        [T,features,~,sleep_info] =  lr_mt(mt_data,3); % the 3 refers to only looking at sleep
    else
        [T,features,~,sleep_info] =  lr_mt(mt_data,[2 3]); % wake and sleep
    end

    % Remove non-temporal
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
    sleep_info(~temporal,:) = [];

    %% Main analysis and good outcome analysis
    for ia = 1:2


        % Remove those without a response (soz_lats is the response variable)
        empty_class = cellfun(@isempty,T.soz_lats);
        T(empty_class,:) = [];
        sleep_info(empty_class,:) = [];
        

        % If doing good outcome analysis, restrict here
        if ia == 1 % all patients
            approach(ia).type = 'all patients';
            approach(ia).nums= [];
            
        else % good outcome only
            approach(ia).type = 'For HUP patients, only allow bilateral or unilateral good outcome (Engel 1)';
            which_year = 1;
            which_outcome = 'engel';

            surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
            outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
            outcome_bin = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
            good_outcome = strcmp(outcome_bin,'good') & surg == 1;

            good_outcome_unilat = good_outcome & (strcmp(T.soz_lats,'left') | strcmp(T.soz_lats,'right'));
            good_outcome_left = good_outcome & (strcmp(T.soz_lats,'left'));
            good_outcome_right = good_outcome & (strcmp(T.soz_lats,'right'));
            bilat = strcmp(T.soz_lats,'bilateral');

            % How many are there? I confirmed that I am removing 15 HUP
            % unilateral patients who either did not have surgery or did
            % not have good outcomes.
            
            hup = (contains(T.names,'HUP'));
            hup_unilat = hup & (strcmp(T.soz_lats,'left') | strcmp(T.soz_lats,'right'));
            hup_left = hup & strcmp(T.soz_lats,'left');
            hup_right = hup & strcmp(T.soz_lats,'right');
            hup_bilat = hup & bilat;
            hup_good_outcome_unilat = hup & good_outcome_unilat;
            hup_good_outcome_left = hup & good_outcome_left;
            hup_good_outcome_right = hup & good_outcome_right;

            if 0
            fprintf(['\n There are %d HUP patients (%d left, %d right, and %d bilat). \n'...
                'Of the unilateral patients, %d had good outcomes (%d left, %d right).\n'],...
                sum(hup),sum(hup_left),sum(hup_right),sum(hup_bilat),...
                sum(hup_good_outcome_unilat),sum(hup_good_outcome_left),sum(hup_good_outcome_right))
            end
            
            % allow if (MUSC) or (bilateral) or (unilateral and good outcome)
            allowed_outcome = contains(T.names,'MP')|good_outcome_unilat|bilat;
            T(~allowed_outcome,:) = [];
            sleep_info(~allowed_outcome,:) = [];
            approach(ia).nums.hup = sum(hup);
            approach(ia).nums.hup_left = sum(hup_left);
            approach(ia).nums.hup_right = sum(hup_right);
            approach(ia).nums.hup_bilat = sum(hup_bilat);
            approach(ia).nums.hup_good_outcome_unilat = sum(hup_good_outcome_unilat);
            approach(ia).nums.hup_good_outcome_left = sum(hup_good_outcome_left);
            approach(ia).nums.hup_good_outcome_right = sum(hup_good_outcome_right);
            approach(ia).sleep_info = sleep_info;
            
        end
        
        % Establish HUP and MUSC as training and testing, respectively
        train  = contains(T.names,'HUP');
        test  = contains(T.names,'MP');
        
       
        % Establish model types
        model(1).type = 'All features';
        model(2).type = 'Spikes';
        model(3).type = 'Binary spikes';
        nmodels = length(model);
        
        for im = 1:nmodels
            model(im).val(1).description = 'HUP cross-validation';
            model(im).val(2).description = 'MUSC external validation';
        end
        
        for im = 1:nmodels
            for iv = 1:2
                model(im).val(iv).side(1).description = 'Left vs right/bilateral';
                model(im).val(iv).side(2).description = 'Right vs left/bilateral';
            end
        end
        
        % Run all the models
        fprintf('\nDoing main models...');
        tic
        % Loop over model types
        for im = 1:nmodels
            just_spikes = im - 1; % if full model, just_spikes = 0, if spikes, just_spikes = 1, if binary, just_spikes = 2
            
            % If doing spikes only, restrict to sleep as well
            if just_spikes == 1 || just_spikes == 2
                curr_features = features(contains(features,'sleep'));
            else
                curr_features = features;
            end
           
            % Loop over internal vs external validation
            for iv = 1:2
        
                % Loop over sides
                for is = 1:2
        
                    % run the internal cross validation
                    if iv == 1
                        out = classifier_wrapper(T(train,:),curr_features ,pca_perc,is,...
                            just_spikes,rm_non_temporal,[],which_refs{ir});
            
                    % run the external test
                    elseif iv == 2
                        out = validation_classifier_wrapper(T,train,test,curr_features ,pca_perc,...
                            is,just_spikes,rm_non_temporal,which_refs{ir});
            
                    end
        
                    % Do the perfcurve
                    [out.X,out.Y,~,out.AUC] = perfcurve(out.class,out.scores,out.pos_class);
        
                    % Fill the appropriate structure entry
                    model(im).val(iv).side(is).result = out;
        
                end
                
            end
        end
        
        
        approach(ia).model = model;
        all.approach(ia) = approach(ia);
        all.which_ref = which_refs{ir};
        save([plot_folder,file_name],'all')
        
        
        % test - show results
        if 0
            figure
            tiledlayout(2,3)
            for iv = 1:2
                for im = 1:nmodels
                    nexttile
                    ll = plot(model(im).val(iv).side(1).result.X,...
                        model(im).val(iv).side(1).result.Y,'linewidth',2);
                    hold on
                    lr = plot(model(im).val(iv).side(2).result.X,...
                        model(im).val(iv).side(2).result.Y,'linewidth',2);
                    plot([0 1],[0 1],'k--','linewidth',2)
                    xlabel('False positive rate')
                    ylabel('True positive rate')
                    legend([ll,lr],{sprintf('%s: AUC = %1.2f',...
                        model(im).val(iv).side(1).description,model(im).val(iv).side(1).result.AUC),...
                        sprintf('%s: AUC = %1.2f',...
                        model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},'fontsize',15,...
                        'location','southeast')
                    title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
                    set(gca,'fontsize',15)
        
        
                end
            end
        end
    
    end
    
    %% Do subsampling analysis
    fprintf('done, took %1.1f seconds',toc);
    fprintf('\nStarting subsampling analysis (this will take several hours to a day).\n')
    tic
    
    % Run the lr_mt to extract features
    [T,features,way,dur,sample,ss,durations] =  lr_mt_multitime(mt_data,[2 3]); 
    empty_class = cellfun(@isempty,T.soz_lats);
    T(empty_class,:) = [];

    % Remove non-temporal
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];

    % Establish HUP and MUSC as training and testing, respectively
    train  = contains(T.names,'HUP');
    test  = contains(T.names,'MP');
    
    % Restrict to car spikes
    car_spikes = contains(features,sprintf('spikes_%s',which_refs{ir}));
    features = features(car_spikes);
    just_spikes = 1; % Just spikes
    
    % Establish what I am varying
    all_durs = unique(dur);
    ndurs = length(all_durs);
    
    all_samples = unique(sample);
    nsamples = length(all_samples);
    
    all_ss = unique(ss);
    nss = length(all_ss);
    
    % Initialize subsampling data
    cv_data = nan(nss,2,ndurs,nsamples);  % 2 sleep stages (wake and sleep); left vs right; which durs; which sample in dur
    ext_data = nan(nss,2,ndurs,nsamples); 
    
    % Loop over nss
    for iss = 1:nss
        
        curr_ss = all_ss(iss);
    
        cv_auc_l = nan(ndurs,nsamples);
        cv_auc_r = nan(ndurs,nsamples);
    
        ext_auc_l = nan(ndurs,nsamples);
        ext_auc_r = nan(ndurs,nsamples);
        
        % Loop over ndurs -  these will be different error bar points
        for id = 1:ndurs
            fprintf('\nDoing ss %d, dur %d...',iss,id);
            curr_dur = all_durs(id);
            
            
            % Loop over samples
            for is = 1:nsamples
                curr_sample = all_samples(is);
    
                % Get the relevant features
                relevant_features = contains(features,'_way1_') & ...
                    contains(features,sprintf('_ss%d_',curr_ss)) & ...
                    contains(features,sprintf('_dur%d_',curr_dur)) & ...
                    contains(features,sprintf('_samp%d_',curr_sample));
                curr_features = features(relevant_features);
    
                % Run the internal CV models
                left_int = classifier_wrapper(T(train,:),curr_features,...
                    pca_perc,1,just_spikes,rm_non_temporal,[],which_refs{ir});
                right_int = classifier_wrapper(T(train,:),curr_features,...
                    pca_perc,2,just_spikes,rm_non_temporal,[],which_refs{ir});
    
                % Get ROC stats
                [~,~,~,AUCL] = perfcurve(left_int.class,left_int.scores,left_int.pos_class);
                [~,~,~,AUCR] = perfcurve(right_int.class,right_int.scores,right_int.pos_class);
                cv_auc_l(id,is) = AUCL;
                cv_auc_r(id,is) = AUCR;
    
                % Run the external validation models
                left_ext = validation_classifier_wrapper(T,train,test,curr_features,...
                    pca_perc,1,just_spikes,rm_non_temporal,which_refs{ir});
                right_ext = validation_classifier_wrapper(T,train,test,curr_features,...
                    pca_perc,2,just_spikes,rm_non_temporal,which_refs{ir});
    
                % Get ROC stats
                [~,~,~,AUCL] = perfcurve(left_ext.class,left_ext.scores,left_ext.pos_class);
                [~,~,~,AUCR] = perfcurve(right_ext.class,right_ext.scores,right_ext.pos_class);
                ext_auc_l(id,is) = AUCL;
                ext_auc_r(id,is) = AUCR;
    
            end
    
        end
    
        % save data
        cv_data(iss,1,:,:) = cv_auc_l;
        cv_data(iss,2,:,:) = cv_auc_r;
    
        ext_data(iss,1,:,:) = ext_auc_l;
        ext_data(iss,2,:,:) = ext_auc_r;
    
    end
    fprintf('\nDone with subsampling analysis, took %1.1f seconds.\n',toc)
    
    all.cv_ss = cv_data;
    all.ext_ss = ext_data;
    all.durations = durations;
    
    save([plot_folder,file_name],'all')

end

end