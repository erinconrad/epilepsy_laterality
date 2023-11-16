function outcome_plots

%% Parameters
which_year = 1;
which_model = 'spikes';
which_refs = {'car','bipolar','machine'};


%% Get file locs
locations = epilepsy_laterality_locs;
plot_folder = locations.el_plots_folder;
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end


%% Load the file containing intermediate data
inter_folder = locations.el_data_folder;
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% Load the MUSC outcome file
muscT = readtable([inter_folder,'LEN patient list research erin.xlsx']);


%% Initialize results file
for ir = 1:length(which_refs)

    if ir == 1
        fname = [plot_folder,'results.html'];
        fid = fopen(fname,'a');
        fprintf(fid,'<br><u><i>Concordance between spike-predicted laterality and surgical laterality is higher for patients with good surgical outcomes</i></u></br>');
        
      
    end
    
    %% Load the model file
    out = load([plot_folder,sprintf('ext_models_%s.mat',which_refs{ir})]);
    out = out.all;
    
    % which model
    switch which_model
        case 'full'
            model = out.approach(1).model(1).val(1);
            musc_model = out.approach(1).model(1).val(2);
        case 'spikes'
            model = out.approach(1).model(2).val(1);
            musc_model = out.approach(1).model(2).val(2);
    end
    
    %% Run the mt_lr again just to get overall outcome stuff
    T =  lr_mt(mt_data,3);
    empty_class = cellfun(@isempty,T.soz_lats);
    T(empty_class,:) = [];
    temporal_loc = contains(T.soz_locs,'temporal');
    T(~temporal_loc,:) = [];
    hup = contains(T.names,'HUP');
    musc = contains(T.names,'MP');

    %% Get MUSC outcomes
    % Loop over patients in T
    if 1
    for ip = 1:size(T,1)

        if ~contains(T.names,'MP')
            continue
        end

        % find matching musc table patient
        musc_row = contains(muscT.LENID_,T.names{ip});
        if sum(musc_row) == 1
            if contains(muscT.EngelYear1{musc_row},'N/a') % don't change format if empty outcome
                T.engel_yr1{ip} = '';
                T.engel_yr2{ip} = '';
                T.ilae_yr1{ip} = '';
                T.ilae_yr2{ip} = '';
                continue
            end

            % get outcomes
            T.engel_yr1{ip} = muscT.EngelYear1{musc_row};
            T.engel_yr2{ip} = muscT.EngelYear2{musc_row};
            T.ilae_yr1{ip} = muscT.ILAEYear1{musc_row};
            T.ilae_yr2{ip} = muscT.ILAEYear2{musc_row};

            % get surgery
            if strcmp(muscT.TypeOfSurgery{musc_row},'ATL')
                T.surgery{ip} = 'Resection';
                T.surg_lat{ip} = lower(muscT.SurgeryLaterality{musc_row});
                T.surg_loc{ip} = 'temporal';
            end

        end


    end
    end


    %% Remove patients (should be one patient for bipolar) with nan feature
    features = T.Properties.VariableNames;
    spike_features = features(contains(features,'spikes') & contains(features,which_refs{ir}));
    nan_feature = isnan(T{:,spike_features});
    T(nan_feature,:) = [];    

     % For ILAE, add "ILAE" to the MUSC folks
    non_empty_musc = cellfun(@(x,y) ~isempty(x) & contains(y,'MP'),T.ilae_yr1,T.names);
    T.ilae_yr1(non_empty_musc) = cellfun(@(x) sprintf('ILAE %s',x),T.ilae_yr1(non_empty_musc),'UniformOutput',false);

    
    %% Initialize figure
    figure
    set(gcf,'position',[1 1 1400 1000])
        tiledlayout(2,3,"TileSpacing",'tight','padding','tight')
    
    % Prep stats for text
    good_bad = nan(2,2); % engel, ilae; good, bad
    prob_stats = nan(2,7); % engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p
    auc_stats = nan(2,8); %engel, ilae; simple, simpleCI, complicated, complicatedCI, delong p-value, N
    lr_stats = nan(2,7); %engel, ilae; mean left, std left, mean right, std right, df, tstat, p

    %% Get surg nums
    surg_nums = nan(2,2); % left,right;laser,resection
    surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
    %{
    surg_nums = nan(4,3); % all, left, right, bilat; all surg, resection, ablation
    surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
    surg_w_out = surg & cellfun(@(x) ~isempty(x),T.engel_yr1);
    surg_nums(1,1) = sum(surg_w_out);
    surg_nums(1,2) = sum(contains(T.surgery,'Resection')& surg_w_out);
    surg_nums(1,3) = sum(contains(T.surgery,'Laser ablation')& surg_w_out);

    surg_nums(2,1) = sum(surg_w_out & strcmp(T.soz_lats,'left'));
    surg_nums(2,2) = sum(contains(T.surgery,'Resection')& surg_w_out& strcmp(T.soz_lats,'left'));
    surg_nums(2,3) = sum(contains(T.surgery,'Laser ablation')& surg_w_out& strcmp(T.soz_lats,'left'));

    surg_nums(3,1) = sum(surg_w_out & strcmp(T.soz_lats,'right'));
    surg_nums(3,2) = sum(contains(T.surgery,'Resection')& surg_w_out& strcmp(T.soz_lats,'right'));
    surg_nums(3,3) = sum(contains(T.surgery,'Laser ablation')& surg_w_out& strcmp(T.soz_lats,'right'));

    surg_nums(4,1) = sum(surg_w_out & strcmp(T.soz_lats,'bilateral'));
    surg_nums(4,2) = sum(contains(T.surgery,'Resection')& surg_w_out& strcmp(T.soz_lats,'bilateral'));
    surg_nums(4,3) = sum(contains(T.surgery,'Laser ablation')& surg_w_out& strcmp(T.soz_lats,'bilateral'));
    
    


    if ir == 1
        fprintf(fid, [' %d patients (%d, %d, and %d clinically-designated as having left, '...
            'right, and bilateral SOZs, respectively) underwent temporal-targeted resection or ablation. '...
            '%d patients (%d, %d, and %d left, right, and bilateral) underwent resection, and '...
            '%d (%d, %d, and %d left, right, and bilateral) underwent ablation.'],...
            surg_nums(1,1),surg_nums(2,1),surg_nums(3,1),surg_nums(4,1),...
            surg_nums(1,2),surg_nums(2,2),surg_nums(3,2),surg_nums(4,2),...
            surg_nums(1,3),surg_nums(2,3),surg_nums(3,3),surg_nums(4,3));
    end
    %}

    
    % Loop over outcome approaches (Engel vs ILAE), each one gets its own row
    for io = 1:2
        if io == 1
            which_outcome = 'engel';
        elseif io == 2
            which_outcome = 'ilae';
        end
    
    
        % Anonymous function to define good outcome for first plot
        switch which_outcome
            case 'engel'
                good_outcome = @(x) strcmp(x(2),'A') | strcmp(x(2),'B') | strcmp(x(2),'C') | strcmp(x(2),'D');
                which_outcome_text = 'Engel';
                good_outcome_text = 'Engel I';
                bad_outcome_text = 'Engel II+';
            case 'ilae'
                good_outcome = @(x) contains(x,'2') | contains(x,'1');
                which_outcome_text = 'ILAE';
                good_outcome_text = 'ILAE 1-2';
                bad_outcome_text = 'ILAE 3+';
        
        end

       
        %% A and E: Show overall outcomes
        % find those who had surgery
        
        outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
        outcome = T.(outcome_name); 
        empty_outcome = cellfun(@isempty,outcome);
        out_cat = categorical(outcome(surg&~empty_outcome));
        cats = unique(out_cat);
        good = arrayfun(@(x) good_outcome(char(x)),cats);
        
        
        if 0
            oT = T(surg,:);
            table(oT.names,oT.(outcome_name),oT.surgery,oT.surg_lat,oT.soz_lats,oT.surg_loc,oT.soz_locs)
        
        end


        nexttile
        histogram(out_cat,cats)
        hold on
        yl = ylim;
        yl_new = [yl(1) (yl(2)-yl(1))*1.3];
        ybar = (yl(2)-yl(1))*1.1;
        ytext = (yl(2)-yl(1))*1.2;
        ylim(yl_new)
        plot([1 sum(good)],[ybar ybar],'Color',[0.4660, 0.6740, 0.1880]	,'linewidth',2)
        text((1+sum(good))/2,ytext,'Good outcome','fontsize',15,'HorizontalAlignment','center',...
            'color',[0.4660, 0.6740, 0.1880])
        plot([sum(good)+1 length(good)],[ybar ybar],'Color',[0.8500, 0.3250, 0.0980],'linewidth',2)
        text((sum(good)+1+length(good))/2,ytext,'Poor outcome','fontsize',15,'HorizontalAlignment','center',...
            'color',[0.8500, 0.3250, 0.0980])
        plot([(sum(good)+sum(good)+1)/2,(sum(good)+sum(good)+1)/2],ylim, 'k--','linewidth',2)
        ylabel('Number of patients')
        title(sprintf('%s outcome',which_outcome_text))
        set(gca,'fontsize',15)
    
        %% B and E: Compare outcomes for left vs right temporal
        % Just get numerical portion of outcome
        outcome_num =  cellfun(@(x) parse_outcome_num(x,which_outcome),outcome);
        out_text = sprintf('%s numerical portion',which_outcome_text);
        ablation = contains(T.surgery,'ablation');
        resection = contains(T.surgery,'Resection');
        surg_text = cell(length(T.surgery),1);
        surg_text(ablation) = {'Ablation'};
        surg_text(resection) = {'Resection'};
        nexttile
        %{
        out_lat_stats = unpaired_plot_special(outcome_num(strcmp(T.surg_lat,'left')),...
            outcome_num(strcmp(T.surg_lat,'right')),...
            surg_text(strcmp(T.surg_lat,'left')),...
             surg_text(strcmp(T.surg_lat,'right')),...
             {'Left','Right'},out_text,'para');
        %}
        out_lat_stats = unpaired_plot(outcome_num(strcmp(T.surg_lat,'left')),...
            outcome_num(strcmp(T.surg_lat,'right')),...
             {'Left','Right'},out_text,'para');
        set(gca().Children(3),'MarkerSize',10)
        set(gca().Children(4),'MarkerSize',10)
        title('Outcome according to side of surgery')
        
        set(gca,'fontsize',15)
        lr_stats(io,:) = [out_lat_stats.means(1) out_lat_stats.sd(1) out_lat_stats.means(2) out_lat_stats.sd(2),...
            out_lat_stats.df out_lat_stats.tstat out_lat_stats.p];

        has_out = cellfun(@(x) ~isempty(x),T.engel_yr1);
        surg_nums(1,1) = sum(strcmp(T.surg_lat,'left')&ablation&has_out);
        surg_nums(1,2) = sum(strcmp(T.surg_lat,'left')&resection&has_out);
        surg_nums(2,1) = sum(strcmp(T.surg_lat,'right')&ablation&has_out);
        surg_nums(2,2) = sum(strcmp(T.surg_lat,'right')&resection&has_out);
        
        %% C and F: See if modeled probability of concordant laterality is higher for good outcome patients
        % Get models
        left = model.side(1).result;
        right = model.side(2).result;

        left_musc = musc_model.side(1).result;
        right_musc = musc_model.side(2).result;


        % combine hup and musc models
        all_left_names = [left.names;left_musc.names];
        all_right_names = [right.names;right_musc.names];
        
        % confirm that the patients all align with the outcome table
        assert(isequal(all_left_names,all_right_names))
        names = all_left_names;
        assert(isequal(T.names,names))
    
        % 224 did not have surgery, just planned to get it
        T.surgery(strcmp(T.names,'HUP224')) = {'none'};
        
        % Get some basic outcome stuff
        surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
        outcome_bin = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
        good_outcome = strcmp(outcome_bin,'good') & surg == 1;
        bad_outcome = strcmp(outcome_bin,'bad') & surg == 1;
        left_surg = surg & strcmp(T.surg_lat,'left');
        right_surg = surg & strcmp(T.surg_lat,'right');
        npts = length(good_outcome);
    
        
    
        % Make sure no one had both left and right surg
        assert(sum(left_surg&right_surg)==0)
    
        % how many good and bad outcome
        good_bad(io,1) = sum(good_outcome==1);
        good_bad(io,2) = sum(bad_outcome==1);
        
        % Get the model scores (combining hup and musc)
        left_scores = [left.scores;left_musc.scores];
        right_scores = [right.scores;right_musc.scores];
    
        % Get the concordant laterality scores for good and bad outcome
        concordant_lat_scores = nan(npts,1);
        concordant_lat_scores(left_surg == 1) = left_scores(left_surg==1);
        concordant_lat_scores(right_surg == 1) = right_scores(right_surg==1);

        if 0
            table(T.names(surg),outcome_bin(surg),T.surgery(surg),T.surg_lat(surg),...
                left_scores(surg),right_scores(surg),concordant_lat_scores(surg),...
                'VariableNames',{'Name','Outcome','Surgery','Lat','Left score',...
                'Righ score','concordant score'})
        end
    
    
        nexttile
        % concordant lateralty: probability of left for those with left
        % surgery; probability of right for those with right surgery
        %
        stats = unpaired_plot(concordant_lat_scores(good_outcome),concordant_lat_scores(bad_outcome),...
            {good_outcome_text,bad_outcome_text},{'Modeled probability of','concordant laterality'},'para');
        set(gca().Children(3),'MarkerSize',10)
        set(gca().Children(4),'MarkerSize',10)
        %}
        %{
        stats = unpaired_plot_special(concordant_lat_scores(good_outcome),...
            concordant_lat_scores(bad_outcome),...
            surg_text(good_outcome),...
             surg_text(bad_outcome),...
             {good_outcome_text,bad_outcome_text},...
             {'Modeled probability of','concordant laterality'},'para');
        %}
        title({'Surgery-model laterality concordance'})
        xlim([0.5 2.5])
        set(gca,'fontsize',15)
    
        % double check some stuff
        if 0
            table(T.names(left_surg),T.surg_lat(left_surg),T.ilae_yr1(left_surg),good_outcome(left_surg),bad_outcome(left_surg),left_scores(left_surg),right_scores(left_surg))
            table(T.names(right_surg),T.surg_lat(right_surg),T.ilae_yr1(right_surg),good_outcome(right_surg),bad_outcome(right_surg),left_scores(right_surg),right_scores(right_surg))
        end
    
        % engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p
        prob_stats(io,:) = [stats.means(1) stats.sd(1) stats.means(2) stats.sd(2),...
            stats.df stats.tstat stats.p];
    
        % investigating the patients with high predicted concordant
        % laterality but poor outcome
        if 0
            table([T.names(bad_outcome&left_surg);T.names((bad_outcome&right_surg))],[left_scores(bad_outcome&left_surg);right_scores(bad_outcome&right_surg)])
        end
        % HUP138 (0.81) had a high modeled
        % probability of concordant laterality but poor outcome. They
        % had an ablation, and no repeat surgical evaluation. 
    
    
        
    
    end
    
    %% Do text
    if ir == 1
        fprintf(fid,[' %d of %d (%1.1f%%) patients had good one-year Engel outcomes (Engel I), '...
            'and %d of %d (%1.1f%%) had poor Engel outcomes (Engel 2+) (Fig. 4A). '...
            '%d of %d (%1.1f%%) patients had good one-year ILAE outcomes (ILAE 1-2), and '...
            '%d of %d (%1.1f%%) had poor ILAE outcomes (ILAE 3+) (Fig. 4D).'],good_bad(1,1),...
            sum(good_bad(1,:)),good_bad(1,1)/sum(good_bad(1,:))*100,...
            good_bad(1,2),...
            sum(good_bad(1,:)),good_bad(1,2)/sum(good_bad(1,:))*100,...
            good_bad(2,1),...
            sum(good_bad(2,:)),good_bad(2,1)/sum(good_bad(2,:))*100,...
            good_bad(2,2),...
            sum(good_bad(2,:)),good_bad(2,2)/sum(good_bad(2,:))*100);

        fprintf(fid,[' The means of the numerical portions of the outcome scales were '...
            'similar for patients who underwent left versus right-sided surgeries '...
            '(Fig. 4B and 4E; Engel: <i>t</i>(%d) = %1.1f, %s; ILAE: <i>t</i>(%d) = %1.1f, %s).'],...
            lr_stats(1,5),lr_stats(1,6),get_p_html_el(lr_stats(1,7)),...
            lr_stats(2,5),lr_stats(2,6),get_p_html_el(lr_stats(2,7)));

        fprintf(fid,[' Left-sided surgeries were disproportionately ablations (%d ablations vs %d '...
            'resections), and right-sided surgeries were more often resections ('...
            '%d ablations vs %d resections).'],...
            surg_nums(1,1),surg_nums(1,2),surg_nums(2,1),surg_nums(2,2));
        
        
        fprintf(fid,[' We hypothesized that patients with a good surgical outcome would have a '...
            'higher modeled probability of SOZ laterality concordant with the side of surgery. '...
            'We identified the spike rate model corresponding to the side of surgery. '...
            'Mean concordant model probability was significantly higher in patients with good Engel '...
            'outcomes (mean (SD) %1.2f (%1.2f)) than in patients with poor '...
            'Engel outcomes (%1.2f (%1.2f)) (<i>t</i>(%d) = %1.1f, %s) (Fig. 4C), and '...
            'in patients with good ILAE outcomes (%1.2f (%1.2f)) than '...
            'in patients with poor ILAE outcomes (%1.2f (%1.2f)) (<i>t</i>(%d) = %1.1f, %s) (Fig. 4F). '...
            'Together, these results suggest that a model trained to predict the SOZ using spike rate '...
            'asymmetry also predicts surgical outcome.'],...
            prob_stats(1,1),prob_stats(1,2),prob_stats(1,3),prob_stats(1,4),prob_stats(1,5),prob_stats(1,6),...
            get_p_html_el(prob_stats(1,7)),...
            prob_stats(2,1),prob_stats(2,2),prob_stats(2,3),prob_stats(2,4),prob_stats(2,5),prob_stats(2,6),...
            get_p_html_el(prob_stats(2,7)));
        
      
        fprintf(fid,[' Results were similar when we used spikes detected in bipolar and machine references (Fig. S6 and S7).</p>']);
    end
    
    %% Add subtitles
    annotation('textbox',[0 0.9 0.1 0.1],'String','A','LineStyle','none','fontsize',25)
    annotation('textbox',[0.34 0.9 0.1 0.1],'String','B','LineStyle','none','fontsize',25)
    annotation('textbox',[0.67 0.9 0.1 0.1],'String','C','LineStyle','none','fontsize',25)
    annotation('textbox',[0 0.4 0.1 0.1],'String','D','LineStyle','none','fontsize',25)
    annotation('textbox',[0.34 0.4 0.1 0.1],'String','E','LineStyle','none','fontsize',25)
    annotation('textbox',[0.67 0.4 0.1 0.1],'String','F','LineStyle','none','fontsize',25)
   
    
    
    if ir == 1
        print(gcf,[plot_folder,'Fig4'],'-dpng')
    elseif ir == 2
        print(gcf,[plot_folder,'FigS6'],'-dpng')
    elseif ir == 3
        print(gcf,[plot_folder,'FigS7'],'-dpng')
    end
    
    
end



end

