function model_plots

%% Parameters
fsize = 20;
use_model_predict_fcn = 1;
which_refs = {'car','bipolar','machine'};

%% Get file locs
locations = epilepsy_laterality_locs;
plot_folder = locations.el_plots_folder;
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end
addpath(genpath(locations.el_script_folder))

%% Load sozT
sozT = readtable('Manual validation.xlsx','Sheet','SOZ');

%% Load preimplant data
piT = readtable('Manual validation.xlsx','Sheet','Pre-implant data');

% Loop over choice of reference
for ia = 1:length(which_refs) 

    %% Load the intermediate file
    out = load([plot_folder,sprintf('ext_models_%s.mat',which_refs{ia})]);
    out = out.all;

    % If first reference, Loop over all patients vs good outcome HUP
    % patients, otherwise just do all patients
    if ia == 1
        no = 2;
    else
        no = 1;
    end

    for io = 1:no

    
        model = out.approach(io).model;
        nmodels = length(model);
        
        %% Initialize results file
        if ia == 1 && io == 1
            fname = [plot_folder,'results.html'];
            fid = fopen(fname,'a');
            fprintf(fid,'<br><u><i>A classifier incorporating interictal EEG features predicts SOZ laterality</i></u></br>');
        elseif ia == 1 && io == 2
            sfname = [plot_folder,'supplemental_results.html'];
            sfid = fopen(sfname,'a');
            fprintf(sfid,'<br><u><i>SOZ laterality prediction - good outcome analysis</i></u></br>');
        end
        
        
        %% Initialize figure
        figure
        set(gcf,'position',[1 1 1400 1000])
        t = tiledlayout(3,3,"TileSpacing",'tight','padding','tight');
        
        %% A-B Do ROC curves for internal validation
        for iv = 1
            for im = 1:3
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
                    model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},'fontsize',fsize,...
                    'location','southeast')
                title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
                set(gca,'fontsize',fsize)
        
        
            end
        end

        if 0
            % plots for graphical abstract
            iv = 1;
            for im = 1:2
                figure
                ll = plot(model(im).val(iv).side(1).result.X,...
                    model(im).val(iv).side(1).result.Y,'linewidth',3);
                hold on
                lr = plot(model(im).val(iv).side(2).result.X,...
                    model(im).val(iv).side(2).result.Y,'linewidth',3);
                plot([0 1],[0 1],'k--','linewidth',3)
                xlabel('False positive rate')
                ylabel('True positive rate')
                legend([ll,lr],{sprintf('%s: AUC = %1.2f',...
                    'Left model',model(im).val(iv).side(1).result.AUC),...
                    sprintf('%s: AUC = %1.2f',...
                    'Right model',model(im).val(iv).side(2).result.AUC)},'fontsize',25,...
                    'location','southeast')
                %title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
  
                set(gca,'fontsize',25)
                print(gcf,[plot_folder,sprintf('abstract_%d',im)],'-depsc')
                close(gcf)
            end
        end
        
        if ia == 1 && io == 1
            fprintf(fid,['We tested whether interictal features predict SOZ laterality in unseen patients. '...
                'The AUCs of the ROC of the left- and right-sided internal cross-validation models trained on all features were %1.2f '...
                'and %1.2f, respectively (Fig. 4A). A model '...
                'trained on only spike rates (restricted to %s reference and sleep '...
                'in order to have a single feature) achieved higher AUCs (%1.2f '...
                'and %1.2f for the left and right models, respectively (Fig. 4B)). Finally, a model trained only on '...
                'binary spike rates indicating whether there were more spikes on the left or the right '...
                'performed poorly (AUC of %1.2f and %1.2f'...
                ', respectively (Fig. 4C)). '],...
                model(1).val(1).side(1).result.AUC,model(1).val(1).side(2).result.AUC,...
                upper(which_refs{ia}),...
                model(2).val(1).side(1).result.AUC,model(2).val(1).side(2).result.AUC,...
                model(3).val(1).side(1).result.AUC,model(3).val(1).side(2).result.AUC);
        elseif ia == 1 && io == 2
            fprintf(sfid,['We again tested the ability of interictal features to predict SOZ laterality in unseen patients, '...
                'now restricting the unilateral patients in the HUP dataset to be those with Engel 1 surgical outcomes. '...
                'The AUCs of the ROC of the left- and right-sided internal cross-validation models trained on all features were %1.2f '...
                'and %1.2f, respectively (Fig. S5A). A model '...
                'trained on only spike rates (%s reference) achieved better AUCs (%1.2f '...
                'and %1.2f for the left and right models, respectively (Fig. S5B)). Finally, a model trained only on '...
                'binary spike rates indicating whether there were more spikes on the left or the right '...
                'performed poorly (AUC of %1.2f and %1.2f'...
                ', respectively (Fig. S5C)). '],...
                model(1).val(1).side(1).result.AUC,model(1).val(1).side(2).result.AUC,...
                upper(which_refs{ia}),...
                model(2).val(1).side(1).result.AUC,model(2).val(1).side(2).result.AUC,...
                model(3).val(1).side(1).result.AUC,model(3).val(1).side(2).result.AUC);
        end
        
        
        
        %% D: ROC for spikes, external validation
        if 1
        
            %% Stem plot of most important features
            %rm_sleep_text = @(x) strrep(x,' sleep','');
            rm_sleep_text = @(x) x;
            shorten_bi_text = @(x) strrep(x,'bipolar','bi');
            shorten_machine_text = @(x) strrep(x,'machine','mac');
            all_shorten = @(x) rm_sleep_text(shorten_machine_text(shorten_bi_text(x)));
        
            left = model(1).val(2).side(1).result;
            right = model(1).val(2).side(2).result;
            n_to_plot = 15;
            left_coefs = left.coefs(1:n_to_plot);
            left_names = left.sorted_features(1:n_to_plot);
            right_coefs = right.coefs(1:n_to_plot);
            right_names = right.sorted_features(1:n_to_plot);
            
            tt = tiledlayout(t,1,1);
            tt.Layout.Tile = 4;
            tt.Layout.TileSpan = [1 1];
            ax1 = axes(tt);
            
            
            pl = plot(ax1,1:n_to_plot,abs(left_coefs),'o','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
            ax1.XAxisLocation = 'top';
            ax1.YAxisLocation = 'left';
            ax1.XColor = [0 0.4470 0.7410];
            ax1.YColor = [0 0.4470 0.7410];
            ax1.Box = 'off';
            ax1.YLim = [min([min(abs(left_coefs)),min(abs(right_coefs))]),...
                max([max(abs(left_coefs)),max(abs(right_coefs))])];
            ax1.XTick = 1:n_to_plot;
            ax1.XTickLabel = cellfun(all_shorten,cellfun(@greek_letters_plots,left_names,'uniformoutput',false),...
                'uniformoutput',false);
            
            ax2 = axes(tt);
            pr = plot(ax2,1:n_to_plot,abs(right_coefs),'+','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
            ax2.XAxisLocation = 'bottom';
            ax2.YAxisLocation = 'right';
            ax2.XColor = [0.8500 0.3250 0.0980];
            ax2.YColor = [0.8500 0.3250 0.0980];
            ax2.Box = 'off';
            ax2.YLim = [min([min(abs(left_coefs)),min(abs(right_coefs))]),...
                max([max(abs(left_coefs)),max(abs(right_coefs))])];
            ax2.XTick = 1:n_to_plot;
            ax2.XTickLabel = cellfun(all_shorten,cellfun(@greek_letters_plots,right_names,'uniformoutput',false),...
                'uniformoutput',false);
            ax2.Color = 'none';
            ax2.Box = 'off';
            
            legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',fsize)
            set(ax1,'fontsize',fsize); set(ax2,'fontsize',fsize)
            ylabel('|Estimated model coefficient|','color','k','fontsize',15)
            
            if ia == 1 && io == 1
                fprintf(fid,['Spike rates, spike timing, and bandpower were the most important '...
                    'features for both the left- and the right-sided models (Fig. 4D).</p>']);
            elseif ia ==1 && io ==2
                fprintf(sfid,['Spike rates, spike timing, and bandpower were the most important '...
                    'features for both the left- and the right-sided models (Fig. 4D).</p>']);
            end
        
        else
        
            iv = 2;
            im = 2;
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
                model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},...
                'location','southeast','fontsize',15)
            title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
            set(gca,'fontsize',15)
            
            fprintf(fid,['<p>Finally, we tested how the spike-only models performed in an '...
                'external dataset. We trained the models using only spike rate data '...
                'on all HUP patients, and tested the models on patients from MUSC. '...
                'The model AUCs on the MUSC patients were %1.2f and %1.2f for '...
                'the left-sided and right-sided models, respectively (Fig. 5G).'],...
                model(2).val(2).side(1).result.AUC,model(2).val(2).side(2).result.AUC);
        
        end
        
        %% E-F Confusion matrix (spikes only, internal cross-validation)
        curr = model(2).val(1);
        both_bal_acc = nan(2,1);
        both_opt_thresh = nan(2,1);
        both_opt = nan(2,2);
        for is = 1:2
            scores = curr.side(is).result.scores;
            class = curr.side(is).result.class;
        
            classes = curr.side(is).result.unique_classes;
            nclasses = length(classes);
            pos_class = curr.side(is).result.pos_class;
            neg_class = classes(~strcmp(classes,pos_class));
        
            % put the positive class on top
            if ~strcmp(curr.side(is).result.pos_class,classes{1})
                classes = classes([2 1]);
            end
        
        
            if use_model_predict_fcn 
                pred = curr.side(is).result.all_pred;
        
            else % find optimal point myself
                % get cost for finding optimal ROC point
                cost = [0 sum(strcmp(class,neg_class));sum(strcmp(class,pos_class)) 0];
            
                % Get optimal ROC point
                
                [X,Y,T,~,opt] = perfcurve(class,scores,pos_class,'cost',cost);
                %[X,Y,T,~,opt] = perfcurve(class,scores,pos_class);
                opt_thresh = T((X==opt(1))&(Y==opt(2)));
                %opt_thresh = my_opt(X,Y,T);
                opt_thresh = 0.5;
                both_opt_thresh(is) = opt_thresh;
                both_opt(is,:) = opt;
                
                pred = cell(length(class),1);
                pred(scores >= opt_thresh) = {pos_class};
                pred(scores < opt_thresh) = neg_class;
            end
        
            % Rebuild C
            positive = strcmp(class,pos_class);
            negative = strcmp(class,neg_class);
            pred_positive = strcmp(pred,pos_class);
            pred_negative = strcmp(pred,neg_class);
            C(1,1) = sum(positive & pred_positive);
            C(1,2) = sum(positive & pred_negative);
            C(2,1) = sum(negative & pred_positive);
            C(2,2) = sum(negative & pred_negative);
            
            % Calculate accuracy
            accuracy = sum(diag(C))/sum(C(:));
            % Balanced accuracy is the average across all classes of the number of 
            % data accurately predicted belonging to class m divided by the number of
            % data belonging to class m
            recall = nan(nclasses,1);
            for i = 1:nclasses
                tp = C(i,i);
                fn = sum(C(i,~ismember(1:nclasses,i))); 
                recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
            end
            balanced_accuracy = mean(recall);
            both_bal_acc(is) = balanced_accuracy;
            
            % Plot
            nexttile(t)
            % Map numbers onto 0 to 1
            new_numbers = map_numbers_onto_range(C,[1 0]);
            Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
            D = diag(new_numbers);
            Dcolor = [repmat(D,1,2),ones(length(D),1)];
            Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
            imagesc(Ccolor)
            
            % replace classnames
            pretty_name = classes;
            pretty_name = strrep(pretty_name,'left','Left');
            pretty_name = strrep(pretty_name,'right','Right');
            pretty_name = strrep(pretty_name,'br','Right/bilateral');
            pretty_name = strrep(pretty_name,'bl','Left/bilateral');
            xticks(1:nclasses)
            xticklabels((pretty_name))
            yticks(1:nclasses)
            yticklabels((pretty_name))
            ytickangle(90)
            xlabel('Predicted')
            ylabel('True')
            hold on
            for i = 1:nclasses
                for j = 1:nclasses
                    text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
                end
            end
            title(sprintf('Spikes (HUP cross-validation)\nBalanced accuracy %1.1f%%',balanced_accuracy*100))
            set(gca,'fontsize',fsize)
        
        end
        
        if ia == 1 && io == 1
            fprintf(fid,['<p>We further probed the accuracy of the spike-rate only model. '...
                'Confusion matrices for the left- and right-sided models at the optimal '...
                'operating points '...
                'are shown in Figs. 4E and 4F. The balanced accuracy '...
                'was %1.1f%% for the model predicting left vs. right/bilateral SOZ, '...
                'and %1.1f%% for the model predicting right vs. left/bilateral SOZ. '],...
                both_bal_acc(1)*100,both_bal_acc(2)*100);
        elseif ia == 1 && io == 2
            fprintf(sfid,['<p>We further probed the accuracy of the spike-rate only model. '...
                'Confusion matrices for the left- and right-sided models at the optimal '...
                'operating points '...
                'are shown in Figs. S3E and S3F. The balanced accuracy '...
                'was %1.1f%% for the model predicting left vs. right/bilateral SOZ, '...
                'and %1.1f%% for the model predicting right vs. left/bilateral SOZ. '],...
                both_bal_acc(1)*100,both_bal_acc(2)*100);
        end
        
        %% G: Subsampling plots, internal validation
        if 1
        sub = out.cv_ss; % cross validation
        
        nexttile(t)
        
        durations = [1 5 10 20 30];
        
        ndurs = length(durations);
        nsamples = size(sub,4);
        
        
        curr_ss = 2; % just do sleep
        auc_l = squeeze(sub(curr_ss,1,:,:));
        auc_r = squeeze(sub(curr_ss,2,:,:));
        
        median_l = nanmedian(auc_l,2);
        median_r = nanmedian(auc_r,2);
        P_l_25 = prctile(auc_l,[25],2);
        P_r_25 = prctile(auc_r,[25],2);
        P_l_75 = prctile(auc_l,[75],2);
        P_r_75 = prctile(auc_r,[75],2);
        
        U_l = P_l_75-median_l;
        U_r = P_r_75-median_r;
        L_l = median_l - P_l_25;
        L_r = median_r - P_r_25;
        
        % Plot it
        
        %
        el = shaded_error_bars_fc_el(1:ndurs,median_l,[P_l_75';P_l_25'],[0, 0.4470, 0.7410]);
        hold on
        er = shaded_error_bars_fc_el(1:ndurs,median_r,[P_r_75';P_r_25'],[0.8500, 0.3250, 0.0980]);
        
        
        errorbar(1:ndurs,median_l,L_l,U_l,'o','color',[0, 0.4470, 0.7410],...
            'LineWidth',2,'MarkerSize',10);
        hold on
        errorbar(1:ndurs,median_r,L_r,U_r,'o','color',[0.8500, 0.3250, 0.0980],...
            'LineWidth',2,'MarkerSize',10);
        %}
        
        ylim([0.4 0.9])
        
        legend([el,er],{'Left vs right/bilateral','Right vs left/bilateral'},...
            'location','southeast','fontsize',fsize)
        xticks(1:ndurs)
        xticklabels(arrayfun(@(x) sprintf('%d min',x),durations,'uniformoutput',false))
        ylabel('Median (IQR) AUC')
        title(sprintf('Spike model accuracy by duration\n(HUP cross-validation)'))
        set(gca,'fontsize',fsize)
        
        if ia == 1 && io == 1
            fprintf(fid,['Model accuracies rise quickly with duration sampled, achieving '...
                'an accuracy similar to the full-duration models with 5 minutes of sampling (Fig. 4G).']);
        elseif ia == 1 && io == 2
            fprintf(sfid,['Model accuracies rise quickly with duration sampled, achieving '...
                'an accuracy similar to the full-duration models with 5 minutes of sampling (Fig. S5G).']);
        end
        
        end
        
        
        %% H-I: Confusion matrix (spikes only, external validation)
        curr = model(2).val(2); % spike  model, external validation
        both_bal_acc = nan(2,1);
        for is = 1:2
            scores = curr.side(is).result.scores;
            class = curr.side(is).result.class;
        
            classes = curr.side(is).result.unique_classes;
            nclasses = length(classes);
            pos_class = curr.side(is).result.pos_class;
            neg_class = classes(~strcmp(classes,pos_class));
        
            % put the positive class on top
            if ~strcmp(curr.side(is).result.pos_class,classes{1})
                classes = classes([2 1]);
            end
        
            if use_model_predict_fcn 
                pred = curr.side(is).result.all_pred;
        
            else
                % Get optimal ROC point
                %[X,Y,T,~,opt] = perfcurve(class,scores,pos_class);
                %opt_thresh = T((X==opt(1))&(Y==opt(2)));
                %opt_thresh = my_opt(X,Y,T);
                % use the optimal threshold from the internal cross validated model
                opt_thresh = both_opt_thresh(is);
                
                pred = cell(length(class),1);
                pred(scores >= opt_thresh) = {pos_class};
                pred(scores < opt_thresh) = neg_class;
            end
    
        
            % Rebuild C
            positive = strcmp(class,pos_class);
            negative = strcmp(class,neg_class);
            pred_positive = strcmp(pred,pos_class);
            pred_negative = strcmp(pred,neg_class);
            C(1,1) = sum(positive & pred_positive);
            C(1,2) = sum(positive & pred_negative);
            C(2,1) = sum(negative & pred_positive);
            C(2,2) = sum(negative & pred_negative);
            
            % Calculate accuracy
            accuracy = sum(diag(C))/sum(C(:));
            % Balanced accuracy is the average across all classes of the number of 
            % data accurately predicted belonging to class m divided by the number of
            % data belonging to class m
            recall = nan(nclasses,1);
            for i = 1:nclasses
                tp = C(i,i);
                fn = sum(C(i,~ismember(1:nclasses,i))); 
                recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
            end
            
            balanced_accuracy = mean(recall);
            both_bal_acc(is) = balanced_accuracy;
            
            % Plot
            nexttile(t)
            % Map numbers onto 0 to 1
            new_numbers = map_numbers_onto_range(C,[1 0]);
            Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
            D = diag(new_numbers);
            Dcolor = [repmat(D,1,2),ones(length(D),1)];
            Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
            imagesc(Ccolor)
            
            % replace classnames
            pretty_name = classes;
            pretty_name = strrep(pretty_name,'left','Left');
            pretty_name = strrep(pretty_name,'right','Right');
            pretty_name = strrep(pretty_name,'br','Right/bilateral');
            pretty_name = strrep(pretty_name,'bl','Left/bilateral');
            xticks(1:nclasses)
            xticklabels((pretty_name))
            yticks(1:nclasses)
            yticklabels((pretty_name))
            ytickangle(90)
            xlabel('Predicted')
            ylabel('True')
            hold on
            for i = 1:nclasses
                for j = 1:nclasses
                    text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
                end
            end
            title(sprintf('Spikes (MUSC external validation)\nBalanced accuracy %1.1f%%',balanced_accuracy*100))
            set(gca,'fontsize',fsize)
        
        end
        
        if ia == 1 && io == 1
            fprintf(fid,[' Finally, we tested how the spike-only models performed in the '...
                    'external MUSC dataset.']);
            fprintf(fid,[' The balanced accuracies were %1.1f%% and %1.1f%% '...
                'for the left-sided and right-sided models, respectively (Fig. 4H and 4I).</p>'],both_bal_acc(1)*100,...
                both_bal_acc(2)*100);
            
            fprintf(fid,['<p>These results suggest that models '...
                'using only spike rate asymmetry accurately distinguished left from right or bilateral SOZs '...
                'in both internal cross-validation and in a separate institution''s test dataset. However, '...
                'although right-sided SOZs could be distinguished from left/bilateral SOZs in the external '...
                'dataset set, they were not well-classified in the internal validation dataset.']);

            fprintf(fid,[' Results were similar, although with higher AUCs '...
                'across all models, when we restricted analysis of unilateral HUP patients to be those with '...
                'Engel 1 surgical outcomes to ' ...
                'build and internally validate the SOZ laterality classifier (Fig. S4).']);
    
            fprintf(fid,[' Results were also similar when we used spikes detected in bipolar and machine references to '...
                'build the SOZ laterality classifier (Fig. S5 and S6).']);
        elseif ia == 1 && io == 2

            fprintf(sfid,[' Finally, we tested how the spike-only models performed in the '...
                    'external MUSC dataset.']);
            fprintf(sfid,[' The balanced accuracies were %1.1f%% and %1.1f%% '...
                'for the left-sided and right-sided models, respectively (Fig. 4H and 4I).</p>'],both_bal_acc(1)*100,...
                both_bal_acc(2)*100);
            
            fprintf(sfid,['<p>These results were overall similar as those when we include '...
                'all patients regardless of surgical outcome, although higher AUCs were '...
                'achieved across all models when we restrict analysis of unilateral '...
                'patients to those with good outcomes.</p>']);
        end
            
        
        %% Add subtitles
        annotation('textbox',[0 0.905 0.1 0.1],'String','A','LineStyle','none','fontsize',20)
        annotation('textbox',[0.36 0.905 0.1 0.1],'String','B','LineStyle','none','fontsize',20)
        annotation('textbox',[0.70 0.905 0.1 0.1],'String','C','LineStyle','none','fontsize',20)
        annotation('textbox',[0 0.56 0.1 0.1],'String','D','LineStyle','none','fontsize',20)
        annotation('textbox',[0.36 0.56 0.1 0.1],'String','E','LineStyle','none','fontsize',20)
        annotation('textbox',[0.70 0.56 0.1 0.1],'String','F','LineStyle','none','fontsize',20)
        annotation('textbox',[0 0.22 0.1 0.1],'String','G','LineStyle','none','fontsize',20)
        annotation('textbox',[0.36 0.22 0.1 0.1],'String','H','LineStyle','none','fontsize',20)
        annotation('textbox',[0.70 0.22 0.1 0.1],'String','I','LineStyle','none','fontsize',20)
        
        
        
        
        if ia == 1 && io == 1
            %print(gcf,[plot_folder,'Fig3'],'-dpng')
            print(gcf,[plot_folder,'Fig4'],'-dtiff')
        elseif ia == 1 && io == 2
            %print(gcf,[plot_folder,'FigS3'],'-dpng')
            print(gcf,[plot_folder,'FigS4'],'-dtiff')
        elseif ia == 2
            %print(gcf,[plot_folder,'FigS4'],'-dpng')
            print(gcf,[plot_folder,'FigS5'],'-dtiff')
        elseif ia == 3
            %print(gcf,[plot_folder,'FigS5'],'-dpng')
            print(gcf,[plot_folder,'FigS6'],'-dtiff')
        end


        %% Examining error sources
        if ia ==1 && io == 1
            error_stats = nan(2,4);
            mri_error_stats = nan(2,4);

            % loop over left and right
            for is = 1:2

                % Find the label of correct vs incorrect in both hup and musc
                musc_agree = cellfun(@(x,y) strcmp(x,y),...
                    model(2).val(2).side(is).result.class,...
                    model(2).val(2).side(is).result.all_pred);
                musc_names = model(2).val(2).side(is).result.names;

                hup_agree = cellfun(@(x,y) strcmp(x,y),...
                    model(2).val(1).side(is).result.class,...
                    model(2).val(1).side(is).result.all_pred);
                hup_names = model(2).val(1).side(is).result.names;

                all_agree = [musc_agree;hup_agree];
                all_names = [musc_names;hup_names];

                if 0
                    table(all_names,all_agree)
                end

                % Find the corresponding soz loc
                soz_loc = cell(length(all_names),1);
                for in = 1:length(all_names)
                    % find the corresponding row
                    r = strcmp(all_names{in},sozT.name);
                    if sum(r) ~= 1, error('what'); end
                    soz_loc{in} = sozT.region{r};
                end
                soz_spec = cell(length(all_names),1);
                soz_spec(contains(soz_loc,'mesial')) = {'mesial temporal'};
                soz_spec(contains(soz_loc,'neocortical')) = {'temporal neocortical'};
                soz_spec(cellfun(@isempty,soz_spec)) = {'other'};
                nother = sum(sum(strcmp(soz_spec,'other')));

                % Find the corresponding lesional status
                mri_loc_table = piT.MRILesionalLocalization_temporal_Frontal_Other_Multifocal_Broad;
                mri_lat_table = piT.MRILesionalLaterality_left_Right_Bilateral_Broad_NA__NAMeansNon;
                mri_name = piT.name;
                mri_loc = cell(length(all_names),1);
                mri_lat = cell(length(all_names),1);
                for in = 1:length(all_names)
                    % find the corresponding row
                    r = strcmp(all_names{in},mri_name);
                    if sum(r) > 1, error('what'); end

                    % don't have musc yet
                    if sum(r) == 0, error('what'); end
                    mri_loc{in} = mri_loc_table{r};
                    mri_lat{in} = mri_lat_table{r};
                end
                mri_spec = cell(length(all_names),1);
                mri_loc(cellfun(@(x) isempty(x),mri_loc)) = {'na'};
                mri_lat(cellfun(@(x) isempty(x),mri_lat)) = {'na'};
                %{
                mri_spec(contains(mri_loc,'temporal','IgnoreCase',true) & ...
                    (strcmpi(mri_lat,'left') | strcmpi(mri_lat,'right'))) = {'unilateral temporal'};
                mri_spec(~contains(mri_loc,'temporal','IgnoreCase',true) | ...
                    ~(strcmpi(mri_lat,'left') | strcmpi(mri_lat,'right'))) = {'other'};
                %}
                mri_spec(contains(mri_loc,'temporal','IgnoreCase',true)) = {'temporal'};
                mri_spec(~contains(mri_loc,'temporal','IgnoreCase',true)) = {'other'};

                % make a table
                errorT = table(all_names,all_agree,soz_spec);
                mri_errorT = table(all_names,all_agree,mri_spec);

                % remove the others
                errorT(strcmp(errorT.soz_spec,'other'),:) = [];

       
                % 2x2 table
                [tbl2x2,~,~,labels] = crosstab(errorT.all_agree,errorT.soz_spec);
                assert(strcmp(labels{1,1},'0'))
                assert(strcmp(labels{1,2},'mesial temporal'))
                [h,p,stats] = fishertest(tbl2x2);
                error_stats(is,:) = [stats.OddsRatio, stats.ConfidenceInterval, p];

                % mri 2x2 table
                [tbl2x2,~,~,labels] = crosstab(mri_errorT.all_agree,mri_errorT.mri_spec);
                assert(strcmp(labels{1,1},'0'))
                assert(strcmp(labels{1,2},'temporal'))
                [h,p,stats] = fishertest(tbl2x2);
                mri_error_stats(is,:) = [stats.OddsRatio, stats.ConfidenceInterval, p];

                if 0
                    heatmap(mri_errorT,'all_agree','mri_spec')
                end

            end

            fprintf(fid,[' There was no significant association between'...
                ' model accuracy and mesial temporal (N = %d) vs. temporal neocortical (N = %d) '...
                'localization for either the left- or right-sided model '...
                '(left-sided model odds-ratio: %1.1f (95%% CI %1.1f-%1.1f), '...
                '<i>p</i> = %1.2f; right-sided model: %1.1f (%1.1f-%1.1f), '...
                '<i>p</i> = %1.2f).'],...
                sum(strcmp(soz_spec,'mesial temporal')),sum(strcmp(soz_spec,'temporal neocortical')),...
                error_stats(1,1),error_stats(1,2),error_stats(1,3),error_stats(1,4),...
                error_stats(2,1),error_stats(2,2),error_stats(2,3),error_stats(2,4));

            fprintf(fid,[' Similarly, there was no significant difference in model accuracy between patients '...
                'with temporal lobe lesions on MRI (N = %d) and '...
                'patients whose MRI had other lesions or no lesion (N = %d) '...
                '(left-sided model odds-ratio: %1.1f (95%% CI %1.1f-%1.1f), '...
                '<i>p</i> = %1.2f; right-sided model: %1.1f (%1.1f-%1.1f), '...
                '<i>p</i> = %1.2f).</p>'],...
                sum(strcmp(mri_spec,'temporal')),sum(strcmp(mri_spec,'other')),...
                mri_error_stats(1,1),mri_error_stats(1,2),mri_error_stats(1,3),mri_error_stats(1,4),...
                mri_error_stats(2,1),mri_error_stats(2,2),mri_error_stats(2,3),mri_error_stats(2,4));
        end
    
    end
end

end