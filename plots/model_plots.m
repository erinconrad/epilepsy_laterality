function model_plots

%% Parameters
use_model_predict_fcn = 1;
which_refs = {'car','bipolar','machine'};

%% Get file locs
locations = epilepsy_laterality_locs;
plot_folder = locations.el_plots_folder;
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end



% Loop over choice of reference
for ia = 1:length(which_refs) 

    %% Load the intermediate file
    out = load([plot_folder,sprintf('ext_models_%s.mat',which_refs{ia})]);
    out = out.all;

    % Loop over all patients vs good outcome HUP patients
    for io = 1:2

    
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
        set(gcf,'position',[1 1 1200 1400])
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
                    model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},'fontsize',15,...
                    'location','southeast')
                title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
                set(gca,'fontsize',15)
        
        
            end
        end
        
        if ia == 1 && io == 1
            fprintf(fid,['We tested whether interictal features predict SOZ laterality in unseen patients. '...
                'The AUCs of the ROC of the left- and right-sided internal cross-validation models trained on all features were %1.2f '...
                'and %1.2f, respectively (Fig. 3A). A model '...
                'trained on only spike rates (restricted to %s reference and sleep '...
                'in order to have a single feature) achieved higher AUCs (%1.2f '...
                'and %1.2f for the left and right models, respectively (Fig. 3B)). Finally, a model trained only on '...
                'binary spike rates indicating whether there were more spikes on the left or the right '...
                'performed poorly (AUC of %1.2f and %1.2f'...
                ', respectively (Fig. 3C)). '],...
                model(1).val(1).side(1).result.AUC,model(1).val(1).side(2).result.AUC,...
                upper(which_refs{ia}),...
                model(2).val(1).side(1).result.AUC,model(2).val(1).side(2).result.AUC,...
                model(3).val(1).side(1).result.AUC,model(3).val(1).side(2).result.AUC);
        elseif ia == 1 && io == 2
            fprintf(sfid,['We again tested the ability of interictal features predict SOZ laterality in unseen patients, '...
                'now restricting the unilateral patients in the HUP dataset to be those with Engel 1 surgical outcomes. '...
                'The AUCs of the ROC of the left- and right-sided internal cross-validation models trained on all features were %1.2f '...
                'and %1.2f, respectively (Fig. S3A). A model '...
                'trained on only spike rates (%s reference) achieved better AUCs (%1.2f '...
                'and %1.2f for the left and right models, respectively (Fig. S3B)). Finally, a model trained only on '...
                'binary spike rates indicating whether there were more spikes on the left or the right '...
                'performed poorly (AUC of %1.2f and %1.2f'...
                ', respectively (Fig. S3C)). '],...
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
            
            legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',15)
            set(ax1,'fontsize',15); set(ax2,'fontsize',15)
            ylabel('|Estimated model coefficient|','color','k','fontsize',15)
            
            if ia == 1 && io == 1
                fprintf(fid,['Spike rates, spike timing, and bandpower were the most important '...
                    'features for both the left- and the right-sided models (Fig. 3D).</p>']);
            elseif ia ==1 && io ==2
                fprintf(sfid,['Spike rates, spike timing, and bandpower were the most important '...
                    'features for both the left- and the right-sided models (Fig. 3D).</p>']);
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
            set(gca,'fontsize',15)
        
        end
        
        if ia == 1 && io == 1
            fprintf(fid,['<p>We further probed the accuracy of the spike-rate only model. '...
                'Confusion matrices for the left- and right-sided models at the optimal '...
                'operating points '...
                'are shown in Figs. 3E and 3F. The balanced accuracy '...
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
            'location','southeast','fontsize',15)
        xticks(1:ndurs)
        xticklabels(arrayfun(@(x) sprintf('%d min',x),durations,'uniformoutput',false))
        ylabel('Median (IQR) AUC')
        title(sprintf('Spike model accuracy by duration\n(HUP cross-validation)'))
        set(gca,'fontsize',15)
        
        if ia == 1 && io == 1
            fprintf(fid,['Model accuracies rise quickly with duration sampled, achieving '...
                'an accuracy similar to the full-duration models with 5 minutes of sampling (Fig. 3G).']);
        elseif ia == 1 && io == 2
            fprintf(sfid,['Model accuracies rise quickly with duration sampled, achieving '...
                'an accuracy similar to the full-duration models with 5 minutes of sampling (Fig. S3G).']);
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
            set(gca,'fontsize',15)
        
        end
        
        if ia == 1 && io == 1
            fprintf(fid,[' Finally, we tested how the spike-only models performed in the '...
                    'external MUSC dataset.']);
            fprintf(fid,[' The balanced accuracies were %1.1f%% and %1.1f%% '...
                'for the left-sided and right-sided models, respectively (Fig. 3H and 3I).</p>'],both_bal_acc(1)*100,...
                both_bal_acc(2)*100);
            
            fprintf(fid,['<p>These results suggest that models '...
                'using only spike rate asymmetry accurately distinguished left from right or bilateral SOZs '...
                'in both internal cross-validation and in a separate institution''s test dataset. However, '...
                'although right-sided SOZs could be distinguished from left/bilateral SOZs in the external '...
                'dataset set, they were not well-classified in the internal validation dataset.']);

            fprintf(fid,[' Results were similar, although with higher AUCs '...
                'across all models, when we restricted analysis of unilateral HUP patients to be those with '...
                'Engel 1 surgical outcomes to ' ...
                'build and internally validate the SOZ laterality classifier (Fig. S3).']);
    
            fprintf(fid,[' Results were also similar when we used spikes detected in bipolar and machine references to '...
                'build the SOZ laterality classifier (Fig. S4 and S5).</p>']);
        elseif ia == 1 && io == 2

            fprintf(sfid,[' Finally, we tested how the spike-only models performed in the '...
                    'external MUSC dataset.']);
            fprintf(sfid,[' The balanced accuracies were %1.1f%% and %1.1f%% '...
                'for the left-sided and right-sided models, respectively (Fig. 3H and 3I).</p>'],both_bal_acc(1)*100,...
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
            print(gcf,[plot_folder,'Fig3'],'-dpng')
        elseif ia == 1 && io == 2
            print(gcf,[plot_folder,'FigS3'],'-dpng')
        elseif ia == 2
            print(gcf,[plot_folder,'FigS4'],'-dpng')
        elseif ia == 3
            print(gcf,[plot_folder,'FigS5'],'-dpng')
        end
    
    end
end

end