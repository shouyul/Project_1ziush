%clear all;
%close all;
config_PIONEER_VEP;
study_config.epochs.limits_wdw = [0,10];
warning('Hardcoding epoch limit window, actually depends on subject (was 30, now 10)');

switch lower(user)
    case 'sl'
        % stop MATLAB from stealing focus
        set(groot, 'DefaultFigureVisible', 'off'); % doesn't work here, set figure visible off below, doesn't work either
end

%% Single subject analysis
params.model = study_config.class.model;

for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    %subject_ind = 2;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    input_filepath = N.searchFolder_3arch_rej_ICcats;
    
    if ~exist(fullfile(study_config.figures_folder,'Classifiers'), 'dir')
        mkdir(fullfile(study_config.figures_folder,'Classifiers'));
    end
    
    [params.saveDataFolder, params.saveFigFolder] = makeClassifierArchitecture(input_filepath, study_config);
    
    if strcmp(study_config.epochs.event, 'TrialStart')
        %% Load data
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath); % For chanlocs
        chanlocs = EEG.chanlocs(~strcmp({EEG.chanlocs.labels},study_config.eog_channels) &...
            ~contains({EEG.chanlocs.type}, 'MOCAP'));
        %clear EEG;
        
        switch study_config.feat.type
            case 'psd'
                opts = study_config.psd;
                opts.SR = EEG.srate;
                folderName = makePSDFolderName(opts);
                input_filepath_feats = fullfile(input_filepath, folderName);
            case 'amp'
                opts = study_config.amps;
                input_filepath_feats = input_filepath;
        end
        
        opts.event = study_config.epochs.event;
        delims = strfind(N.epochedFile,'_');
        opts.epochWdW = N.epochedFile(delims(end-1)+1:delims(end)-1);
        fileName = makePSDFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo
        fileName = makePSDFileName(study_config.feat.type, subject, opts);
        load(fullfile(input_filepath_feats, fileName)); % Creates PSD_norm or Amps
        
        delims = strfind(fileName, '_');
        params.suffix = fileName(delims(4)+1:end-4);
        
        [~, chan_labels, dim2_labels] = select_features(PSD_norm, chanlocs,...
            study_config.psd.FoI, study_config.feat);
        clear PSD_norm
        
        
        %% Create labels
        switch study_config.class.condAnalysis
            case ''
                if contains(study_config.class.contrast, 'Pairwise')
                    all_Labels = zeros(size(TrialsInfo,1),1);
                    all_Labels(strcmp(TrialsInfo.TrialType, 'Disc')) = 1;
                            if study_config.subjects(study_config.current_subject).date < datetime(2023,7,30)
                                all_Labels(strcmp(TrialsInfo.TrialType, '3Rings')) = 2;
                                all_Labels(strcmp(TrialsInfo.TrialType, '5Rings')) = 3;
                                all_Labels(strcmp(TrialsInfo.TrialType, '45Bar')) = 4;
                                all_Labels(strcmp(TrialsInfo.TrialType, '-45Bar')) = 5;
                                all_classes = {'Control','Disc','3Rings','5Rings','45Bar','-45Bar'};
                            else
                                all_Labels(strcmp(TrialsInfo.TrialType, '0Bar')) = 2;
                                all_Labels(strcmp(TrialsInfo.TrialType, '90Bar')) = 3;                                
                                all_classes = {'Control','Disc','0Bar','90Bar'};
                            end
                    n_classes = numel(all_classes);
                    
                    NrepRandom = 1000;
                    n_pairs = nchoosek(n_classes,2);
                    RandomErr = nan(NrepRandom, n_classes, n_classes);
                    %TestingErrAll = nan(length(study_config.class.feat2test), study_config.class.nbfolds, 2*numel(conds));
                    TestingErrAve = nan(length(study_config.class.feat2test), n_classes, n_classes);
                    TrainingErrAve = nan(length(study_config.class.feat2test), n_classes, n_classes);
                    TestingErrSE = nan(length(study_config.class.feat2test), n_classes, n_classes);
                    TrainingErrSE = nan(length(study_config.class.feat2test), n_classes, n_classes);
                    
                    optErrors = zeros(2,n_classes,n_classes);
                    optSEs = zeros(2,n_classes,n_classes);
                    %             if fold < length(nb_folds)+1
                    %                 All_samples = zeros(size(TestingErrAve,2),nb_folds(fold));
                    %             end
                    optNfeats = zeros(n_classes,n_classes);
                    
                    %plot_params.error = true;
                    %plot_params.errorBal = false;
                    %plot_params.CM = true;
                    
                    n_chans = numel(chan_labels);
                    n_freqs = numel(dim2_labels);
                    feat_occurences = zeros(n_chans, n_freqs, n_classes, n_classes);
                    
                    for c1 = 1:n_classes-1
                        for c2 = c1+1:n_classes
                            params.name = sprintf('%s_%s_%svs%s',...
                                subject, opts.epochWdW, all_classes{c2}, all_classes{c1});
                            
                            trials_selection = (all_Labels == c1-1) | (all_Labels == c2-1);
                            Labels_pair = all_Labels(trials_selection,:);
                            % Replace by [0,1] values for classification
                            Labels_pair(Labels_pair == c1-1) = 0;
                            Labels_pair(Labels_pair == c2-1) = 1;
                            
                            params.classes = all_classes([c1,c2]);
                            
                            if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                                params.CV = 'LOO';
                                %plot_params.error = true;
                            elseif ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo_tr')...
                                    && study_config.psd.chunks > 0
                                params.CV = 'LOOTR';
                                %params.ChunksDependency = 100.*TrialsInfo.Block(trials_selection)+...
                                %    TrialsInfo.Trial(trials_selection);
                                %plot_params.error = true;
                            else
                                if study_config.class.nbfolds < length(Labels_pair)
                                    params.CV = 'folds';
                                else
                                    params.CV = 'LOO';
                                    %plot_params.error = true;
                                end
                            end
                            
                            % Random prediction
                            RandomRes = randomPrediction(Labels_pair, NrepRandom);
                            RandomErr(:,c1,c2) = RandomRes.Err_bal;
                            
                            switch params.CV
                                case 'folds'
                                    results = load(sprintf('%s%s_results_%dfolds_%s', params.saveDataFolder, params.name, study_config.class.nbfolds, params.suffix));
                                    TestingErrAve(:,c1,c2) = 100*mean(results.CVstats.testing.Err_bal,2);
                                    TestingErrSE(:,c1,c2) = 100*std(results.CVstats.testing.Err_bal,[],2)/sqrt(size(results.CVstats.testing.Err_bal,2));
                                case 'LOO'
                                    results = load(sprintf('%s%s_results_LOO_%s', params.saveDataFolder, params.name, params.suffix));
                                    TestingErrAve(:,c1,c2) = 100*(1-(sum(squeeze(results.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,1,:)),[2,3])...
                                        +sum(squeeze(results.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
                                case 'LOOTR'
                                    results = load(sprintf('%s%s_results_LOOTR_%s', params.saveDataFolder, params.name, params.suffix));
                                    TestingErrAve(:,c1,c2) = 100*(1-(sum(squeeze(results.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,1,:)),[2,3])...
                                        +sum(squeeze(results.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
                                    
                            end
                            
                            %TestingErrAll(:,:,(c-1)*2+1) = 100*results_EC.CVstats.testing.Err;
                            %TestingErrAve(:,(c-1)*2+1) = 100*mean(results_EC.CVstats.testing.Err_bal,2);
                            TrainingErrAve(:,c1,c2) = 100*mean(results.CVstats.training.Err_bal,2);
                            %TestingErrSE(:,(c-1)*2+1) = 100*std(results_EC.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EC.CVstats.testing.Err_bal,2));
                            TrainingErrSE(:,c1,c2) = 100*std(results.CVstats.training.Err_bal,[],2)/sqrt(size(results.CVstats.training.Err_bal,2));
                            
                            [optTestErr, optInd] = min(TestingErrAve(:,c1,c2));
                            optNfeats(c1,c2) = study_config.class.feat2test(optInd);
                            optErrors(:,c1,c2) = [TrainingErrAve(optInd,c1,c2); optTestErr];
                            optSEs(:,c1,c2) = [TrainingErrSE(optInd,c1,c2), TestingErrSE(optInd,c1,c2)];
                            
                            sel = results.classifierInfo.Nb_feat == optNfeats(c1,c2);
                            all_feats_used = cell2mat(table2array(results.classifierInfo(sel,3)));
                            [chans, freqs] = ind2sub([n_chans, n_freqs], all_feats_used);
                            n_trials = size(all_feats_used,1);
                            for tr = 1:n_trials
                                for f = 1:optNfeats(c1,c2)
                                    feat_occurences(chans(tr,f),freqs(tr,f),c1,c2) = feat_occurences(chans(tr,f),freqs(tr,f),c1,c2)+1;
                                end
                            end
                            feat_occurences(:,:,c1,c2) = feat_occurences(:,:,c1,c2)./n_trials;
                        end
                    end
                    
                    %% Stats
                    err_distribution = RandomErr(:);
                    err_distribution(isnan(err_distribution)) = [];
                    chancelevel = mean(err_distribution);
                    signiflevel = quantile(err_distribution,0.05);
                    
                    % Computing significance according to
                    % Benjamini-Hochberg procedure for false discovery rate
                    pvalues_test = squeeze(sum(100*err_distribution < optErrors(2,:,:),1))/length(err_distribution);
                    pvalues_rank = nan(size(pvalues_test));
                    for c1 = 1:n_classes-1
                        for c2 = c1+1:n_classes
                            if c1 == 1 && c2 == 2
                                pvalues_rank(c1,c2) = 1;
                            else
                                seen_pvals = ~isnan(pvalues_rank);
                                sorted_seen_pvals = sort(pvalues_test(seen_pvals), 'ascend');
                                
                                if any(pvalues_test(c1,c2) == sorted_seen_pvals)
                                    new_rank = find(pvalues_test(c1,c2) == sorted_seen_pvals,1);
                                else
                                    new_rank = find(pvalues_test(c1,c2) < sorted_seen_pvals,1);
                                    if isempty(new_rank)
                                        % Largest p_val found so far
                                        new_rank = length(sorted_seen_pvals)+1;
                                    end
                                end
                                pvalues_rank(c1,c2) = new_rank;
                                ranks2update = pvalues_test > pvalues_test(c1,c2);
                                pvalues_rank(ranks2update) = pvalues_rank(ranks2update)+1;
                            end
                        end
                    end
                    
                    FDR_signif = pvalues_test < 0.05*pvalues_rank./sum(~isnan(pvalues_rank), 'all');
                    largest_rank_accepted = max(pvalues_rank(FDR_signif),[],'all');
                    if isempty(largest_rank_accepted)
                        signif_pvals = false(size(pvalues_rank));
                    else
                        signif_pvals = pvalues_rank <= largest_rank_accepted;
                    end
                    
                    %% Figure
                    %leg = [];
                    %leg_labels = {};
                    switch user
                        case 'JB'
                            figure('Position', [100 100 1200 700],'Visible','Off')
                        otherwise
                            figure
                    end
                    for c1 = 1:n_classes-1
                        for c2 = c1+1:n_classes
                            subplot(n_classes, n_classes, n_classes*(c1-1) + c2);
                            set(gca,'YAxisLocation','right')
                            
                            hold on
                            %yline(mean(err_distribution, 'omitnan'), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
                            yl1 = yline(100*(1-chancelevel), 'k--', 'LineWidth', 1.5);
                            %yline(100*quantile(err_distribution,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05');
                            yl2 = yline(100*(1-signiflevel), 'k:', 'LineWidth', 2);
                            
                            bar_plot = bar(100-optErrors(:,c1,c2), 0.8);
                            
                            ylim([0,100]);
                            
                            if signif_pvals(c1,c2)
                                ylims=get(gca,'ylim');
                                xlims=get(gca,'xlim');
                                text(xlims(1)+0.5*(xlims(2)-xlims(1)),ylims(1)+0.85*(ylims(2)-ylims(1)),...
                                    sprintf('*p=%.2g',pvalues_test(c1,c2)), 'Fontsize',10);
                                
                                %a = annotation(sp, 'textbox', [.45,.9,.1,.1], 'String', '*');
                                %a.LineStyle = 'none';
                            end
                            
                            xticks(1:2);
                            xticklabels({'TRAIN','TEST'});
                            
                            if c2 == n_classes
                                ylabel(all_classes{c1},'FontWeight', 'bold');
                            end
                            
                            if c1 == 1
                                title(all_classes{c2});
                            end
                            
                            subplot(n_classes, n_classes, n_classes*(c2-1) + c1);
                            hold on
                            imagesc(100*feat_occurences(:,:,c1,c2), [0,100])
                            if contains(study_config.feat.folderName, 'occChansExt2')
                                yline(16.5, 'w:', 'LineWidth', 1.5);
                                yline(22.5, 'w:', 'LineWidth', 1.5);
                                if c1 == 1
                                    yticks([8.5,19.5,30.5]);
                                    yticklabels({'L','Z','R'});
                                    ylabel(all_classes{c2},'FontWeight', 'bold');
                                else
                                    yticks([]);
                                end
                            elseif contains(study_config.feat.folderName, 'occChansExt')
                                yline(5.5, 'w:', 'LineWidth', 1.5);
                                yline(18.5, 'w:', 'LineWidth', 1.5);
                                if c1 == 1
                                    yticks([2.5,12,24.5]);
                                    yticklabels({'Z','L','R'});
                                    ylabel(all_classes{c2},'FontWeight', 'bold');
                                else
                                    yticks([]);
                                end
                            else
                            end
                            xlim([dim2_labels(1)-0.5,dim2_labels(end)+0.5])
                            
                            if c2 == n_classes
                                xlabel(all_classes{c1},'FontWeight', 'bold');
                            else
                                xticks([]);
                            end
                            
                            
                            
                            if signif_pvals(c1,c2)
                                title(sprintf('Nfeats=%d *',optNfeats(c1,c2)));
                            else
                                title(sprintf('Nfeats=%d',optNfeats(c1,c2)));
                            end
                            
                            if c1 == 1 && c2 == 2
                                cb = colorbar('Location', 'northoutside');
                                cb.Position = [0.125,0.825,0.1,0.02];
                                cb.Label.String = 'Feature occurence (%)';
                                cb.Label.FontSize = 10;
                            end
                            
                            if c1 == n_classes-1 && c2 == n_classes
                                lgd = legend([yl1,yl2], {sprintf('Chance level: %.1f%%',100*(1-chancelevel)),...
                                    sprintf('Significant level (p<0.05): %.1f%%',100*(1-signiflevel))},...
                                    'Location', 'none');
                                lgd.Position = [0.8,0.1,0.1,0.1];
                                lgd.FontSize = 10;
                            end
                        end
                    end
                    
                    switch params.CV
                        case 'folds'
                            suptitle(['Classification accuracy for pairwise comparisons (mean across ', study_config.class.nbfolds,' folds)'])
                            saveCurrentFig(params.saveFigFolder,...
                                sprintf('%s_%s_PairwiseClassAccSummary_%dfolds_%s', subject, opts.epochWdW, study_config.class.nbfolds, params.suffix), {'fig','png','svg'}, []);
                        case 'LOO'
                            suptitle('Classification accuracy for pairwise comparisons (Leave-one-out)')
                            saveCurrentFig(params.saveFigFolder,...
                                sprintf('%s_%s_PairwiseClassAccSummary_LOO_%s', subject, opts.epochWdW, params.suffix), {'fig','png','svg'}, []);
                        case 'LOOTR'
                            suptitle('Classification accuracy for pairwise comparisons (Leave-one-out trial-based)')
                            saveCurrentFig(params.saveFigFolder,...
                                sprintf('%s_%s_PairwiseClassAccSummary_LOOTR_%s', subject, opts.epochWdW, params.suffix), {'fig','png','svg'}, []);
                    end
                    
                    
                else
                    NrepRandom = 10000;
                    RandomErr = nan(NrepRandom, 1);
                    %TestingErrAll = nan(length(study_config.class.feat2test), study_config.class.nbfolds, 2*numel(conds));
                    TestingErrAve = nan(length(study_config.class.feat2test), 1);
                    TrainingErrAve = nan(length(study_config.class.feat2test), 1);
                    TestingErrSE = nan(length(study_config.class.feat2test), 1);
                    TrainingErrSE = nan(length(study_config.class.feat2test), 1);
                    
                    params.name = sprintf('%s_%s', subject, opts.epochWdW);
                    %% Create labels
                    switch study_config.class.contrast
                        case 'TrialType'
                            Labels = zeros(size(TrialsInfo,1),1);
                            Labels(strcmp(TrialsInfo.TrialType, 'Disc')) = 1;
                            Labels(strcmp(TrialsInfo.TrialType, '3Rings')) = 2;
                            Labels(strcmp(TrialsInfo.TrialType, '5Rings')) = 3;
                            Labels(strcmp(TrialsInfo.TrialType, '45Bar')) = 4;
                            Labels(strcmp(TrialsInfo.TrialType, '-45Bar')) = 5;
                            
                            params.classes = {'Control','Disc','3Rings','5Rings','45Bar','-45Bar'};
                            %plot_params.error = true;
                            %plot_params.errorBal = false;
                            
                        case 'AllvsControl'
                            Labels = ones(size(TrialsInfo,1),1);
                            Labels(strcmp(TrialsInfo.TrialType, 'Control')) = 0;
                            
                            params.classes = {'Control','Stimulus'};
                            %plot_params.error = false;
                            %plot_params.errorBal = true;
                            %plot_params.CM = true;
                            
                        case 'DiscvsControl'
                            Labels = -ones(size(TrialsInfo,1),1);
                            Labels(strcmp(TrialsInfo.TrialType, 'Control')) = 0;
                            Labels(strcmp(TrialsInfo.TrialType, 'Disc')) = 1;
                            
                            trials_selection = ~(Labels == -1);
                            Labels = Labels(trials_selection,:);
                            Features = Features(:,trials_selection);
                            
                            params.classes = {'Control','Disc'};
                            %plot_params.error = true;
                            %plot_params.errorBal = false;
                            %plot_params.CM = true;
                        otherwise
                            error('Not coded yet')
                    end
                    
                    if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                        params.CV = 'LOO';
                        %plot_params.error = true;
                    else
                        if study_config.class.nbfolds < length(Labels)
                            params.CV = 'folds';
                        else
                            params.CV = 'LOO';
                            %plot_params.error = true;
                        end
                    end
                    
                    error('todo');
                end
            otherwise
                error('Not coded yet')
        end
        
        
        
        
        
        %         conds = unique(TrialsInfo_EC.Condition);
        %
        %         NrepRandom = 10000;
        %         RandomErr = nan(NrepRandom, 2*numel(conds));
        %         %TestingErrAll = nan(length(study_config.class.feat2test), study_config.class.nbfolds, 2*numel(conds));
        %         TestingErrAve = nan(length(study_config.class.feat2test), 2*numel(conds));
        %         TrainingErrAve = nan(length(study_config.class.feat2test), 2*numel(conds));
        %         TestingErrSE = nan(length(study_config.class.feat2test), 2*numel(conds));
        %         TrainingErrSE = nan(length(study_config.class.feat2test), 2*numel(conds));
        %
        %         for c = 1:numel(conds)
        %             if contains(study_config.class.contrast, 'Visibility') && contains(conds{c}, 'OFF')
        %                 continue
        %             end
        %             fprintf('Classification in %s condition\n', conds{c});
        %             params.name = sprintf('%s-%s', subject, conds{c});
        %
        %             trials_sel_EC = strcmp(TrialsInfo_EC.Condition, conds{c});
        %             TrialsInfo_EC_cond = TrialsInfo_EC(trials_sel_EC,:);
        %             trials_sel_EO = strcmp(TrialsInfo_EO.Condition, conds{c});
        %             TrialsInfo_EO_cond = TrialsInfo_EO(trials_sel_EO,:);
        %
        %             switch study_config.class.contrast
        %                 case 'TrialType'
        %                     Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
        %                     Labels_binary_EC(strcmp(TrialsInfo_EC_cond.TrialType, 'WithObject')) = 1;
        %
        %                     Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
        %                     Labels_binary_EO(strcmp(TrialsInfo_EO_cond.TrialType, 'WithObject')) = 1;
        %
        %                     params.classes = {'WithoutObject','WithObject'};
        %                 case 'Answer'
        %                     Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
        %                     Labels_binary_EC(strcmp(TrialsInfo_EC_cond.Answer, 'Present')) = 1;
        %
        %                     Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
        %                     Labels_binary_EO(strcmp(TrialsInfo_EO_cond.Answer, 'Present')) = 1;
        %
        %                     params.classes = {'Absent','Present'};
        %                 case 'Visibility_binary'
        %                     Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
        %                     Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
        %
        %                     Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
        %                     Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
        %
        %                     params.classes = {'TumblerFilmed','TumblerMissed'};
        %                 case 'Visibility_binary2'
        %                     trials_sel_EC = TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin |...
        %                         strcmp(TrialsInfo_EC_cond.TrialType, 'WithoutObject');
        %                     TrialsInfo_EC_cond = TrialsInfo_EC_cond(trials_sel_EC,:);
        %
        %                     trials_sel_EO = TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin |...
        %                         strcmp(TrialsInfo_EO_cond.TrialType, 'WithoutObject');
        %                     TrialsInfo_EO_cond = TrialsInfo_EO_cond(trials_sel_EO,:);
        %
        %                     Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
        %                     Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
        %
        %                     Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
        %                     Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
        %
        %                     plot_params.classes = {'TumblerFilmed','NoTumbler'};
        %                 otherwise
        %                     error('Not coded yet')
        %             end
        %
        %             if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
        %                 params.CV = 'LOO';
        %             else
        %                 if study_config.class.nbfolds < length(Labels_binary_EC)
        %                     params.CV = 'folds';
        %                 else
        %                     params.CV = 'LOO';
        %                 end
        %             end
        %
        %             %%%%%%%%%%%%%%%%%%%%%% EYES CLOSED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             if abs(length(Labels_binary_EC)/2 - sum(Labels_binary_EC)) <= length(Labels_binary_EC)/3
        %                 disp('Eyes closed phase');
        %
        %                 % Random prediction
        %                 RandomRes = randomPrediction(Labels_binary_EC, NrepRandom);
        %                 RandomErr(:,(c-1)*2+1) = RandomRes.Err_bal;
        %
        %                 params.phase = 'EC';
        %                 switch params.CV
        %                     case 'folds'
        %                         results = load(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, study_config.class.nbfolds, params.suffix));
        %                         TestingErrAve(:,(c-1)*2+1) = 100*mean(results.CVstats.testing.Err_bal,2);
        %                         TestingErrSE(:,(c-1)*2+1) = 100*std(results.CVstats.testing.Err_bal,[],2)/sqrt(size(results.CVstats.testing.Err_bal,2));
        %                     case 'LOO'
        %                         results = load(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix));
        %                         TestingErrAve(:,(c-1)*2+1) = 100*(1-(sum(squeeze(results.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,1,:)),[2,3])...
        %                             +sum(squeeze(results.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
        %                 end
        %
        %                 %TestingErrAll(:,:,(c-1)*2+1) = 100*results_EC.CVstats.testing.Err;
        %                 %TestingErrAve(:,(c-1)*2+1) = 100*mean(results_EC.CVstats.testing.Err_bal,2);
        %                 TrainingErrAve(:,(c-1)*2+1) = 100*mean(results.CVstats.training.Err_bal,2);
        %                 %TestingErrSE(:,(c-1)*2+1) = 100*std(results_EC.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EC.CVstats.testing.Err_bal,2));
        %                 TrainingErrSE(:,(c-1)*2+1) = 100*std(results.CVstats.training.Err_bal,[],2)/sqrt(size(results.CVstats.training.Err_bal,2));
        %             else
        %                 disp('Too unbalanced binary labels, skipping eyes closed phase');
        %             end
        %
        %             %%%%%%%%%%%%%%%%%%%%%% EYES OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             if abs(length(Labels_binary_EO)/2 - sum(Labels_binary_EO)) <= length(Labels_binary_EO)/3
        %                 disp('Eyes open phase');
        %
        %                 % Random prediction
        %                 RandomRes_EO = randomPrediction(Labels_binary_EO, NrepRandom);
        %                 RandomErr(:,(c-1)*2+2) = RandomRes_EO.Err_bal;
        %
        %                 params.phase = 'EO';
        %                 switch params.CV
        %                     case 'folds'
        %                         results_EO = load(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, study_config.class.nbfolds, params.suffix));
        %                         TestingErrAve(:,(c-1)*2+2) = 100*mean(results_EO.CVstats.testing.Err_bal,2);
        %                         TestingErrSE(:,(c-1)*2+2) = 100*std(results_EO.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EO.CVstats.testing.Err_bal,2));
        %                     case 'LOO'
        %                         results_EO = load(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix));
        %                         TestingErrAve(:,(c-1)*2+2) = 100*(1-(sum(squeeze(results_EO.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results_EO.CVstats.testing.CM(:,:,1,:)),[2,3])...
        %                             +sum(squeeze(results_EO.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results_EO.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
        %
        %                 end
        %                 %TestingErrAll(:,:,(c-1)*2+2) = 100*results_EO.CVstats.testing.Err;
        %                 %TestingErrAve(:,(c-1)*2+2) = 100*mean(results_EO.CVstats.testing.Err_bal,2);
        %                 TrainingErrAve(:,(c-1)*2+2) = 100*mean(results_EO.CVstats.training.Err_bal,2);
        %                 %TestingErrSE(:,(c-1)*2+2) = 100*std(results_EO.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EO.CVstats.testing.Err_bal,2));
        %                 TrainingErrSE(:,(c-1)*2+2) = 100*std(results_EO.CVstats.training.Err_bal,[],2)/sqrt(size(results_EO.CVstats.training.Err_bal,2));
        %             else
        %                 disp('Too unbalanced binary labels, skipping eyes open phase');
        %             end
        %         end
        %
        %         otherwise
        %             error('Not coded yet')
        %     end
        
    else
        error('Event not suitable for this analysis')
    end
    
    %     if ~all(isnan(TestingErrAve),'all')
    %         Ncomps = size(TestingErrAve,2);
    %         leg = [];
    %         leg_labels = {};
    %         figure
    %         hold on
    %         for col = 1:Ncomps
    %             if all(isnan(TestingErrAve(:,col)))
    %                 continue
    %             else
    %                 switch params.CV
    %                     case 'folds'
    %                         pl = errorbar(study_config.class.feat2test, TestingErrAve(:,col), TestingErrSE(:,col), 'Linewidth',2);
    %                         %errorbar(study_config.class.feat2test, TrainingErrAve(:,col), TrainingErrSE(:,col), '--', 'Linewidth',2, 'Color', pl.Color);
    %                     case 'LOO'
    %                         pl = plot(study_config.class.feat2test, TestingErrAve(:,col), '-o', 'Linewidth',2);
    %                         %plot(study_config.class.feat2test, TrainingErrAve(:,col), '--', 'Linewidth',2, 'Color', pl.Color);
    %                 end
    %                 errorbar(study_config.class.feat2test, TrainingErrAve(:,col), TrainingErrSE(:,col), '--', 'Linewidth',2, 'Color', pl.Color);
    %                 leg = [leg, pl];
    %                 if mod(col,2)==1
    %                     label = sprintf('%s-EC',conds{ceil(col/2)});
    %                 else
    %                     label = sprintf('%s-EO',conds{ceil(col/2)});
    %                 end
    %                 leg_labels = [leg_labels,label];
    %             end
    %         end
    %
    %         %if strcmp(study_config.class.contrast, 'TrialType')
    %         %    err_distribution = RandomRes_EO.Err;
    %         %else
    %         %    err_distribution = [RandomRes_EC.Err, RandomRes_EO.Err];
    %         %end
    %         %         yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
    %         %     patch([study_config.class.feat2test(1),study_config.class.feat2test(end),study_config.class.feat2test(end),study_config.class.feat2test(1)],...
    %         %         100*[quantile(err_distribution,0.05),quantile(err_distribution,0.05),quantile(err_distribution,0.95),quantile(err_distribution,0.95)],...
    %         %         [0.5,0.5,0.5], 'FaceAlpha', 0.2);
    %
    %         yline(50, 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
    %         switch study_config.class.contrast
    %             case 'TrialType'
    %                 err_distribution = RandomErr(:);
    %                 %yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
    %                 yline(100*quantile(err_distribution,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05');
    %             case 'Answer'
    %                 for c = 1:numel(conds)
    %                     err_distribution = RandomErr(:,(c-1)*2+[1:2]);
    %                     err_distribution = err_distribution(:);
    %                     if all(isnan(err_distribution))
    %                         continue
    %                     end
    %                     if c == 2
    %                         horAlign = 'Left';
    %                     else
    %                         horAlign = 'Right';
    %                     end
    %                     %yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', sprintf('Chance level %s',conds{c}),...
    %                     %    'LabelHorizontalAlignment', horAlign);
    %                     yline(100*quantile(err_distribution,0.05), 'k:', 'LineWidth', 2, 'Label', sprintf('p=0.05 %s',conds{c}),...
    %                         'LabelHorizontalAlignment', horAlign);
    %                 end
    %             case 'Visibility_binary'
    %                 err_distribution_EC = RandomErr(:,3);
    %                 %             yline(100*mean(err_distribution_EC), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EC',...
    %                 %                 'LabelHorizontalAlignment', 'Right');
    %                 yline(100*quantile(err_distribution_EC,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EC',...
    %                     'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Right');
    %
    %                 err_distribution_EO = RandomErr(:,4);
    %                 %             yline(100*mean(err_distribution_EO), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EO',...
    %                 %                 'LabelHorizontalAlignment', 'Left');
    %                 yline(100*quantile(err_distribution_EO,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EO',...
    %                     'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Left');
    %
    %             case 'Visibility_binary2'
    %                 err_distribution_EC = RandomErr(:,3);
    %                 %             yline(100*mean(err_distribution_EC), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EC',...
    %                 %                 'LabelHorizontalAlignment', 'Right');
    %                 yline(100*quantile(err_distribution_EC,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EC',...
    %                     'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Right');
    %
    %                 err_distribution_EO = RandomErr(:,4);
    %                 %             yline(100*mean(err_distribution_EO), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EO',...
    %                 %                 'LabelHorizontalAlignment', 'Left');
    %                 yline(100*quantile(err_distribution_EO,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EO',...
    %                     'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Left');
    %         end
    %         ylim([0,100])
    %         xlabel('Number of features used')
    %         ylabel('Percentage of error')
    %         legend(leg, leg_labels, 'Location', 'best')
    %
    %         switch params.CV
    %             case 'folds'
    %                 title(['Classification errors (mean across ', study_config.class.nbfolds,' folds)'])
    %                 saveCurrentFig(params.saveFigFolder,...
    %                     sprintf('%s_ClassErrBalSummary_%dfolds_%s', subject, study_config.class.nbfolds, params.suffix), {'png'}, []);
    %             case 'LOO'
    %                 title('Classification errors (Leave-one-out)')
    %                 saveCurrentFig(params.saveFigFolder,...
    %                     sprintf('%s_ClassErrBalSummary_LOO_%s', subject, params.suffix), {'png'}, []);
    %         end
    %
    %         Ncomps = size(TestingErrAve,2);
    %         Errors = zeros(Ncomps,2);
    %         SEs = zeros(Ncomps,2);
    %         %             if fold < length(nb_folds)+1
    %         %                 All_samples = zeros(size(TestingErrAve,2),nb_folds(fold));
    %         %             end
    %         optNfeats = zeros(1,Ncomps);
    %         for col = 1:Ncomps
    %             [optTestErr, optNfeats(col)] = min(TestingErrAve(:,col));
    %             Errors(col,:) = [TrainingErrAve(optNfeats(col),col), optTestErr];
    %             SEs(col,:) = [TrainingErrSE(optNfeats(col),col), TestingErrSE(optNfeats(col),col)];
    %
    %             %                 if fold < length(nb_folds)+1
    %             %                     All_samples(col,:) = TestingErrAll(optNfeats(col),:,col);
    %             %                 end
    %
    %             %         if mod(col,2)==1
    %             %             [Feats, AppearanceRate, ~] = calcStabilityIndex(results_EO{ceil(col/2)}.classifierInfo,features2test);
    %             %         else
    %             %             [Feats, AppearanceRate, ~] = calcStabilityIndex(results_EC{ceil(col/2)}.classifierInfo,features2test);
    %             %         end
    %             %         figName = classes_names{col};
    %             %
    %             %         freq_vect = freq_sel(1):freq_sel(2);
    %             %         n_feats = length(Feats{optNfeats(col)});
    %             %         featNames = cell(n_feats,1);
    %             %         for f = 1:n_feats
    %             %             if mod(Feats{optNfeats(col)}(f),numel(chan_sel)) == 0
    %             %                 featNames{f} = [chan_sel{end},'-'];
    %             %             else
    %             %                 featNames{f} = [chan_sel{mod(Feats{optNfeats(col)}(f),numel(chan_sel))},'-'];
    %             %             end
    %             %
    %             %             featNames{f} = [featNames{f},...
    %             %                 num2str(freq_vect(ceil(Feats{optNfeats(col)}(f)/numel(chan_sel)))), 'Hz'];
    %             %         end
    %             %
    %             %         subplot(size(TestingErrAve,2)/2,2,col)
    %             %         bar(AppearanceRate{optNfeats(col)})
    %             %         xticks(1:n_feats)
    %             %         xticklabels(featNames)
    %             %         xtickangle(45)
    %             %         ylim([0,100])
    %             %         ylabel('Frequency of appearance across folds (%)')
    %             %         if optNfeats(col) == 1
    %             %             title([figName, ' (1 feature)'])
    %             %         else
    %             %             title([figName, ' (', num2str(optNfeats(col)), ' features)'])
    %             %         end
    %         end
    %
    %         %     if fold == length(nb_folds)+1
    %         %         suptitle('LOO - Features contribution to the optimal classification accuracy')
    %         %         saveCurrentFig([config.workingFolder, 'ClassifierResults', filesep,'NormPerFeat' filesep, save_folderName, filesep],...
    %         %             ['Features_appearance_LOO'], {'fig', 'svg', 'png'}, [1000,800])
    %         %     else
    %         %         suptitle([num2str(nb_folds(fold)),'folds - Features contribution to the optimal classification accuracy'])
    %         %         saveCurrentFig([config.workingFolder, 'ClassifierResults', filesep,'NormPerFeat' filesep, save_folderName, filesep],...
    %         %             ['Features_appearance_', num2str(nb_folds(fold)),'folds'], {'fig', 'svg', 'png'}, [1000,800])
    %         %     end
    %
    %         figure;
    %         hold on;
    %         %     yline(100*(1-mean(err_distribution)), 'k--', 'LineWidth', 1.5, 'Label', 'Mean chance lvl', 'LabelHorizontalAlignment', 'left');
    %         %     yline(100*(1-quantile(err_distribution,0.05)), 'k--', 'LineWidth', 1, 'Label', '95th quantile chance lvl', 'LabelHorizontalAlignment', 'left');
    %
    %
    %         pl_mean = plot([0,numel(conds)+1],[50,50], 'k--', 'LineWidth', 1.5);
    %         switch study_config.class.contrast
    %             case 'TrialType'
    %                 err_distribution = RandomErr(:);
    %                 pl_pval = plot([0,numel(conds)+1],100*(1-quantile(err_distribution,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
    %             case 'Answer'
    %                 err_distribution_OFF = RandomErr(:,1:2);
    %                 err_distribution_OFF = err_distribution_OFF(:);
    %                 pl_pvalOFF = plot([0,numel(conds)+1],100*(1-quantile(err_distribution_OFF,0.05))*ones(1,2), ':', 'LineWidth', 1.5);
    %
    %                 err_distribution_ON = RandomErr(:,3:4);
    %                 err_distribution_ON = err_distribution_ON(:);
    %                 pl_pvalON = plot([0,numel(conds)+1],100*(1-quantile(err_distribution_ON,0.05))*ones(1,2), ':', 'LineWidth', 1.5);
    %
    %             case 'Visibility_binary'
    %                 err_distribution_EC = RandomErr(:,3);
    %                 plot([0.6,1.4], 100*(1-quantile(err_distribution_EC,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
    %
    %                 err_distribution_EO = RandomErr(:,4);
    %                 pl_pval = plot([1.6,2.4], 100*(1-quantile(err_distribution_EO,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
    %
    %             case 'Visibility_binary2'
    %                 err_distribution_EC = RandomErr(:,3);
    %                 plot([0.6,1.4], 100*(1-quantile(err_distribution_EC,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
    %
    %                 err_distribution_EO = RandomErr(:,4);
    %                 pl_pval = plot([1.6,2.4], 100*(1-quantile(err_distribution_EO,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
    %         end
    %
    %         bar_plot = bar(100-reshape(Errors(:,2),[2,numel(conds)]), 0.8);
    %         xlim([0,numel(conds)+1])
    %         ylim([0,100])
    %         xticks(1:numel(conds))
    %         xticklabels({'Eyes closed','Eyes Open'})
    %         %xtickangle(45)
    %         if strcmp(study_config.class.contrast,'Answer')
    %             bar_plot(1).FaceColor = pl_pvalOFF.Color;
    %             bar_plot(2).FaceColor = pl_pvalON.Color;
    %             legend([bar_plot,pl_mean,pl_pvalOFF,pl_pvalON], [conds', 'Mean chance level','p=0.05 level for OFF','p=0.05 level for ON'],'Location','best')
    %         else
    %             legend([bar_plot,pl_mean,pl_pval], [conds', 'Mean chance level','p=0.05 level'],'Location','best')
    %         end
    %
    %         switch params.CV
    %             case 'folds'
    %                 title(['Best test classification accuracy (mean across ', study_config.class.nbfolds,' folds)'])
    %                 saveCurrentFig(params.saveFigFolder,...
    %                     sprintf('%s_BestClassAccBal_%dfolds_%s', subject, study_config.class.nbfolds, params.suffix), {'png'}, []);
    %             case 'LOO'
    %                 title('Best test classification accuracy (Leave-one-out)')
    %                 saveCurrentFig(params.saveFigFolder,...
    %                     sprintf('%s_BestClassAccBal_LOO_%s', subject, params.suffix), {'png'}, []);
    %         end
    %     end
end

switch lower(user)
    case 'jb'
        set(groot, 'DefaultFigureVisible', 'on');
end