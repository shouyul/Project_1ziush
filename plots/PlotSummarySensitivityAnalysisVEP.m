clear all;
close all;
config_PIONEER_VEP;

%% Single subject analysis
params.model = study_config.class.model;

for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        switch user
            case 'sl'
                launchEEGLAB
            case 'JB'
                eeglab
        end
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
        %load(fullfile(input_filepath_feats, fileName)); % Creates PSD_norm or Amps
        
        delims = strfind(fileName, '_');
        params.suffix = fileName(delims(4)+1:end-4);
        
        %         [~, chan_labels, dim2_labels] = select_features(PSD_norm, chanlocs,...
        %             study_config.psd.FoI, study_config.feat);
        clear PSD_norm
        
        
        %% Create labels
        switch study_config.class.condAnalysis
            case ''
                if contains(study_config.class.contrast, 'Pairwise')
                    all_Labels = zeros(size(TrialsInfo,1),1);
                    all_Labels(strcmp(TrialsInfo.TrialType, 'Disc')) = 1;
                    all_Labels(strcmp(TrialsInfo.TrialType, '3Rings')) = 2;
                    all_Labels(strcmp(TrialsInfo.TrialType, '5Rings')) = 3;
                    all_Labels(strcmp(TrialsInfo.TrialType, '45Bar')) = 4;
                    all_Labels(strcmp(TrialsInfo.TrialType, '-45Bar')) = 5;
                    
                    all_classes = {'Control','Disc','3Rings','5Rings','45Bar','-45Bar'};
                    n_classes = numel(all_classes);
                    
                    NrepRandom = 1000;
                    n_pairs = nchoosek(n_classes,2);
                    
                    timelims = 0:5:30;
                    n_lims = length(timelims);
                    
                    TestingErrOpt = nan(n_classes, n_classes, n_lims, n_lims);
                    Pvalues_test = nan(n_classes, n_classes, n_lims, n_lims);
                    signif_pvals = nan(n_classes, n_classes, n_lims, n_lims);
                    
                    %                                         n_chans = numel(chan_labels);
                    %                     n_freqs = numel(dim2_labels);
                    %                     feat_occurences = zeros(n_chans, n_freqs, n_classes, n_classes);
                    for t1 = 1:n_lims
                        for t2 = 1:n_lims
                            if t1 < t2 %&& (t1 == 1 || t2 == 7 || t2-t1 == 1)
                                RandomErr = nan(NrepRandom, n_classes, n_classes);
                                
                                study_config.epochs.limits_wdw = timelims([t1,t2]);
                                N = makeFolderFileNames(study_config, subject);
                                delims = strfind(N.epochedFile,'_');
                                opts.epochWdW = N.epochedFile(delims(end-1)+1:delims(end)-1);
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
                                                TestingErrAve = 100*mean(results.CVstats.testing.Err_bal,2);
                                            case 'LOO'
                                                results = load(sprintf('%s%s_results_LOO_%s', params.saveDataFolder, params.name, params.suffix));
                                                TestingErrAve = 100*(1-(sum(squeeze(results.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,1,:)),[2,3])...
                                                    +sum(squeeze(results.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
                                            case 'LOOTR'
                                                results = load(sprintf('%s%s_results_LOOTR_%s', params.saveDataFolder, params.name, params.suffix));
                                                TestingErrAve = 100*(1-(sum(squeeze(results.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,1,:)),[2,3])...
                                                    +sum(squeeze(results.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
                                        end
                                        [TestingErrOpt(c1,c2,t1,t2), ~] = min(TestingErrAve);
                                    end
                                end
                                
                                
                                %% Stats
                                err_distribution = RandomErr(:);
                                err_distribution(isnan(err_distribution)) = [];
                                chancelevel = mean(err_distribution);
                                signiflevel = quantile(err_distribution,0.05);
                                
                                % Computing significance according to
                                % Benjamini-Hochberg procedure for false discovery rate
                                
                                pvalues_test = squeeze(sum(100*err_distribution < reshape(TestingErrOpt(:,:,t1,t2),[1,n_classes,n_classes]),1))/length(err_distribution);
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
                                    
                                    FDR_signif = pvalues_test < 0.05*pvalues_rank./sum(~isnan(pvalues_rank), 'all');
                                    largest_rank_accepted = max(pvalues_rank(FDR_signif),[],'all');
                                    if isempty(largest_rank_accepted)
                                        signif_pvals(:,:,t1,t2) = false(size(pvalues_rank));
                                    else
                                        signif_pvals(:,:,t1,t2) = pvalues_rank <= largest_rank_accepted;
                                    end
                                    Pvalues_test(:,:,t1,t2) = pvalues_test;
                                end
                            else
                                continue
                            end
                        end
                    end
                           
                    
                    %% Figure
                    %leg = [];
                    %leg_labels = {};
                    figure
                    for c1 = 1:n_classes-1
                        for c2 = c1+1:n_classes
                            %% Plot accuracy
                            acc_plot = subplot(n_classes, n_classes, n_classes*(c1-1) + c2);
                            imagesc(100-squeeze(TestingErrOpt(c1,c2,:,:)), [0,100]);
                            colormap(acc_plot, 'hot');
                            
                            set(gca,'XAxisLocation','top')
                            xticks(1:n_lims);
                            xticklabels(timelims);
                            
                            set(gca,'YAxisLocation','right')
                            yticks(1:n_lims);
                            yticklabels(timelims);
                            
                            if c1 == 1 && c2 == 2
                                cb = colorbar(acc_plot, 'Location', 'northoutside');
                                cb.Position = [0.125,0.865,0.1,0.02];
                                cb.Label.String = 'Test Class. Accuracy (%)';
                                cb.Label.FontSize = 10;
                                cb.TickDirection = 'out';
                                cb.TickLength = 0.02;
                            end
                            
                            for t1 = 1:n_lims
                                for t2 = 1:n_lims
                                    if t1 < t2 %&& (t1 == 1 || t2 == 7 || t2-t1 == 1)
                                        if signif_pvals(c1,c2,t1,t2)
                                            ylims=get(gca,'ylim');
                                            xlims=get(gca,'xlim');
                                            text(xlims(1)+((t2-0.6)/n_lims)*(xlims(2)-xlims(1)),ylims(1)+((t1-0.25)/n_lims)*(ylims(2)-ylims(1)),...
                                                '*', 'Fontsize',10);
                                        end
                                    end
                                end
                            end
                            
                            % Legends
                            if c1 == 1
                                xlabel('Stop (s)');
                                title(all_classes{c2});
                            end
                            
                            if c2 == n_classes
                                ylabel('Start (s)');
                                ylims=get(gca,'ylim');
                                xlims=get(gca,'xlim');
                                text(xlims(2)+0.25*(xlims(2)-xlims(1)),ylims(1)+0.5*(ylims(2)-ylims(1)),...
                                    all_classes{c1}, 'Fontsize', 9, 'FontWeight', 'bold');
                            end
                            
                            %% Plot p-values
                            pval_plot = subplot(n_classes, n_classes, n_classes*(c2-1) + c1);
                            pvals = squeeze(Pvalues_test(c1,c2,:,:));
                            pvals(isnan(pvals)) = 1;
                            imagesc(log10(pvals),[-4,-1]);
                            colormap(pval_plot,flipud(parula));
                            
                            %set(gca,'XAxisLocation','top')
                            xticks(1:n_lims);
                            xticklabels(timelims);
                            
                            %set(gca,'YAxisLocation','right')
                            yticks(1:n_lims);
                            yticklabels(timelims);
                            
                            if c1 == c2-1 && c2 == n_classes
                                cb = colorbar(pval_plot,'Ticks',[-4,-3,-2,-1.3010],...
                                    'TickLabels',{'0.0001','0.001','0.01','0.05'},...
                                    'Location', 'northoutside');
                                cb.Position = [0.8,0.15,0.1,0.02];
                                cb.Label.String = 'p-value';
                                cb.Label.FontSize = 10;
                                cb.TickDirection = 'out';
                                cb.TickLength = 0.02;
                            end
                            
                            % Legends
                            if c1 == 1
                                ylabel('Start (s)');
                                ylims=get(gca,'ylim');
                                xlims=get(gca,'xlim');
                                text(xlims(1)-0.2*(xlims(2)-xlims(1)),ylims(1)+0.5*(ylims(2)-ylims(1)),...
                                    all_classes{c2}, 'Fontsize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                            end
                            
                            if c2 == n_classes
                                xlabel('Stop (s)');
                                ylims=get(gca,'ylim');
                                xlims=get(gca,'xlim');
                                text(xlims(1)+0.5*(xlims(2)-xlims(1)),ylims(2)+0.5*(ylims(2)-ylims(1)),...
                                    all_classes{c1}, 'Fontsize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                            end
                        end
                    end
                    
                    switch params.CV
                        case 'folds'
                            suptitle(['Best test classification accuracy for pairwise comparisons (mean across ', study_config.class.nbfolds,' folds)'])
                            saveCurrentFig(params.saveFigFolder,...
                                sprintf('%s_WindowSensitivityAnalysis_%dfolds_%s', subject, study_config.class.nbfolds, params.suffix), {'fig'}, []);
                        case 'LOO'
                            suptitle('Best test classification accuracy for pairwise comparisons (Leave-one-out)')
                            saveCurrentFig(params.saveFigFolder,...
                                sprintf('%s_WindowSensitivityAnalysis_LOO_%s', subject, params.suffix), {'fig'}, []);
                        case 'LOOTR'
                            suptitle('Best test classification accuracy for pairwise comparisons (Leave-one-out trial-based)')
                            saveCurrentFig(params.saveFigFolder,...
                                sprintf('%s_WindowSensitivityAnalysis_LOOTR_%s', subject, params.suffix), {'fig'}, []);
                    end
                    
                else                    
                    error('todo');
                end
            otherwise
                error('Not coded yet')
        end
        
        
    else
        error('Event not suitable for this analysis')
    end
    
    
end