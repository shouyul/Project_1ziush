%clear all;
%close all;
config_PIONEER_Tumbler;

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
    
    if any(strcmp({'P1001', 'P1004'},subject)) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,20]; % in seconds
    elseif any(strcmp({'P1001old-1', 'P1001old-2','P1001-2', 'P1001-3', 'P1002-2', 'P1002-3','P1009','P1001-4','P1004-2'},subject)) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,15]; % in seconds
    end
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    input_filepath = N.searchFolder_3arch_rej_ICcats;
    
    if ~exist(fullfile(study_config.figures_folder,'Classifiers'), 'dir')
        mkdir(fullfile(study_config.figures_folder,'Classifiers'));
    end
    
    [params.saveDataFolder, params.saveFigFolder] = makeClassifierArchitecture(input_filepath, study_config);
    
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Load data
        %         EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath); % For chanlocs
        %         chanlocs = EEG.chanlocs(~strcmp({EEG.chanlocs.labels},study_config.eog_channels) &...
        %             ~contains({EEG.chanlocs.type}, 'MOCAP'));
        
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
        opts.phase = 'EC';
        
        fileName = makePSDFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EC
        
        opts.phase = 'EO';
        fileName = makePSDFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EO
        fileName = makePSDFileName(study_config.feat.type, subject, opts);
        
        
        delims = strfind(fileName, '_');
        params.suffix = fileName(delims(4)+1:end-4);
        
        %% Create labels
        switch study_config.class.condAnalysis
            case 'separate'
                conds = unique(TrialsInfo_EC.Condition);
                
                NrepRandom = 10000;
                RandomErr = nan(NrepRandom, 2*numel(conds));
                %TestingErrAll = nan(length(study_config.class.feat2test), study_config.class.nbfolds, 2*numel(conds));
                TestingErrAve = nan(length(study_config.class.feat2test), 2*numel(conds));
                TrainingErrAve = nan(length(study_config.class.feat2test), 2*numel(conds));
                TestingErrSE = nan(length(study_config.class.feat2test), 2*numel(conds));
                TrainingErrSE = nan(length(study_config.class.feat2test), 2*numel(conds));
                
                for c = 1:numel(conds)
                    if contains(study_config.class.contrast, 'Visibility') && contains(conds{c}, 'OFF')
                        continue
                    end
                    fprintf('Classification in %s condition\n', conds{c});
                    params.name = sprintf('%s-%s', subject, conds{c});
                    
                    trials_sel_EC = strcmp(TrialsInfo_EC.Condition, conds{c});
                    TrialsInfo_EC_cond = TrialsInfo_EC(trials_sel_EC,:);
                    trials_sel_EO = strcmp(TrialsInfo_EO.Condition, conds{c});
                    TrialsInfo_EO_cond = TrialsInfo_EO(trials_sel_EO,:);
                    
                    switch study_config.class.contrast
                        case 'TrialType'
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(strcmp(TrialsInfo_EC_cond.TrialType, 'WithObject')) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(strcmp(TrialsInfo_EO_cond.TrialType, 'WithObject')) = 1;
                            
                            params.classes = {'WithoutObject','WithObject'};
                        case 'Answer'
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(strcmp(TrialsInfo_EC_cond.Answer, 'Present')) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(strcmp(TrialsInfo_EO_cond.Answer, 'Present')) = 1;
                            
                            params.classes = {'Absent','Present'};
                        case 'Visibility_binary'
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            params.classes = {'TumblerFilmed','TumblerMissed'};
                        case 'Visibility_binary2'
                            trials_sel_EC = TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin |...
                                strcmp(TrialsInfo_EC_cond.TrialType, 'WithoutObject');
                            TrialsInfo_EC_cond = TrialsInfo_EC_cond(trials_sel_EC,:);
                            
                            trials_sel_EO = TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin |...
                                strcmp(TrialsInfo_EO_cond.TrialType, 'WithoutObject');
                            TrialsInfo_EO_cond = TrialsInfo_EO_cond(trials_sel_EO,:);
                            
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            plot_params.classes = {'TumblerFilmed','NoTumbler'};
                        otherwise
                            error('Not coded yet')
                    end
                    
                    if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                        params.CV = 'LOO';
                    else
                        if study_config.class.nbfolds < length(Labels_binary_EC)
                            params.CV = 'folds';
                        else
                            params.CV = 'LOO';
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES CLOSED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if abs(length(Labels_binary_EC)/2 - sum(Labels_binary_EC)) <= length(Labels_binary_EC)/3
                        disp('Eyes closed phase');
                        
                        % Random prediction
                        RandomRes_EC = randomPrediction(Labels_binary_EC, NrepRandom);
                        RandomErr(:,(c-1)*2+1) = RandomRes_EC.Err_bal;
                        
                        params.phase = 'EC';
                        switch params.CV
                            case 'folds'
                                results_EC = load(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, study_config.class.nbfolds, params.suffix));
                                TestingErrAve(:,(c-1)*2+1) = 100*mean(results_EC.CVstats.testing.Err_bal,2);
                                TestingErrSE(:,(c-1)*2+1) = 100*std(results_EC.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EC.CVstats.testing.Err_bal,2));
                            case 'LOO'
                                results_EC = load(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix));
                                TestingErrAve(:,(c-1)*2+1) = 100*(1-(sum(squeeze(results_EC.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results_EC.CVstats.testing.CM(:,:,1,:)),[2,3])...
                                    +sum(squeeze(results_EC.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results_EC.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
                        end
                        
                        %TestingErrAll(:,:,(c-1)*2+1) = 100*results_EC.CVstats.testing.Err;
                        %TestingErrAve(:,(c-1)*2+1) = 100*mean(results_EC.CVstats.testing.Err_bal,2);
                        TrainingErrAve(:,(c-1)*2+1) = 100*mean(results_EC.CVstats.training.Err_bal,2);
                        %TestingErrSE(:,(c-1)*2+1) = 100*std(results_EC.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EC.CVstats.testing.Err_bal,2));
                        TrainingErrSE(:,(c-1)*2+1) = 100*std(results_EC.CVstats.training.Err_bal,[],2)/sqrt(size(results_EC.CVstats.training.Err_bal,2));
                    else
                        disp('Too unbalanced binary labels, skipping eyes closed phase');
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if abs(length(Labels_binary_EO)/2 - sum(Labels_binary_EO)) <= length(Labels_binary_EO)/3
                        disp('Eyes open phase');
                        
                        % Random prediction
                        RandomRes_EO = randomPrediction(Labels_binary_EO, NrepRandom);
                        RandomErr(:,(c-1)*2+2) = RandomRes_EO.Err_bal;
                        
                        params.phase = 'EO';
                        switch params.CV
                            case 'folds'
                                results_EO = load(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, study_config.class.nbfolds, params.suffix));
                                TestingErrAve(:,(c-1)*2+2) = 100*mean(results_EO.CVstats.testing.Err_bal,2);
                                TestingErrSE(:,(c-1)*2+2) = 100*std(results_EO.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EO.CVstats.testing.Err_bal,2));
                            case 'LOO'
                                results_EO = load(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix));
                                TestingErrAve(:,(c-1)*2+2) = 100*(1-(sum(squeeze(results_EO.CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(results_EO.CVstats.testing.CM(:,:,1,:)),[2,3])...
                                    +sum(squeeze(results_EO.CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(results_EO.CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
                                
                        end
                        %TestingErrAll(:,:,(c-1)*2+2) = 100*results_EO.CVstats.testing.Err;
                        %TestingErrAve(:,(c-1)*2+2) = 100*mean(results_EO.CVstats.testing.Err_bal,2);
                        TrainingErrAve(:,(c-1)*2+2) = 100*mean(results_EO.CVstats.training.Err_bal,2);
                        %TestingErrSE(:,(c-1)*2+2) = 100*std(results_EO.CVstats.testing.Err_bal,[],2)/sqrt(size(results_EO.CVstats.testing.Err_bal,2));
                        TrainingErrSE(:,(c-1)*2+2) = 100*std(results_EO.CVstats.training.Err_bal,[],2)/sqrt(size(results_EO.CVstats.training.Err_bal,2));
                    else
                        disp('Too unbalanced binary labels, skipping eyes open phase');
                    end
                end
                
            otherwise
                error('Not coded yet')
        end
        
    else
        error('Event not suitable for this analysis')
    end
    
    if ~all(isnan(TestingErrAve),'all')
        Ncomps = size(TestingErrAve,2);
        leg = [];
        leg_labels = {};
        figure
        hold on
        for col = 1:Ncomps
            if all(isnan(TestingErrAve(:,col)))
                continue
            else
                switch params.CV
                    case 'folds'
                        pl = errorbar(study_config.class.feat2test, TestingErrAve(:,col), TestingErrSE(:,col), 'Linewidth',2);
                        %errorbar(study_config.class.feat2test, TrainingErrAve(:,col), TrainingErrSE(:,col), '--', 'Linewidth',2, 'Color', pl.Color);
                    case 'LOO'
                        pl = plot(study_config.class.feat2test, TestingErrAve(:,col), '-o', 'Linewidth',2);
                        %plot(study_config.class.feat2test, TrainingErrAve(:,col), '--', 'Linewidth',2, 'Color', pl.Color);
                end
                errorbar(study_config.class.feat2test, TrainingErrAve(:,col), TrainingErrSE(:,col), '--', 'Linewidth',2, 'Color', pl.Color);
                leg = [leg, pl];
                if mod(col,2)==1
                    label = sprintf('%s-EC',conds{ceil(col/2)});
                else
                    label = sprintf('%s-EO',conds{ceil(col/2)});
                end
                leg_labels = [leg_labels,label];
            end
        end
        
        %if strcmp(study_config.class.contrast, 'TrialType')
        %    err_distribution = RandomRes_EO.Err;
        %else
        %    err_distribution = [RandomRes_EC.Err, RandomRes_EO.Err];
        %end
        %         yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
        %     patch([study_config.class.feat2test(1),study_config.class.feat2test(end),study_config.class.feat2test(end),study_config.class.feat2test(1)],...
        %         100*[quantile(err_distribution,0.05),quantile(err_distribution,0.05),quantile(err_distribution,0.95),quantile(err_distribution,0.95)],...
        %         [0.5,0.5,0.5], 'FaceAlpha', 0.2);
        
        yline(50, 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
        switch study_config.class.contrast
            case 'TrialType'
                err_distribution = RandomErr(:);
                %yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
                yline(100*quantile(err_distribution,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05');
            case 'Answer'
                for c = 1:numel(conds)
                    err_distribution = RandomErr(:,(c-1)*2+[1:2]);                    
                    err_distribution = err_distribution(:);
                    if all(isnan(err_distribution))
                        continue
                    end
                    if c == 2
                        horAlign = 'Left';
                    else
                        horAlign = 'Right';
                    end
                    %yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', sprintf('Chance level %s',conds{c}),...
                    %    'LabelHorizontalAlignment', horAlign);
                    yline(100*quantile(err_distribution,0.05), 'k:', 'LineWidth', 2, 'Label', sprintf('p=0.05 %s',conds{c}),...
                        'LabelHorizontalAlignment', horAlign);
                end
            case 'Visibility_binary'
                err_distribution_EC = RandomErr(:,3);
                %             yline(100*mean(err_distribution_EC), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EC',...
                %                 'LabelHorizontalAlignment', 'Right');
                yline(100*quantile(err_distribution_EC,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EC',...
                    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Right');
                
                err_distribution_EO = RandomErr(:,4);
                %             yline(100*mean(err_distribution_EO), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EO',...
                %                 'LabelHorizontalAlignment', 'Left');
                yline(100*quantile(err_distribution_EO,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EO',...
                    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Left');
                
            case 'Visibility_binary2'
                err_distribution_EC = RandomErr(:,3);
                %             yline(100*mean(err_distribution_EC), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EC',...
                %                 'LabelHorizontalAlignment', 'Right');
                yline(100*quantile(err_distribution_EC,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EC',...
                    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Right');
                
                err_distribution_EO = RandomErr(:,4);
                %             yline(100*mean(err_distribution_EO), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level EO',...
                %                 'LabelHorizontalAlignment', 'Left');
                yline(100*quantile(err_distribution_EO,0.05), 'k:', 'LineWidth', 2, 'Label', 'p=0.05 EO',...
                    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment', 'Left');
        end
        ylim([0,100])
        xlabel('Number of features used')
        ylabel('Percentage of error')
        legend(leg, leg_labels, 'Location', 'best')
        
        switch params.CV
            case 'folds'
                title(['Classification errors (mean across ', study_config.class.nbfolds,' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBalSummary_%dfolds_%s', subject, study_config.class.nbfolds, params.suffix), {'png','svg'}, []);
            case 'LOO'
                title('Classification errors (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBalSummary_LOO_%s', subject, params.suffix), {'png','svg'}, []);
        end
        
        Ncomps = size(TestingErrAve,2);
        Errors = zeros(Ncomps,2);
        SEs = zeros(Ncomps,2);
        %             if fold < length(nb_folds)+1
        %                 All_samples = zeros(size(TestingErrAve,2),nb_folds(fold));
        %             end
        optNfeats = zeros(1,Ncomps);
        for col = 1:Ncomps
            [optTestErr, optNfeats(col)] = min(TestingErrAve(:,col));
            Errors(col,:) = [TrainingErrAve(optNfeats(col),col), optTestErr];
            SEs(col,:) = [TrainingErrSE(optNfeats(col),col), TestingErrSE(optNfeats(col),col)];
            
            %                 if fold < length(nb_folds)+1
            %                     All_samples(col,:) = TestingErrAll(optNfeats(col),:,col);
            %                 end
            
            %         if mod(col,2)==1
            %             [Feats, AppearanceRate, ~] = calcStabilityIndex(results_EO{ceil(col/2)}.classifierInfo,features2test);
            %         else
            %             [Feats, AppearanceRate, ~] = calcStabilityIndex(results_EC{ceil(col/2)}.classifierInfo,features2test);
            %         end
            %         figName = classes_names{col};
            %
            %         freq_vect = freq_sel(1):freq_sel(2);
            %         n_feats = length(Feats{optNfeats(col)});
            %         featNames = cell(n_feats,1);
            %         for f = 1:n_feats
            %             if mod(Feats{optNfeats(col)}(f),numel(chan_sel)) == 0
            %                 featNames{f} = [chan_sel{end},'-'];
            %             else
            %                 featNames{f} = [chan_sel{mod(Feats{optNfeats(col)}(f),numel(chan_sel))},'-'];
            %             end
            %
            %             featNames{f} = [featNames{f},...
            %                 num2str(freq_vect(ceil(Feats{optNfeats(col)}(f)/numel(chan_sel)))), 'Hz'];
            %         end
            %
            %         subplot(size(TestingErrAve,2)/2,2,col)
            %         bar(AppearanceRate{optNfeats(col)})
            %         xticks(1:n_feats)
            %         xticklabels(featNames)
            %         xtickangle(45)
            %         ylim([0,100])
            %         ylabel('Frequency of appearance across folds (%)')
            %         if optNfeats(col) == 1
            %             title([figName, ' (1 feature)'])
            %         else
            %             title([figName, ' (', num2str(optNfeats(col)), ' features)'])
            %         end
        end
        
        %     if fold == length(nb_folds)+1
        %         suptitle('LOO - Features contribution to the optimal classification accuracy')
        %         saveCurrentFig([config.workingFolder, 'ClassifierResults', filesep,'NormPerFeat' filesep, save_folderName, filesep],...
        %             ['Features_appearance_LOO'], {'fig', 'svg', 'png'}, [1000,800])
        %     else
        %         suptitle([num2str(nb_folds(fold)),'folds - Features contribution to the optimal classification accuracy'])
        %         saveCurrentFig([config.workingFolder, 'ClassifierResults', filesep,'NormPerFeat' filesep, save_folderName, filesep],...
        %             ['Features_appearance_', num2str(nb_folds(fold)),'folds'], {'fig', 'svg', 'png'}, [1000,800])
        %     end
        
        figure;
        hold on;
        %     yline(100*(1-mean(err_distribution)), 'k--', 'LineWidth', 1.5, 'Label', 'Mean chance lvl', 'LabelHorizontalAlignment', 'left');
        %     yline(100*(1-quantile(err_distribution,0.05)), 'k--', 'LineWidth', 1, 'Label', '95th quantile chance lvl', 'LabelHorizontalAlignment', 'left');
        
        
        pl_mean = plot([0,numel(conds)+1],[50,50], 'k--', 'LineWidth', 1.5);
        switch study_config.class.contrast
            case 'TrialType'
                err_distribution = RandomErr(:);
                pl_pval = plot([0,numel(conds)+1],100*(1-quantile(err_distribution,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
            case 'Answer'
                err_distribution_OFF = RandomErr(:,1:2);
                err_distribution_OFF = err_distribution_OFF(:);
                pl_pvalOFF = plot([0,numel(conds)+1],100*(1-quantile(err_distribution_OFF,0.05))*ones(1,2), ':', 'LineWidth', 1.5);
                
                err_distribution_ON = RandomErr(:,3:4);
                err_distribution_ON = err_distribution_ON(:);
                pl_pvalON = plot([0,numel(conds)+1],100*(1-quantile(err_distribution_ON,0.05))*ones(1,2), ':', 'LineWidth', 1.5);
                
            case 'Visibility_binary'
                err_distribution_EC = RandomErr(:,3);
                plot([0.6,1.4], 100*(1-quantile(err_distribution_EC,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
                
                err_distribution_EO = RandomErr(:,4);
                pl_pval = plot([1.6,2.4], 100*(1-quantile(err_distribution_EO,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
                
            case 'Visibility_binary2'
                err_distribution_EC = RandomErr(:,3);
                plot([0.6,1.4], 100*(1-quantile(err_distribution_EC,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
                
                err_distribution_EO = RandomErr(:,4);
                pl_pval = plot([1.6,2.4], 100*(1-quantile(err_distribution_EO,0.05))*ones(1,2), 'k:', 'LineWidth', 1.5);
        end
        
        bar_plot = bar(100-reshape(Errors(:,2),[2,numel(conds)]), 0.8);
        xlim([0,numel(conds)+1])
        ylim([0,100])
        xticks(1:numel(conds))
        xticklabels({'Eyes closed','Eyes Open'})
        %xtickangle(45)
        if strcmp(study_config.class.contrast,'Answer')
            bar_plot(1).FaceColor = pl_pvalOFF.Color;
            bar_plot(2).FaceColor = pl_pvalON.Color;
            legend([bar_plot,pl_mean,pl_pvalOFF,pl_pvalON], [conds', 'Mean chance level','p=0.05 level for OFF','p=0.05 level for ON'],'Location','best')
        else
            legend([bar_plot,pl_mean,pl_pval], [conds', 'Mean chance level','p=0.05 level'],'Location','best')
        end
        
       switch params.CV
            case 'folds'
                title(['Best test classification accuracy (mean across ', study_config.class.nbfolds,' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_BestClassAccBal_%dfolds_%s', subject, study_config.class.nbfolds, params.suffix), {'png','svg'}, []);
            case 'LOO'
                title('Best test classification accuracy (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_BestClassAccBal_LOO_%s', subject, params.suffix), {'png','svg'}, []);
        end
    end
end