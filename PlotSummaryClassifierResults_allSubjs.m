clear all;
%close all;
config_PIONEER_Tumbler;

%% All subjects analysis
params.model = study_config.class.model;
MSname = 'Allsubjs';
%MSname = 'P1001all';

if ~exist('EEG', 'var')
    switch user
        case 'sl'
            launchEEGLAB
        case 'JB'
            eeglab
    end
end

subject = study_config.subjects(1).id;

if strcmp(task,'Tumbler2020')
    study_config.epochs.limits_wdw = [-5,15]; % in seconds
end

%% Folders and Files names:
N = makeFolderFileNames(study_config, subject);
input_filepath = N.searchFolder_4arch_rej_ICcats;

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
    
    fileName = makePSDFileName('labels', MSname, opts);
    load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EC_all
    
    opts.phase = 'EO';
    fileName = makePSDFileName('labels', MSname, opts);
    load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EO_all
    fileName = makePSDFileName(study_config.feat.type, MSname, opts);    
    
    delims = strfind(fileName, '_');
    params.suffix = fileName(delims(4)+1:end-4);    
    
    %% Create labels
    switch study_config.class.condAnalysis
        case 'separate'
            conds = unique(TrialsInfo_EC_all.Condition);
            
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
                params.name = sprintf('%s-%s', MSname, conds{c});
                
                trials_sel_EC_all = strcmp(TrialsInfo_EC_all.Condition, conds{c});
                TrialsInfo_EC_all_cond = TrialsInfo_EC_all(trials_sel_EC_all,:);
                trials_sel_EO_all = strcmp(TrialsInfo_EO_all.Condition, conds{c});
                TrialsInfo_EO_all_cond = TrialsInfo_EO_all(trials_sel_EO_all,:);
                
                switch study_config.class.contrast
                    case 'TrialType'
                        Labels_binary_EC_all = zeros(size(TrialsInfo_EC_all_cond,1),1);
                        Labels_binary_EC_all(strcmp(TrialsInfo_EC_all_cond.TrialType, 'WithObject')) = 1;
                        
                        Labels_binary_EO_all = zeros(size(TrialsInfo_EO_all_cond,1),1);
                        Labels_binary_EO_all(strcmp(TrialsInfo_EO_all_cond.TrialType, 'WithObject')) = 1;
                        
                        params.classes = {'WithoutObject','WithObject'};
                    case 'Answer'
                        Labels_binary_EC_all = zeros(size(TrialsInfo_EC_all_cond,1),1);
                        Labels_binary_EC_all(strcmp(TrialsInfo_EC_all_cond.Answer, 'Present')) = 1;
                        
                        Labels_binary_EO_all = zeros(size(TrialsInfo_EO_all_cond,1),1);
                        Labels_binary_EO_all(strcmp(TrialsInfo_EO_all_cond.Answer, 'Present')) = 1;
                        
                        params.classes = {'Absent','Present'};
                    case 'Visibility_binary'
                        Labels_binary_EC_all = zeros(size(TrialsInfo_EC_all_cond,1),1);
                        Labels_binary_EC_all(TrialsInfo_EC_all_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                        
                        Labels_binary_EO_all = zeros(size(TrialsInfo_EO_all_cond,1),1);
                        Labels_binary_EO_all(TrialsInfo_EO_all_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                        
                        params.classes = {'TumblerFilmed','TumblerMissed'};
                    otherwise
                        error('Not coded yet')
                end
                
                if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                    params.CV = 'LOO';
                else
                    if study_config.class.nbfolds < length(Labels_binary_EC_all)
                        params.CV = 'folds';
                    else
                        params.CV = 'LOO';
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%% EYES CLOSED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if abs(length(Labels_binary_EC_all)/2 - sum(Labels_binary_EC_all)) <= length(Labels_binary_EC_all)/3
                    disp('Eyes closed phase');
                    
                    % Random prediction
                    RandomRes_EC_all = randomPrediction(Labels_binary_EC_all, 1000);
                    
                    params.phase = 'EC';
                    switch params.CV
                        case 'folds'
                            results_EC_all = load(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, study_config.class.nbfolds, params.suffix));
                        case 'LOO'
                            results_EC_all = load(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix));
                    end
                    
                    %TestingErrAll(:,:,(c-1)*2+1) = 100*results_EC.CVstats.testing.Err;
                    TestingErrAve(:,(c-1)*2+1) = 100*mean(results_EC_all.CVstats.testing.Err,2);
                    TrainingErrAve(:,(c-1)*2+1) = 100*mean(results_EC_all.CVstats.training.Err,2);
                    if strcmp(params.CV,'folds')
                        TestingErrSE(:,(c-1)*2+1) = 100*std(results_EC_all.CVstats.testing.Err,[],2)/sqrt(size(results_EC_all.CVstats.testing.Err,2));
                        TrainingErrSE(:,(c-1)*2+1) = 100*std(results_EC_all.CVstats.training.Err,[],2)/sqrt(size(results_EC_all.CVstats.training.Err,2));
                    end
                else
                    disp('Too unbalanced binary labels, skipping eyes closed phase');
                end
                
                %%%%%%%%%%%%%%%%%%%%%% EYES OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if abs(length(Labels_binary_EO_all)/2 - sum(Labels_binary_EO_all)) <= length(Labels_binary_EO_all)/3
                    disp('Eyes open phase');
                    
                    % Random prediction
                    RandomRes_EO_all = randomPrediction(Labels_binary_EO_all, 1000);
                    
                    params.phase = 'EO';
                    switch params.CV
                        case 'folds'
                            results_EO_all = load(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, study_config.class.nbfolds, params.suffix));
                        case 'LOO'
                            results_EO_all = load(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix));
                    end
                    %TestingErrAll(:,:,(c-1)*2+2) = 100*results_EO.CVstats.testing.Err;
                    TestingErrAve(:,(c-1)*2+2) = 100*mean(results_EO_all.CVstats.testing.Err,2);
                    TrainingErrAve(:,(c-1)*2+2) = 100*mean(results_EO_all.CVstats.training.Err,2);
                    if strcmp(params.CV,'folds')
                        TestingErrSE(:,(c-1)*2+2) = 100*std(results_EO_all.CVstats.testing.Err,[],2)/sqrt(size(results_EO_all.CVstats.testing.Err,2));
                        TrainingErrSE(:,(c-1)*2+2) = 100*std(results_EO_all.CVstats.training.Err,[],2)/sqrt(size(results_EO_all.CVstats.training.Err,2));
                    end
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
        %             switch col
        %                 case 1
        %                     leg_labels = [leg_labels,'GogglesOFF-EC'];
        %                 case 2
        %                     leg_labels = [leg_labels,'GogglesOFF-EO'];
        %                 case 3
        %                     leg_labels = [leg_labels,'GogglesON-EC'];
        %                 case 4
        %                     leg_labels = [leg_labels,'GogglesON-EO'];
        %             end
    end
end
if strcmp(study_config.class.contrast, 'TrialType')
    err_distribution = RandomRes_EO_all.Err;
else
    err_distribution = [RandomRes_EC_all.Err, RandomRes_EO_all.Err];
end
    yline(100*mean(err_distribution), 'k--', 'LineWidth', 1.5, 'Label', 'Chance level');
    patch([study_config.class.feat2test(1),study_config.class.feat2test(end),study_config.class.feat2test(end),study_config.class.feat2test(1)],...
        100*[quantile(err_distribution,0.05),quantile(err_distribution,0.05),quantile(err_distribution,0.95),quantile(err_distribution,0.95)],...
        [0.5,0.5,0.5], 'FaceAlpha', 0.2);
ylim([0,100])
xlabel('Number of features used')
ylabel('Percentage of error')
legend(leg, leg_labels)

switch params.CV
    case 'folds'
        title(['Classification errors (mean across ', study_config.class.nbfolds,' folds)'])
        saveCurrentFig(params.saveFigFolder,...
            sprintf('%s_ClassErrSummary_%dfolds_%s', MSname, study_config.class.nbfolds, params.suffix), {'png'}, []);
    case 'LOO'
        title('Classification errors (Leave-one-out)')
        saveCurrentFig(params.saveFigFolder,...
            sprintf('%s_ClassErrSummary_LOO_%s', MSname, params.suffix), {'png'}, []);
end


    Ncomps = size(TestingErrAve,2);
     Errors_all = zeros(Ncomps,2);
            SEs_all = zeros(Ncomps,2);
%             if fold < length(nb_folds)+1
%                 All_samples = zeros(size(TestingErrAve,2),nb_folds(fold));
%             end
            optNfeats = zeros(1,Ncomps);
            %figure
            for col = 1:Ncomps
                [optTestErr, optNfeats(col)] = min(TestingErrAve(:,col));
                Errors_all(col,:) = [TrainingErrAve(optNfeats(col),col), optTestErr];
                SEs_all(col,:) = [TrainingErrSE(optNfeats(col),col), TestingErrSE(optNfeats(col),col)];
%                 if fold < length(nb_folds)+1
%                     All_samples(col,:) = TestingErrAll(optNfeats(col),:,col);
%                 end

%                 if mod(col,2)==1
%                     [Feats, AppearanceRate, ~] = calcStabilityIndex(results_EO{ceil(col/2)}.classifierInfo,features2test);
%                 else
%                     [Feats, AppearanceRate, ~] = calcStabilityIndex(results_EC{ceil(col/2)}.classifierInfo,features2test);
%                 end
%                 figName = classes_names{col};
% 
%                 freq_vect = freq_sel(1):freq_sel(2);
%                 n_feats = length(Feats{optNfeats(col)});
%                 featNames = cell(n_feats,1);
%                 for f = 1:n_feats
%                     if mod(Feats{optNfeats(col)}(f),numel(chan_sel)) == 0
%                         featNames{f} = [chan_sel{end},'-'];
%                     else
%                         featNames{f} = [chan_sel{mod(Feats{optNfeats(col)}(f),numel(chan_sel))},'-'];
%                     end
% 
%                     featNames{f} = [featNames{f},...
%                         num2str(freq_vect(ceil(Feats{optNfeats(col)}(f)/numel(chan_sel)))), 'Hz'];
%                 end
% 
%                 subplot(size(TestingErrAve,2)/2,2,col)
%                 bar(AppearanceRate{optNfeats(col)})
%                 xticks(1:n_feats)
%                 xticklabels(featNames)
%                 xtickangle(45)
%                 ylim([0,100])
%                 ylabel('Frequency of appearance across folds (%)')
%                 if optNfeats(col) == 1
%                     title([figName, ' (1 feature)'])
%                 else
%                     title([figName, ' (', num2str(optNfeats(col)), ' features)'])
%                 end
            end

%             if fold == length(nb_folds)+1
%                 suptitle('LOO - Features contribution to the optimal classification accuracy')
%                 saveCurrentFig([config.workingFolder, 'ClassifierResults', filesep,'NormPerFeat' filesep, save_folderName, filesep],...
%                     ['Features_appearance_LOO'], {'fig', 'svg', 'png'}, [1000,800])
%             else
%                 suptitle([num2str(nb_folds(fold)),'folds - Features contribution to the optimal classification accuracy'])
%                 saveCurrentFig([config.workingFolder, 'ClassifierResults', filesep,'NormPerFeat' filesep, save_folderName, filesep],...
%                     ['Features_appearance_', num2str(nb_folds(fold)),'folds'], {'fig', 'svg', 'png'}, [1000,800])
%             end

    figure;
    hold on;
    yline(100*(1-mean(err_distribution)), 'k--', 'LineWidth', 1.5, 'Label', 'Mean chance lvl', 'LabelHorizontalAlignment', 'left');
    yline(100*(1-quantile(err_distribution,0.05)), 'k--', 'LineWidth', 1, 'Label', '95th quantile chance lvl', 'LabelHorizontalAlignment', 'left');
    %bar([0.75,1.25,2.75,3.25],100-reshape(Errors(:,2),[2,2]), 1)
    bar_plot = bar(100-reshape(Errors_all(:,2),[2,numel(conds)]), 0.8);
    xlim([-1,numel(conds)+1])
    ylim([0,100])
    xticks(1:numel(conds))
    xticklabels({'Eyes closed','Eyes Open'})
    %xtickangle(45)
    legend(bar_plot, conds)
    
    switch params.CV
        case 'folds'
            title(['Best test classification accuracy (mean across ', study_config.class.nbfolds,' folds)'])
            saveCurrentFig(params.saveFigFolder,...
                sprintf('%s_BestClassAcc_%dfolds_%s', subject, study_config.class.nbfolds, params.suffix), {'png'}, []);
        case 'LOO'
            title('Best test classification accuracy (Leave-one-out)')
            saveCurrentFig(params.saveFigFolder,...
                sprintf('%s_BestClassAcc_LOO_%s', subject, params.suffix), {'png'}, []);
    end    