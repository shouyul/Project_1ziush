clear all;
%close all;
config_PIONEER_Tumbler;

%% Single subject analysis
plot_params.model = study_config.class.model;
plot_params.AUC = true;
plot_params.CM = true;
%MSname = 'Allsubjs';
MSname = 'P1001all';

if ~exist('EEG', 'var')
    switch user
        case 'Alex'
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

[plot_params.saveDataFolder, plot_params.saveFigFolder] = makeClassifierArchitecture(input_filepath, study_config);

if strcmp(study_config.epochs.event, 'EyesOpening')
    %% Load data
    EEG = pop_loadset('filename', N.epochedFile, 'filepath', N.searchFolder_3arch_rej_ICcats); % For chanlocs
    chanlocs = EEG.chanlocs(~strcmp({EEG.chanlocs.labels},study_config.eog_channels) &...
        ~contains({EEG.chanlocs.type}, 'MOCAP'));
    
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
    fileName = makePSDFileName(study_config.feat.type, MSname, opts);
    load(fullfile(input_filepath_feats, fileName)); % Creates PSD_EC_all_norm or Amps_EC_all
    
    opts.phase = 'EO';
    fileName = makePSDFileName('labels', MSname, opts);
    load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EO_all
    fileName = makePSDFileName(study_config.feat.type, MSname, opts);
    load(fullfile(input_filepath_feats, fileName)); % Creates PSD_EO_all_norm or Amps_EO_all
    
    delims = strfind(fileName, '_');
    plot_params.suffix = fileName(delims(4)+1:end-4);
    
    switch study_config.feat.type
        case 'psd'
            switch study_config.psd.normModel
                case 'gain'
                    % Convert to dB
                    PSD_EC_all_dB = 10*log10(PSD_EC_all_norm);
                    PSD_EO_all_dB = 10*log10(PSD_EO_all_norm);
                    
                    %% Convert PSD to Features
                    disp('Selecting Features');
                    [Features_EC_all, ~, dim2_labels_EC] = select_features(PSD_EC_all_dB, chanlocs,...
                        study_config.psd.FoI, study_config.feat);
                    [Features_EO_all, chan_labels, dim2_labels_EO] = select_features(PSD_EO_all_dB, chanlocs,...
                        study_config.psd.FoI, study_config.feat);
                case 'additive'
                    %% Convert PSD to Features
                    disp('Selecting Features');
                    [Features_EC_all, ~, dim2_labels_EC] = select_features(PSD_EC_all_norm, chanlocs,...
                        study_config.psd.FoI, study_config.feat);
                    [Features_EO_all, chan_labels, dim2_labels_EO] = select_features(PSD_EO_all_norm, chanlocs,...
                        study_config.psd.FoI, study_config.feat);
            end
        case 'amp'
            %% Convert Amps to Features
            disp('Selecting Features');
            [Features_EC_all, ~, dim2_labels_EC] = select_features(Amps_EC_all, chanlocs,...
                study_config.psd.FoI, study_config.feat);
            [Features_EO_all, chan_labels, dim2_labels_EO] = select_features(Amps_EO_all, chanlocs,...
                study_config.psd.FoI, study_config.feat);
    end
    
    %% Create labels
    switch study_config.class.condAnalysis
        case 'separate'
            %                 if ~any(strcmp(study_config.psd.normTrialsGroup, 'perCondition'))
            %                     error('Normalization per condition required for that case');
            %                 end
            
            conds = unique(TrialsInfo_EC_all.Condition);
            for c = 1:numel(conds)
                if contains(study_config.class.contrast, 'Visibility') && contains(conds{c}, 'OFF')
                    continue
                end
                fprintf('Classification in %s condition\n', conds{c});
                plot_params.name = sprintf('%s-%s', MSname, conds{c});
                
                trials_sel_EC_all = strcmp(TrialsInfo_EC_all.Condition, conds{c});
                TrialsInfo_EC_all_cond = TrialsInfo_EC_all(trials_sel_EC_all,:);
                Features_EC_all_cond = Features_EC_all(:,trials_sel_EC_all);
                trials_sel_EO_all = strcmp(TrialsInfo_EO_all.Condition, conds{c});
                TrialsInfo_EO_all_cond = TrialsInfo_EO_all(trials_sel_EO_all,:);
                Features_EO_all_cond = Features_EO_all(:,trials_sel_EO_all);
                
                switch study_config.class.contrast
                    case 'TrialType'
                        Labels_binary_EC_all = zeros(size(TrialsInfo_EC_all_cond,1),1);
                        Labels_binary_EC_all(strcmp(TrialsInfo_EC_all_cond.TrialType, 'WithObject')) = 1;
                        
                        Labels_binary_EO_all = zeros(size(TrialsInfo_EO_all_cond,1),1);
                        Labels_binary_EO_all(strcmp(TrialsInfo_EO_all_cond.TrialType, 'WithObject')) = 1;
                        
                        plot_params.classes = {'WithoutObject','WithObject'};
                        plot_params.error = true;
                        plot_params.errorBal = false;
                    case 'Answer'
                        Labels_binary_EC_all = zeros(size(TrialsInfo_EC_all_cond,1),1);
                        Labels_binary_EC_all(strcmp(TrialsInfo_EC_all_cond.Answer, 'Present')) = 1;
                        
                        Labels_binary_EO_all = zeros(size(TrialsInfo_EO_all_cond,1),1);
                        Labels_binary_EO_all(strcmp(TrialsInfo_EO_all_cond.Answer, 'Present')) = 1;
                        
                        plot_params.classes = {'Absent','Present'};
                        plot_params.error = false;
                        plot_params.errorBal = true;
                    case 'Visibility_binary'
                        Labels_binary_EC_all = zeros(size(TrialsInfo_EC_all_cond,1),1);
                        Labels_binary_EC_all(TrialsInfo_EC_all_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                        
                        Labels_binary_EO_all = zeros(size(TrialsInfo_EO_all_cond,1),1);
                        Labels_binary_EO_all(TrialsInfo_EO_all_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                        
                        plot_params.classes = {'TumblerFilmed','TumblerMissed'};
                        plot_params.error = false;
                        plot_params.errorBal = true;
                    otherwise
                        error('Not coded yet')
                end
                
                if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                    plot_params.CV = 'LOO';
                    plot_params.error = true;
                else
                    if study_config.class.nbfolds < length(Labels_binary_EC_all)
                        plot_params.CV = 'folds';
                    else
                        plot_params.CV = 'LOO';
                        plot_params.error = true;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%% EYES CLOSED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if abs(length(Labels_binary_EC_all)/2 - sum(Labels_binary_EC_all)) <= length(Labels_binary_EC_all)/3
                    disp('Eyes closed phase');
                    plot_params.phase = 'EC';
                    %% Fisher score to get the most discriminative features
                    doFisherAnalysis(Features_EC_all_cond', Labels_binary_EC_all, chan_labels, dim2_labels_EC,...
                        study_config.class.feat2test(end), plot_params);
                    %% Classification
                    CVstats = classifyWithFisher(Features_EC_all_cond', Labels_binary_EC_all, study_config.class, plot_params);
                    plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                else
                    disp('Too unbalanced binary labels, skipping eyes closed phase');
                end
                
                %%%%%%%%%%%%%%%%%%%%%% EYES OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if abs(length(Labels_binary_EO_all)/2 - sum(Labels_binary_EO_all)) <= length(Labels_binary_EO_all)/3
                    disp('Eyes open phase');
                    plot_params.phase = 'EO';
                    %% Fisher score to get the most discriminative features
                    doFisherAnalysis(Features_EO_all_cond', Labels_binary_EO_all, chan_labels, dim2_labels_EO,...
                        study_config.class.feat2test(end), plot_params);
                    %% Classification
                    CVstats = classifyWithFisher(Features_EO_all_cond', Labels_binary_EO_all, study_config.class, plot_params);
                    plotCVstats(CVstats, study_config.class.feat2test, plot_params);
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