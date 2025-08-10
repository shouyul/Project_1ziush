%clear all;
%close all;
%config_PIONEER_Tumbler;
switch lower(user)
    case 'alex'
        % Correct bug happening when running this function after
        % PlotSummaryClassifierResults
        addpath('C:\Program Files\MATLAB\R2019a\toolbox\stats\stats\', '-end');
end

%% Single subject analysis
plot_params.model = study_config.class.model;
plot_params.AUC = true;
plot_params.CM = true;

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
    elseif any(strcmp({'P1001old-1', 'P1001old-2','P1001-2', 'P1001-3', 'P1002-2', 'P1002-3', 'P1009','P1001-4','P1004-2'},subject)) &&...
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
    
    [plot_params.saveDataFolder, plot_params.saveFigFolder] = makeClassifierArchitecture(input_filepath, study_config);
    
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Load data
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath); % For chanlocs
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
        
        fileName = makePSDFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EC
        fileName = makePSDFileName(study_config.feat.type, subject, opts);
        load(fullfile(input_filepath_feats, fileName)); % Creates PSD_EC_norm or Amps_EC
        
        opts.phase = 'EO';
        fileName = makePSDFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EO
        fileName = makePSDFileName(study_config.feat.type, subject, opts);
        load(fullfile(input_filepath_feats, fileName)); % Creates PSD_EO_norm or Amps_EO
        
        % Define chan order for plots
        if any(strcmp({'P1002-2', 'P1002-3'},subject))
            new_chan_order = {'AF7','AF3','Fp1','Fpz','Fp2','AF4','AF8',...
                'F7','F5','F3','F1','Fz','F2','F4','F6','F8',...
                'FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8',...
                'M1','T7','C5','C3','C1','Cz','C2','C4','C6','T8','M2',...
                'TP7','CP5','CP3','CP1','CP2','CP4','CP6','TP8',...
                'P7','P5','P3','P1','Pz','P2','P4','P6','P8',...
                'PO7','PO5','PO3','POz','PO4','PO6','PO8',...
                'O1','Oz','O2'};
            chan_permut = 1:numel(new_chan_order);
            
            for ch = 1:numel(new_chan_order)
                chan_permut(ch) = find(strcmp({chanlocs.labels},new_chan_order{ch}));
            end
            
            chanlocs = chanlocs(chan_permut);
            PSD_EC_norm = PSD_EC_norm(chan_permut,:,:);
            PSD_EO_norm = PSD_EO_norm(chan_permut,:,:);
        end
        
        delims = strfind(fileName, '_');
        plot_params.suffix = fileName(delims(4)+1:end-4);
        
        switch study_config.feat.type
            case 'psd'
        switch study_config.psd.normModel
            case 'gain'
                % Convert to dB
                PSD_EC_dB = 10*log10(PSD_EC_norm);
                PSD_EO_dB = 10*log10(PSD_EO_norm);
                
                %% Convert PSD to Features
                disp('Selecting Features');
                [Features_EC, ~, dim2_labels_EC] = select_features(PSD_EC_dB, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
                [Features_EO, chan_labels, dim2_labels_EO] = select_features(PSD_EO_dB, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
                %             [Features_EO1, ~, ~] = select_features(PSD_EO1_dB, chanlocs,...
                %                 study_config.psd.FoI, study_config.feat);
                %             [Features_EO2, ~, ~] = select_features(PSD_EO2_dB, chanlocs,...
                %                 study_config.psd.FoI, study_config.feat);
                %             [Features_EO3, ~, ~] = select_features(PSD_EO3_dB, chanlocs,...
                %                 study_config.psd.FoI, study_config.feat);
                %             [Features_EO4, chan_labels, freq_labels] = select_features(PSD_EO4_dB, chanlocs,...
                %                 study_config.psd.FoI, study_config.feat);
            case 'additive'
                %% Convert PSD to Features
                disp('Selecting Features');
                [Features_EC, ~, dim2_labels_EC] = select_features(PSD_EC_norm, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
                [Features_EO, chan_labels, dim2_labels_EO] = select_features(PSD_EO_norm, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
        end
            case 'amp'
                %% Convert Amps to Features
                disp('Selecting Features');
                [Features_EC, ~, dim2_labels_EC] = select_features(Amps_EC, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
                [Features_EO, chan_labels, dim2_labels_EO] = select_features(Amps_EO, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
        end
        
        %% Create labels
        switch study_config.class.condAnalysis
            case 'separate'
                %                 if ~any(strcmp(study_config.psd.normTrialsGroup, 'perCondition'))
                %                     error('Normalization per condition required for that case');
                %                 end
                
                conds = unique(TrialsInfo_EC.Condition);
                for c = 1:numel(conds)
                    if contains(study_config.class.contrast, 'Visibility') && contains(conds{c}, 'OFF')
                        continue
                    end
                    fprintf('Classification in %s condition\n', conds{c});
                    plot_params.name = sprintf('%s-%s', subject, conds{c});
                    
                    trials_sel_EC = strcmp(TrialsInfo_EC.Condition, conds{c});
                    TrialsInfo_EC_cond = TrialsInfo_EC(trials_sel_EC,:);
                    Features_EC_cond = Features_EC(:,trials_sel_EC);
                    trials_sel_EO = strcmp(TrialsInfo_EO.Condition, conds{c});
                    TrialsInfo_EO_cond = TrialsInfo_EO(trials_sel_EO,:);
                    Features_EO_cond = Features_EO(:,trials_sel_EO);
                    %                     Features_EO1_cond = Features_EO1(:,trials_sel);
                    %                     Features_EO2_cond = Features_EO2(:,trials_sel);
                    %                     Features_EO3_cond = Features_EO3(:,trials_sel);
                    %                     Features_EO4_cond = Features_EO4(:,trials_sel);
                    
                    switch study_config.class.contrast
                        case 'TrialType'
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(strcmp(TrialsInfo_EC_cond.TrialType, 'WithObject')) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(strcmp(TrialsInfo_EO_cond.TrialType, 'WithObject')) = 1;
                            
                            plot_params.classes = {'WithoutObject','WithObject'};
                            plot_params.error = true;
                            plot_params.errorBal = false;
                        case 'Answer'
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(strcmp(TrialsInfo_EC_cond.Answer, 'Present')) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(strcmp(TrialsInfo_EO_cond.Answer, 'Present')) = 1;
                            
                            plot_params.classes = {'Absent','Present'};
                            plot_params.error = false;
                            plot_params.errorBal = true;
                        case 'Visibility_binary'
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            plot_params.classes = {'TumblerFilmed','TumblerMissed'};
                            plot_params.error = false;
                            plot_params.errorBal = true;
                        case 'Visibility_binary2'                            
                            trials_sel_EC = TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin |...
                                strcmp(TrialsInfo_EC_cond.TrialType, 'WithoutObject');
                            TrialsInfo_EC_cond = TrialsInfo_EC_cond(trials_sel_EC,:);
                            Features_EC_cond = Features_EC_cond(:,trials_sel_EC);
                            trials_sel_EO = TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin |...
                                strcmp(TrialsInfo_EO_cond.TrialType, 'WithoutObject');
                            TrialsInfo_EO_cond = TrialsInfo_EO_cond(trials_sel_EO,:);
                            Features_EO_cond = Features_EO_cond(:,trials_sel_EO);
                            
                            Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            plot_params.classes = {'TumblerFilmed','NoTumbler'};
                            plot_params.error = false;
                            plot_params.errorBal = true;
                        otherwise
                            error('Not coded yet')
                    end
                    
                    if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                        plot_params.CV = 'LOO';
                        %plot_params.error = true;
                    else
                        if study_config.class.nbfolds < length(Labels_binary_EC)
                            plot_params.CV = 'folds';
                        else
                            plot_params.CV = 'LOO';
                            %plot_params.error = true;
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES CLOSED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if abs(length(Labels_binary_EC)/2 - sum(Labels_binary_EC)) <= length(Labels_binary_EC)/3
                        disp('Eyes closed phase');
                        plot_params.phase = 'EC';
                        %% Fisher score to get the most discriminative features
                        doFisherAnalysis(Features_EC_cond', Labels_binary_EC, chan_labels, dim2_labels_EC,...
                            study_config.class.feat2test(end), plot_params);
                        %% Classification
                        CVstats = classifyWithFisher(Features_EC_cond', Labels_binary_EC, study_config.class, plot_params);
                        plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                    else
                        disp('Too unbalanced binary labels, skipping eyes closed phase');
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if abs(length(Labels_binary_EO)/2 - sum(Labels_binary_EO)) <= length(Labels_binary_EO)/3
                        disp('Eyes open phase');
                        plot_params.phase = 'EO';
                        %% Fisher score to get the most discriminative features
                        doFisherAnalysis(Features_EO_cond', Labels_binary_EO, chan_labels, dim2_labels_EO,...
                            study_config.class.feat2test(end), plot_params);
                        %% Classification
                        CVstats = classifyWithFisher(Features_EO_cond', Labels_binary_EO, study_config.class, plot_params);
                        plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                    else
                        disp('Too unbalanced binary labels, skipping eyes open phase');
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES OPEN 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                     disp('Eyes open 1 phase');
                    %                     plot_params.phase = 'EO1';
                    %                     %% Fisher score to get the most discriminative features
                    %                     doFisherAnalysis(Features_EO1_cond', Labels_binary, chan_labels, freq_labels,...
                    %                         study_config.class.feat2test(end), plot_params);
                    %                     %% Classification
                    %                     CVstats = classify(Features_EO1_cond', Labels_binary, study_config.class, plot_params);
                    %                     plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES OPEN 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                     disp('Eyes open 2 phase');
                    %                     plot_params.phase = 'EO2';
                    %                     %% Fisher score to get the most discriminative features
                    %                     doFisherAnalysis(Features_EO2_cond', Labels_binary, chan_labels, freq_labels,...
                    %                         study_config.class.feat2test(end), plot_params);
                    %                     %% Classification
                    %                     CVstats = classify(Features_EO2_cond', Labels_binary, study_config.class, plot_params);
                    %                     plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES OPEN 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                     disp('Eyes open 3 phase');
                    %                     plot_params.phase = 'EO3';
                    %                     %% Fisher score to get the most discriminative features
                    %                     doFisherAnalysis(Features_EO3_cond', Labels_binary, chan_labels, freq_labels,...
                    %                         study_config.class.feat2test(end), plot_params);
                    %                     %% Classification
                    %                     CVstats = classify(Features_EO3_cond', Labels_binary, study_config.class, plot_params);
                    %                     plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                    
                    %%%%%%%%%%%%%%%%%%%%%% EYES OPEN 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                     disp('Eyes open 4 phase');
                    %                     plot_params.phase = 'EO4';
                    %                     %% Fisher score to get the most discriminative features
                    %                     doFisherAnalysis(Features_EO4_cond', Labels_binary, chan_labels, freq_labels,...
                    %                         study_config.class.feat2test(end), plot_params);
                    %                     %% Classification
                    %                     CVstats = classify(Features_EO4_cond', Labels_binary, study_config.class, plot_params);
                    %                     plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                end
                
            otherwise
                error('Not coded yet')
        end
        
    else
        error('Event not suitable for this analysis')
    end
end