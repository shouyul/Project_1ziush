%clear all;
%close all;
config_PIONEER_VEP;
study_config.epochs.limits_wdw = [0,10];
warning('Hardcoding epoch limit window, actually depends on subject (was 30, now 10)');

switch lower(user)
    case 'alex'
        % Correct bug happening when running this function after
        % PlotSummaryClassifierResults
        addpath('C:\Program Files\MATLAB\R2019a\toolbox\stats\stats\', '-end');
    case 'jb'
        % stop MATLAB from stealing focus
        set(groot, 'DefaultFigureVisible', 'off');
        % don't know if necessary in my case
        addpath('/Applications/MATLAB_R2021a.app/toolbox/stats/stats', '-end');
end

%% Single subject analysis
plot_params.model = study_config.class.model;
plot_params.AUC = false;
plot_params.CM = false;

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
    
    [plot_params.saveDataFolder, plot_params.saveFigFolder] = makeClassifierArchitecture(input_filepath, study_config);
    
    if strcmp(study_config.epochs.event, 'TrialStart')
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
        
        delims = strfind(N.epochedFile,'_');
        opts.epochWdW = N.epochedFile(delims(end-1)+1:delims(end)-1);
        fileName = makePSDFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo
        fileName = makePSDFileName(study_config.feat.type, subject, opts);
        load(fullfile(input_filepath_feats, fileName)); % Creates PSD_norm or Amps
        
        
        delims = strfind(fileName, '_');
        plot_params.suffix = fileName(delims(4)+1:end-4);
        
        switch study_config.feat.type
            case 'psd'
                switch study_config.psd.normModel
                    case 'gain'
                        % Convert to dB
                        PSD_dB = 10*log10(PSD_norm);
                        %% Convert PSD to Features
                        disp('Selecting Features');
                        [Features, chan_labels, dim2_labels] = select_features(PSD_dB, chanlocs,...
                            study_config.psd.FoI, study_config.feat);
                    case 'additive'
                        %% Convert PSD to Features
                        disp('Selecting Features');
                        [Features, chan_labels, dim2_labels] = select_features(PSD_norm, chanlocs,...
                            study_config.psd.FoI, study_config.feat);
                end
            case 'amp'
                %% Convert Amps to Features
                disp('Selecting Features');
                [Features, chan_labels, dim2_labels] = select_features(Amps, chanlocs,...
                    study_config.psd.FoI, study_config.feat);
        end
        
        
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
                    
                    plot_params.error = true;
                    plot_params.errorBal = false;
                    %plot_params.CM = true;
                    
                    for c1 = 1:numel(all_classes)-1
                        for c2 = c1+1:numel(all_classes)
                            plot_params.name = sprintf('%s_%s_%svs%s',...
                                subject, opts.epochWdW, all_classes{c2}, all_classes{c1});
                            
                            trials_selection = (all_Labels == c1-1) | (all_Labels == c2-1);
                            Labels_pair = all_Labels(trials_selection,:);
                            % Replace by [0,1] values for classification
                            Labels_pair(Labels_pair == c1-1) = 0;
                            Labels_pair(Labels_pair == c2-1) = 1;
                            Features_pair = Features(:,trials_selection);
                            
                            plot_params.classes = all_classes([c1,c2]);
                            
                            if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                                plot_params.CV = 'LOO';
                                %plot_params.error = true;
                            elseif ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo_tr')...
                                    && study_config.psd.chunks > 0
                                plot_params.CV = 'LOOTR';
                                plot_params.ChunksDependency = 100.*TrialsInfo.Block(trials_selection)+...
                                    TrialsInfo.Trial(trials_selection);
                                %plot_params.error = true;
                            else
                                if study_config.class.nbfolds < length(Labels_pair)
                                    plot_params.CV = 'folds';
                                else
                                    plot_params.CV = 'LOO';
                                    %plot_params.error = true;
                                end
                            end
                            
                            %% Classifier
                            %% Fisher score to get the most discriminative features
                            doFisherAnalysis(Features_pair', Labels_pair, chan_labels, dim2_labels,...
                                study_config.class.feat2test(end), plot_params);
                            %% Classification
                            CVstats = classifyWithFisher(Features_pair', Labels_pair, study_config.class, plot_params);
                            plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                        end
                    end
                else
                    plot_params.name = sprintf('%s_%s', subject, opts.epochWdW);
                    %% Create labels
                    switch study_config.class.contrast
                        case 'TrialType'
                            Labels = zeros(size(TrialsInfo,1),1);
                            Labels(strcmp(TrialsInfo.TrialType, 'Disc')) = 1;
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
                            
                            plot_params.error = true;
                            plot_params.errorBal = false;
                            
                        case 'AllvsControl'
                            Labels = ones(size(TrialsInfo,1),1);
                            Labels(strcmp(TrialsInfo.TrialType, 'Control')) = 0;
                            
                            plot_params.classes = {'Control','Stimulus'};
                            plot_params.error = false;
                            plot_params.errorBal = true;
                            plot_params.CM = true;
                            
                        case 'DiscvsControl'
                            Labels = -ones(size(TrialsInfo,1),1);
                            Labels(strcmp(TrialsInfo.TrialType, 'Control')) = 0;
                            Labels(strcmp(TrialsInfo.TrialType, 'Disc')) = 1;
                            
                            trials_selection = ~(Labels == -1);
                            Labels = Labels(trials_selection,:);
                            Features = Features(:,trials_selection);
                            
                            plot_params.classes = {'Control','Disc'};
                            plot_params.error = true;
                            plot_params.errorBal = false;
                            plot_params.CM = true;
                        otherwise
                            error('Not coded yet')
                    end
                    
                    if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                        plot_params.CV = 'LOO';
                        %plot_params.error = true;
                    elseif ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo_tr')...
                            && study_config.psd.chunks > 0
                        plot_params.CV = 'LOOTR';
                        plot_params.ChunksDependency = 100.*TrialsInfo.Block+TrialsInfo.Trial;
                        %plot_params.error = true;
                    else
                        if study_config.class.nbfolds < length(Labels)
                            plot_params.CV = 'folds';
                        else
                            plot_params.CV = 'LOO';
                            %plot_params.error = true;
                        end
                    end
                    
                    %% Classifier
                    %% Fisher score to get the most discriminative features
                    doFisherAnalysis(Features', Labels, chan_labels, dim2_labels,...
                        study_config.class.feat2test(end), plot_params);
                    %% Classification
                    CVstats = classifyWithFisher(Features', Labels, study_config.class, plot_params);
                    plotCVstats(CVstats, study_config.class.feat2test, plot_params);
                end
            otherwise
                error('Not coded yet')
        end
        
    else
        error('Event not suitable for this analysis')
    end
end

switch lower(user)
    case 'jb'
        set(groot, 'DefaultFigureVisible', 'on');
end