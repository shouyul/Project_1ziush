clear all;
%close all;
config_PIONEER_Tumbler;

%% Single subject
plot_params.model = study_config.class.model;

for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        switch user
            case 'SL'
                launchEEGLAB
            case 'JB'
                eeglab
        end
    end
    
    subject_ind = 1;
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
    
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Load data
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath); % For chanlocs
        chanlocs = EEG.chanlocs(~strcmp({EEG.chanlocs.labels},study_config.eog_channels) &...
            ~contains({EEG.chanlocs.type}, 'MOCAP'));
        clear EEG
        
        opts = study_config.psd;
        opts.event = study_config.epochs.event;
        opts.phase = 'EC';
        
        opts.chunks = 0;        
        fileName = makeFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EC
        TrialsInfo_EC_noChunk = TrialsInfo_EC;
        fileName = makeFileName('psd', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates PSD_EC_norm
        PSD_EC_norm_noChunk = PSD_EC_norm;
        
                opts.chunks = 2.5;        
        fileName = makeFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EC
        TrialsInfo_EC_2_5Chunks = TrialsInfo_EC;
        fileName = makeFileName('psd', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates PSD_EC_norm
        PSD_EC_norm_2_5Chunks = PSD_EC_norm;
        
                        opts.chunks = 1.25;        
        fileName = makeFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EC
        TrialsInfo_EC_1_25Chunks = TrialsInfo_EC;
        fileName = makeFileName('psd', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates PSD_EC_norm
        PSD_EC_norm_1_25Chunks = PSD_EC_norm;
        
        for ch = 1:size(PSD_EC_norm_noChunk,1)
            for tr = 1:size(PSD_EC_norm_noChunk,3)
                figure;
                hold on;
                plot(1:40, 10*log10(PSD_EC_norm_noChunk(ch,:,tr)))
                plot(1:40, 10*log10(PSD_EC_norm_2_5Chunks(ch,:,1+2*(tr-1))), '--')
                plot(1:40, 10*log10(PSD_EC_norm_2_5Chunks(ch,:,2*tr)), '--')
                plot(1:40, 10*log10(PSD_EC_norm_1_25Chunks(ch,:,1+4*(tr-1))), ':')
                plot(1:40, 10*log10(PSD_EC_norm_1_25Chunks(ch,:,2+4*(tr-1))), ':')
                plot(1:40, 10*log10(PSD_EC_norm_1_25Chunks(ch,:,3+4*(tr-1))), ':')
                plot(1:40, 10*log10(PSD_EC_norm_1_25Chunks(ch,:,4*tr)), ':')
            end
        end
        
                for fr = 1:size(PSD_EC_norm_noChunk,2)
            for tr = 1:size(PSD_EC_norm_noChunk,3)
                figure;
                hold on;
                plot(10*log10(PSD_EC_norm_noChunk(:,fr,tr)))
                plot(10*log10(PSD_EC_norm_2_5Chunks(:,fr,1+2*(tr-1))), '--')
                plot(10*log10(PSD_EC_norm_2_5Chunks(:,fr,2*tr)), '--')
                plot(10*log10(PSD_EC_norm_1_25Chunks(:,fr,1+4*(tr-1))), ':')
                plot(10*log10(PSD_EC_norm_1_25Chunks(:,fr,2+4*(tr-1))), ':')
                plot(10*log10(PSD_EC_norm_1_25Chunks(:,fr,3+4*(tr-1))), ':')
                plot(10*log10(PSD_EC_norm_1_25Chunks(:,fr,4*tr)), ':')
            end
        end
        
        
        opts.phase = 'EO';
        fileName = makeFileName('labels', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates TrialsInfo_EO
        fileName = makeFileName('psd', subject, opts);
        load(fullfile(input_filepath, fileName)); % Creates PSD_EO_norm
        
        %                 delims = strfind(fileName, '_');
        %         plot_params.suffix = fileName(delims(4)+1:end-4);
        
        if strcmp(study_config.psd.normModel, 'gain')
            % Convert to dB
            PSD_EC_dB = 10*log10(PSD_EC_norm);
            PSD_EO_dB = 10*log10(PSD_EO_norm);
            
            %             %% Convert PSD to Features
            %             disp('Selecting Features');
            %             [Features_EC, ~, ~] = select_features(PSD_EC_dB, chanlocs,...
            %                 study_config.psd.FoI, study_config.feat);
            %             [Features_EO, chan_labels, freq_labels] = select_features(PSD_EO_dB, chanlocs,...
            %                 study_config.psd.FoI, study_config.feat);
        end
        
        %% Create labels
        switch study_config.class.condAnalysis
            case 'separate'
                if ~any(strcmp(study_config.psd.normTrialsGroup, 'perCondition'))
                    error('Normalization per condition required for that case');
                end
                
                conds = unique(TrialsInfo_EC.Condition);
                for c = 1:numel(conds)
                    if contains(study_config.class.contrast, 'Visibility') && contains(conds{c}, 'OFF')
                        continue
                    end
%                     if contains(conds{c}, 'OFF')
%                         continue
%                     end
                    
                    fprintf('%s condition\n', conds{c});
                    plot_params.name = sprintf('%s-%s', subject, conds{c});
                    
                    trials_sel_EC = strcmp(TrialsInfo_EC.Condition, conds{c});
                    %                     Features_EC_cond = Features_EC(:,trials_sel_EC);
                    trials_sel_EO = strcmp(TrialsInfo_EO.Condition, conds{c});
                    %                     Features_EO_cond = Features_EO(:,trials_sel_EO);
                    
                    switch study_config.class.contrast
                        case 'TrialType'
                            trials_cat_EC = strcmp(TrialsInfo_EC.TrialType, 'WithObject');
                            trials_cat_EO = strcmp(TrialsInfo_EO.TrialType, 'WithObject');
                            
                            plot_params.classes = {'WithoutObject','WithObject'};
                        case 'Answer'
                            trials_cat_EC = strcmp(TrialsInfo_EC.Answer, 'Present');
                            trials_cat_EO = strcmp(TrialsInfo_EO.Answer, 'Present');
                            
                            plot_params.classes = {'Absent','Present'};
                        case 'Visibility_binary'
                            trials_cat_EC = TrialsInfo_EC.TumbVis>1;
                            trials_cat_EO = TrialsInfo_EO.TumbVis>1;
                            
                            plot_params.classes = {'TumblerFilmed','TumblerMissed'};
                        otherwise
                            error('Not coded yet')
                    end
                    
                    %% Topoplots
                    plot_params.freqs2plot = [16];% in Hz
                    plot_params.phase = 'EC';
                    topoplotClassComparison(PSD_EC_norm, trials_sel_EC, trials_cat_EC, chanlocs, plot_params);
                    
                    %plot_params.freqs2plot = [7];% in Hz
                    plot_params.phase = 'EO';
                    topoplotClassComparison(PSD_EO_norm, trials_sel_EO, trials_cat_EO, chanlocs, plot_params);
                end
        end         
    end
end