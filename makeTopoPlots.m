clear all;
%close all;
config_PIONEER_Tumbler;

%% Single subject
plot_params.model = study_config.class.model;

for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        switch user
            case 'sl'
                launchEEGLAB
            case 'JB'
                eeglab
        end
    end
    
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    if any(strcmp({'P1001', 'P1004'},subject)) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,20]; % in seconds
    elseif any(strcmp({'P1001old-1', 'P1001old-2','P1001-2', 'P1001-3', 'P1002-2', 'P1002-3'},subject)) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,15]; % in seconds
    end
    
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    input_filepath = N.searchFolder_3arch_rej_ICcats;
    
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Load data
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath);
        
        % Select channels
        if ~isempty(study_config.eog_channels)
            EEG = pop_select(EEG, 'nochannel', study_config.eog_channels);
        end
        
        mocapChans = contains({EEG.chanlocs.type}, 'MOCAP');
        nb_chans = EEG.nbchan - sum(mocapChans);
        chanlocs = EEG.chanlocs(~mocapChans);
        
        nb_trials = EEG.trials;
        freqs = study_config.psd.FoI;
        nb_freqs = length(freqs);
        
        TrialsInfo = EEG.etc.epochInfo.Trials;
        TrialsInfo = TrialsInfo(:,1:8);
        
        [PSD_EC, TumbVis_EC, nb_chunks_EC] = computePSD_TumbVis(EEG, 'EC', study_config);
        [PSD_EO, TumbVis_EO, nb_chunks_EO] = computePSD_TumbVis(EEG, 'EO', study_config);
        
        TrialsInfo_EC = adaptTrialsInfo(TrialsInfo, TumbVis_EC, nb_chunks_EC);
        PSD_EC = reshape(permute(PSD_EC,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks_EC]);
        TrialsInfo_EO = adaptTrialsInfo(TrialsInfo, TumbVis_EO, nb_chunks_EO);
        PSD_EO = reshape(permute(PSD_EO,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks_EO]);
        
        %delims = strfind(fileName, '_');
        %plot_params.suffix = fileName(delims(4)+1:end-4);
        
        %         switch study_config.feat.type
        %             case 'psd'
        %                 switch study_config.psd.normModel
        %                     case 'gain'
        %                         % Convert to dB
        %                         PSD_EC_dB = 10*log10(PSD_EC);
        %                         PSD_EO_dB = 10*log10(PSD_EO);
        %
        %                         %                 %% Convert PSD to Features
        %                         %                 disp('Selecting Features');
        %                         %                 [Features_EC, ~, dim2_labels_EC] = select_features(PSD_EC_dB, chanlocs,...
        %                         %                     study_config.psd.FoI, study_config.feat);
        %                         %                 [Features_EO, chan_labels, dim2_labels_EO] = select_features(PSD_EO_dB, chanlocs,...
        %                         %                     study_config.psd.FoI, study_config.feat);
        %                     case 'additive'
        %                         %                 %% Convert PSD to Features
        %                         %                 disp('Selecting Features');
        %                         %                 [Features_EC, ~, dim2_labels_EC] = select_features(PSD_EC_norm, chanlocs,...
        %                         %                     study_config.psd.FoI, study_config.feat);
        %                         %                 [Features_EO, chan_labels, dim2_labels_EO] = select_features(PSD_EO_norm, chanlocs,...
        %                         %                     study_config.psd.FoI, study_config.feat);
        %                 end
        %             case 'amp'
        %                 %                 %% Convert Amps to Features
        %                 %                 disp('Selecting Features');
        %                 %                 [Features_EC, ~, dim2_labels_EC] = select_features(Amps_EC, chanlocs,...
        %                 %                     study_config.psd.FoI, study_config.feat);
        %                 %                 [Features_EO, chan_labels, dim2_labels_EO] = select_features(Amps_EO, chanlocs,...
        %                 %                     study_config.psd.FoI, study_config.feat);
        %         end
        
        %% Create labels
        switch study_config.class.condAnalysis
            case 'separate'
                conds = unique(TrialsInfo_EC.Condition);
                for c = 1:numel(conds)
                    if contains(study_config.class.contrast, 'Visibility') && contains(conds{c}, 'OFF')
                        continue
                    end
                    fprintf('Classification in %s condition\n', conds{c});
                    plot_params.name = sprintf('%s-%s', subject, conds{c});
                    
                    trials_sel_EC = strcmp(TrialsInfo_EC.Condition, conds{c});
                    %TrialsInfo_EC_cond = TrialsInfo_EC(trials_sel_EC,:);
                    %Features_EC_cond = Features_EC(:,trials_sel_EC);
                    trials_sel_EO = strcmp(TrialsInfo_EO.Condition, conds{c});
                    %TrialsInfo_EO_cond = TrialsInfo_EO(trials_sel_EO,:);
                    %Features_EO_cond = Features_EO(:,trials_sel_EO);
                    
                    switch study_config.class.contrast
                        case 'TrialType'
                            trials_cat_EC = strcmp(TrialsInfo_EC.TrialType, 'WithObject');
                            trials_cat_EO = strcmp(TrialsInfo_EO.TrialType, 'WithObject');
                            %                             Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            %                             Labels_binary_EC(strcmp(TrialsInfo_EC_cond.TrialType, 'WithObject')) = 1;
                            %
                            %                             Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            %                             Labels_binary_EO(strcmp(TrialsInfo_EO_cond.TrialType, 'WithObject')) = 1;
                            
                            plot_params.classes = {'WithoutObject','WithObject'};
                            %plot_params.error = true;
                            %plot_params.errorBal = false;
                        case 'Answer'
                            trials_cat_EC = strcmp(TrialsInfo_EC.Answer, 'Present');
                            trials_cat_EO = strcmp(TrialsInfo_EO.Answer, 'Present');
                            %                             Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            %                             Labels_binary_EC(strcmp(TrialsInfo_EC_cond.Answer, 'Present')) = 1;
                            %
                            %                             Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            %                             Labels_binary_EO(strcmp(TrialsInfo_EO_cond.Answer, 'Present')) = 1;
                            
                            plot_params.classes = {'Absent','Present'};
                            %plot_params.error = false;
                            %plot_params.errorBal = true;
                        case 'Visibility_binary'
                            trials_cat_EC = TrialsInfo_EC.TumbVis>=study_config.class.thresh_vis_bin;
                            trials_cat_EO = TrialsInfo_EO.TumbVis>=study_config.class.thresh_vis_bin;
                            %                             Labels_binary_EC = zeros(size(TrialsInfo_EC_cond,1),1);
                            %                             Labels_binary_EC(TrialsInfo_EC_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            %
                            %                             Labels_binary_EO = zeros(size(TrialsInfo_EO_cond,1),1);
                            %                             Labels_binary_EO(TrialsInfo_EO_cond.TumbVis>=study_config.class.thresh_vis_bin) = 1;
                            
                            plot_params.classes = {'TumblerFilmed','TumblerMissed'};
                            %plot_params.error = false;
                            %plot_params.errorBal = true;
                        otherwise
                            error('Not coded yet')
                    end
                    
                    %                     if ischar(study_config.class.nbfolds) && strcmp(study_config.class.nbfolds,'loo')
                    %                         plot_params.CV = 'LOO';
                    %                         plot_params.error = true;
                    %                     else
                    %                         if study_config.class.nbfolds < length(Labels_binary_EC)
                    %                             plot_params.CV = 'folds';
                    %                         else
                    %                             plot_params.CV = 'LOO';
                    %                             plot_params.error = true;
                    %                         end
                    %                     end
                    
                    %% Topoplots
                    if any(trials_sel_EO & trials_cat_EO)
                        %plot_params.freqs2plot = [8];% in Hz
                        %plot_params.chans2plot = {'Z19Z'}; % cell of channels
                        %plot_params.phase = 'EC';
                        %topoplotClassComparison(PSD_EC, trials_sel_EC, trials_cat_EC, chanlocs, plot_params);
                        
                        plot_params.freqs2plot = [5];% in Hz
                        plot_params.chans2plot = {'O1','Oz','O2','F1','Fz','F2'}; % cell of channels
                        plot_params.phase = 'EO';
                        topoplotClassComparison(PSD_EO, trials_sel_EO, trials_cat_EO, chanlocs, plot_params);
                        
                        plot_params.freqs2plot = [22];% in Hz
                        plot_params.chans2plot = {'O1','Oz','O2','F1','Fz','F2'}; % cell of channels
                        plot_params.phase = 'EO';
                        topoplotClassComparison(PSD_EO, trials_sel_EO, trials_cat_EO, chanlocs, plot_params);
                    end
                end
                
            otherwise
                error('Not coded yet')
        end
    end
end