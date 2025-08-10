clear all;
config_PIONEER_Tumbler;

% Choose MultiSubject name:
if isempty(study_config.channel_selection) && all(any(subject_inds' == [1,4],1))
    MSname = 'Allsubjs';
elseif ~isempty(study_config.channel_selection) && all(any(subject_inds' == [1:4],1))
    MSname = 'Allsubjs';
elseif length(subject_inds) == 3 && all(contains({study_config.subjects(subject_inds).id},'P1001'))
    MSname = 'P1001all';
else
    MSname = '';
end

MSname = ''; % Force not to compute MS data because of decrepancy in trial time duration

Amps_EC_all = [];
Amps_EO_all = [];
TrialsInfo_EC_all = [];
TrialsInfo_EO_all = [];
for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        switch user
            case 'Alex'
                launchEEGLAB
            case 'JB'
                eeglab
        end
    end
    
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    if any(strcmp({'P1001old-1', 'P1001old-2','P1001-2', 'P1001-3', 'P1002-2'},subject)) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,15]; % in seconds
    end
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    input_filepath = N.searchFolder_3arch_rej_ICcats;
    
    if strcmp(study_config.epochs.event, 'EyesOpening')
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath);
        
        % Select channels
        if ~isempty(study_config.eog_channels)
            EEG = pop_select(EEG, 'nochannel', study_config.eog_channels);
        end
        
        %         % Select time
        %         switch study_config.epochs.window
        %             case 'fixed'
        %                 % Eyes-closed period serves as a control
        %                 EEG_EC = pop_select(EEG, 'time', [study_config.epochs.limits_wdw(1),0]);
        %                 EEG_EO = pop_select(EEG, 'time', [0,study_config.epochs.limits_wdw(2)]);
        %             otherwise
        %                 error('Not coded yet')
        %         end
        
        mocapChans = contains({EEG.chanlocs.type}, 'MOCAP');
        nb_chans = EEG.nbchan - sum(mocapChans);
        nb_trials = EEG.trials;
        
        TrialsInfo = EEG.etc.epochInfo.Trials;
        TrialsInfo = TrialsInfo(:,1:8);
        
        [Amps_EC, TumbVis_EC] = computeAmps_TumbVis(EEG, 'EC', study_config);
        nb_times_EC = size(Amps_EC,2);
        nb_chunks_EC = size(Amps_EC,4);
        [Amps_EO, TumbVis_EO] = computeAmps_TumbVis(EEG, 'EO', study_config);
        nb_times_EO = size(Amps_EO,2);
        nb_chunks_EO = size(Amps_EO,4);
        
        TrialsInfo_EC = adaptTrialsInfo(TrialsInfo, TumbVis_EC, nb_chunks_EC);
        Amps_EC = reshape(permute(Amps_EC,[1,2,4,3]),[nb_chans, nb_times_EC, nb_trials*nb_chunks_EC]);
        TrialsInfo_EO = adaptTrialsInfo(TrialsInfo, TumbVis_EO, nb_chunks_EO);
        Amps_EO = reshape(permute(Amps_EO,[1,2,4,3]),[nb_chans, nb_times_EO, nb_trials*nb_chunks_EO]);
        
        % Concatenate over subjects:
        if ~strcmp(MSname,'')
            TrialsInfo_EC_all = cat(1, TrialsInfo_EC_all, TrialsInfo_EC);
            Amps_EC_all = cat(3,Amps_EC_all,Amps_EC);
            TrialsInfo_EO_all = cat(1, TrialsInfo_EO_all, TrialsInfo_EO);
            Amps_EO_all = cat(3,Amps_EO_all,Amps_EO);
        end
        
        %% Normalize data according to options
        %[Amps_EC_norm] = normalizePSD(Amps_EC, TrialsInfo_EC, study_config.psd);
        %[Amps_EO_norm] = normalizePSD(Amps_EO, TrialsInfo_EO, study_config.psd);
        
        %% Save data
        opts = study_config.amps;
        opts.SR = EEG.srate;
        opts.event = study_config.epochs.event;
        opts.phase = 'EC';
        
        fileName = makePSDFileName('labels', subject, opts);
        save(fullfile(input_filepath, fileName), sprintf('TrialsInfo_%s',opts.phase));
        
        fileName = makePSDFileName('amp', subject, opts);
        save(fullfile(input_filepath, fileName), sprintf('Amps_%s',opts.phase), '-v7.3');
        
        opts.phase = 'EO';
        fileName = makePSDFileName('labels', subject, opts);
        save(fullfile(input_filepath, fileName), sprintf('TrialsInfo_%s',opts.phase));
        
        fileName = makePSDFileName('amp', subject, opts);
        save(fullfile(input_filepath, fileName), sprintf('Amps_%s',opts.phase), '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end

%% Multi-subjects data
if ~strcmp(MSname,'')
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Normalize data according to options
        %[Amps_EC_all_norm] = normalizePSD(Amps_EC_all, TrialsInfo_EC_all, study_config.psd);
        %[Amps_EO_all_norm] = normalizePSD(Amps_EO_all, TrialsInfo_EO_all, study_config.psd);
        
        %% Save data
        opts.phase = 'EC';
        fileName = makePSDFileName('labels', 'AllSubjs', opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('TrialsInfo_%s_all',opts.phase));
        
        fileName = makePSDFileName('amp', 'AllSubjs', opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('Amps_%s_all', opts.phase), '-v7.3');
        
        opts.phase = 'EO';
        fileName = makePSDFileName('labels', 'AllSubjs', opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('TrialsInfo_%s_all',opts.phase));
        
        fileName = makePSDFileName('amp', 'AllSubjs', opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('Amps_%s_all', opts.phase), '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end


