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

PSD_EC_all = [];
PSD_EO_all = [];
TrialsInfo_EC_all = [];
TrialsInfo_EO_all = [];
for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    if any(strcmp({'P1001', 'P1004'},subject)) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,20]; % in seconds
    elseif any(strcmp({'P1001old-1', 'P1001old-2','P1001-2', 'P1001-3', 'P1002-2', 'P1002-3','P1009','P1001-4','P1004-2'},...
            subject)) &&...
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
        %         if any(mocapChans)
        %             nb_chans = EEG.nbchan - sum(mocapChans);
        %             EEG = pop_select(EEG, 'nochannel', mocapChans);
        %         else
        %             nb_chans = EEG.nbchan;
        %         end
        
        nb_trials = EEG.trials;
        freqs = study_config.psd.FoI;
        nb_freqs = length(freqs);
        %SR = EEG.srate;
        
        TrialsInfo = EEG.etc.epochInfo.Trials;
        TrialsInfo = TrialsInfo(:,1:8);
        
        [PSD_EC, TumbVis_EC, nb_chunks_EC] = computePSD_TumbVis(EEG, 'EC', study_config);
        [PSD_EO, TumbVis_EO, nb_chunks_EO] = computePSD_TumbVis(EEG, 'EO', study_config);
        
        %         if study_config.psd.chunks == 0
        %             nb_chunks_EC = 1;
        %             nb_chunks_EO = 1;
        %         else
        %             nb_chunks_EC = abs(study_config.epochs.limits_wdw(1))/study_config.psd.chunks;
        %             nb_chunks_EO = study_config.epochs.limits_wdw(2)/study_config.psd.chunks;
        %         end
        
        %         PSD_EC = zeros(nb_chans, nb_freqs, nb_trials, nb_chunks_EC);
        %         PSD_EO = zeros(nb_chans, nb_freqs, nb_trials, nb_chunks_EO);
        %         if any(mocapChans)
        %             TumbVisChan = strcmp({EEG.chanlocs.labels}, 'TumblerVisibility');
        %             TumbVis_EC = zeros(nb_trials, nb_chunks_EC);
        %             TumbVis_EO = zeros(nb_trials, nb_chunks_EO);
        %         end
        
        %         switch study_config.psd.method
        %             case 'pwelch'
        %                 % Compute PSD with pwelch method
        %                 overlap = SR/25;
        %                 if nb_chunks_EC == 1
        %                     for tr = 1:nb_trials
        %                         PSD_EC(:,:,tr,1) = pwelch(squeeze(EEG_EC.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %
        %                         if any(mocapChans)
        %                             TumbVis_EC(tr,1) = mean(squeeze(EEG_EC.data(TumbVisChan,:,tr)));
        %                         end
        %                     end
        %                 else
        %                     for ch = 1:nb_chunks_EC
        %                         EEG_EC_sel = pop_select(EEG, 'time', study_config.epochs.limits_wdw(1)+...
        %                             [ch-1,ch].*study_config.psd.chunks);
        %                         for tr = 1:nb_trials
        %                             PSD_EC(:,:,tr,ch) = pwelch(squeeze(EEG_EC_sel.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %                             if any(mocapChans)
        %                                 TumbVis_EC(tr,ch) = mean(squeeze(EEG_EC_sel.data(TumbVisChan,:,tr)));
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %                 if nb_chunks_EO == 1
        %                     for tr = 1:nb_trials
        %                         PSD_EO(:,:,tr,1) = pwelch(squeeze(EEG_EO.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %                         if any(mocapChans)
        %                             TumbVis_EO(tr,1) = mean(squeeze(EEG_EO.data(TumbVisChan,:,tr)));
        %                         end
        %                     end
        %                 else
        %                     for ch = 1:nb_chunks_EO
        %                         EEG_EO_sel = pop_select(EEG, 'time', [ch-1,ch].*study_config.psd.chunks);
        %                         for tr = 1:nb_trials
        %                             PSD_EO(:,:,tr,ch) = pwelch(squeeze(EEG_EO_sel.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %                             if any(mocapChans)
        %                                 TumbVis_EO(tr,ch) = mean(squeeze(EEG_EO_sel.data(TumbVisChan,:,tr)));
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %             otherwise
        %                 error('Not coded yet')
        %         end
        
        TrialsInfo_EC = adaptTrialsInfo(TrialsInfo, TumbVis_EC, nb_chunks_EC);
        PSD_EC = reshape(permute(PSD_EC,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks_EC]);
        TrialsInfo_EO = adaptTrialsInfo(TrialsInfo, TumbVis_EO, nb_chunks_EO);
        PSD_EO = reshape(permute(PSD_EO,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks_EO]);
        
        % Concatenate over subjects:
        TrialsInfo_EC_all = cat(1, TrialsInfo_EC_all, TrialsInfo_EC);
        PSD_EC_all = cat(3,PSD_EC_all,PSD_EC);
        TrialsInfo_EO_all = cat(1, TrialsInfo_EO_all, TrialsInfo_EO);
        PSD_EO_all = cat(3,PSD_EO_all,PSD_EO);
        
        %% Normalize data according to options
        [PSD_EC_norm] = normalizePSD(PSD_EC, TrialsInfo_EC, study_config.psd);
        [PSD_EO_norm] = normalizePSD(PSD_EO, TrialsInfo_EO, study_config.psd);
        
        %% Save data
        opts = study_config.psd;
        opts.SR = EEG.srate;
        folderName = makePSDFolderName(opts);
        if ~exist(fullfile(input_filepath,folderName),'dir')
            mkdir(fullfile(input_filepath,folderName));
        end
        opts.event = study_config.epochs.event;
        opts.phase = 'EC';
        
        fileName = makePSDFileName('labels', subject, opts);
        save(fullfile(input_filepath, fileName), sprintf('TrialsInfo_%s',opts.phase));
        
        fileName = makePSDFileName('psd', subject, opts);
        save(fullfile(input_filepath, folderName, fileName), sprintf('PSD_%s_norm',opts.phase), '-v7.3');
        
        opts.phase = 'EO';
        fileName = makePSDFileName('labels', subject, opts);
        save(fullfile(input_filepath, fileName), sprintf('TrialsInfo_%s',opts.phase));
        
        fileName = makePSDFileName('psd', subject, opts);
        save(fullfile(input_filepath, folderName, fileName), sprintf('PSD_%s_norm',opts.phase), '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end

%% Multi-subjects data
if ~strcmp(MSname,'')
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Normalize data according to options
        [PSD_EC_all_norm] = normalizePSD(PSD_EC_all, TrialsInfo_EC_all, study_config.psd);
        [PSD_EO_all_norm] = normalizePSD(PSD_EO_all, TrialsInfo_EO_all, study_config.psd);
        
        %% Save data
        if ~exist(fullfile(N.searchFolder_4arch_rej_ICcats,folderName),'dir')
            mkdir(fullfile(N.searchFolder_4arch_rej_ICcats,folderName));
        end
        opts.phase = 'EC';
        fileName = makePSDFileName('labels', MSname, opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('TrialsInfo_%s_all',opts.phase));
        
        fileName = makePSDFileName('psd', MSname, opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, folderName, fileName), sprintf('PSD_%s_all_norm', opts.phase), '-v7.3');
        
        opts.phase = 'EO';
        fileName = makePSDFileName('labels', MSname, opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('TrialsInfo_%s_all',opts.phase));
        
        fileName = makePSDFileName('psd', MSname, opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, folderName, fileName), sprintf('PSD_%s_all_norm', opts.phase), '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end
