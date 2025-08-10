%clear all;
config_PIONEER_VEP;
study_config.epochs.limits_wdw = [0,10];
warning('Hardcoding epoch limit window, actually depends on subject (was 30, now 10)');

% Choose MultiSubject name:
% if isempty(study_config.channel_selection) && all(any(subject_inds' == [1,4],1))
%     MSname = 'Allsubjs';
% elseif ~isempty(study_config.channel_selection) && all(any(subject_inds' == [1:4],1))
%     MSname = 'Allsubjs';
% elseif length(subject_inds) == 3 && all(contains({study_config.subjects(subject_inds).id},'P1001'))
%     MSname = 'P1001all';
% else
%     MSname = '';
% end
MSname = '';

PSD_all = [];
TrialsInfo_all = [];
for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    input_filepath = N.searchFolder_3arch_rej_ICcats;
    
    if strcmp(study_config.epochs.event, 'TrialStart')
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath);
        
        % Select channels
        if ~isempty(study_config.eog_channels)
            EEG = pop_select(EEG, 'nochannel', study_config.eog_channels);
        end       
        
        nb_chans = EEG.nbchan;
        nb_trials = EEG.trials;
        freqs = study_config.psd.FoI;
        nb_freqs = length(freqs);
        
        TrialsInfo = EEG.etc.epochInfo.Trials;
        TrialsInfo = TrialsInfo(:,1:5);        
        
        [PSD, nb_chunks] = computePSD(EEG, study_config);
        TrialsInfo = adaptTrialsInfo(TrialsInfo, [], nb_chunks);
        PSD = reshape(permute(PSD,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks]);
        
        % Concatenate over subjects:
        TrialsInfo_all = cat(1, TrialsInfo_all, TrialsInfo);
        PSD_all = cat(3,PSD_all,PSD);
        
        %% Normalize data according to options
        [PSD_norm] = normalizePSD(PSD, TrialsInfo, study_config.psd);
        
        %% Save data
        opts = study_config.psd;
        opts.SR = EEG.srate;
        folderName = makePSDFolderName(opts);
        if ~exist(fullfile(input_filepath,folderName),'dir')
            mkdir(fullfile(input_filepath,folderName));
        end
        opts.event = study_config.epochs.event;
        
        delims = strfind(N.epochedFile,'_');
        opts.epochWdW = N.epochedFile(delims(end-1)+1:delims(end)-1);
        fileName = makePSDFileName('labels', subject, opts);
        save(fullfile(input_filepath, fileName), 'TrialsInfo');
        
        fileName = makePSDFileName('psd', subject, opts);
        save(fullfile(input_filepath, folderName, fileName), 'PSD_norm', '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end

%% Multi-subjects data
if ~strcmp(MSname,'')
    if strcmp(study_config.epochs.event, 'EyesOpening')
        %% Normalize data according to options
        [PSD_all_norm] = normalizePSD(PSD_all, TrialsInfo_all, study_config.psd);
        
        %% Save data
        if ~exist(fullfile(N.searchFolder_4arch_rej_ICcats,folderName),'dir')
            mkdir(fullfile(N.searchFolder_4arch_rej_ICcats,folderName));
        end
        
        fileName = makePSDFileName('labels', MSname, opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), 'TrialsInfo_all');
        
        fileName = makePSDFileName('psd', MSname, opts);
        save(fullfile(N.searchFolder_4arch_rej_ICcats, folderName, fileName), 'PSD_all_norm', '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end
