clear all;
config_PIONEER_Tumbler;

%% Epoch the data for further analyses
% Options for this script
recompute = true;
do_baseline = false;

N = makeFolderFileNames(study_config, study_config.subjects(study_config.current_subject).id);
fname = 'MergedDataBase_forEEG.mat';
load(fullfile(N.searchFolder_2, fname));

if ~exist('EEG', 'var')
    launchEEGLAB
end

%% Goal epoch structure
% Define a common epoch window for all:
switch lower(study_config.epochs.window)
    case 'full'
        error('Notready')
        % Definition depending on the durations (assumes timewarping later)
        obsDurations_ms_usedTrials = obsDurations_ms([allObservations.completeEEG]);
        
        % Take a 1s margin with respect to min and max values
        epoch_wdw = [-1,1].*study_config.epochs.bumper +...
            [-max(obsDurations_ms_usedTrials(:,1)), max(sum(obsDurations_ms_usedTrials(:,2:end),2))]/1000;
        
    case 'fixed'
        % Definition with no timewarping implied
        epoch_wdw = [-1,1].*study_config.epochs.bumper + study_config.epochs.limits_wdw;
    otherwise
        error('Unknown window type')
end

%% Load the data, filter the data, epoch the data, keep only the right epochs, save the data
subjects = {study_config.subjects(subject_inds).id};
n_sbjs = length(subjects);
for s = 1:n_sbjs
    disp(['Subject ' subjects{s}]);
    
    % Special case for P1001-3 (non homogeneous type of recording)
    %     if strcmp(subjects{s},'P1001-3')
    %         study_config.eog_channels = {};
    %         study_config.channel_selection = setdiff(study_config.channel_selection,{'VEOGR'});
    %     end
    
    if any(strcmp({'P1001', 'P1004'},subjects{s})) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,20]; % in seconds
    elseif any(strcmp({'P1001old-1', 'P1001old-2','P1001-2', 'P1001-3', 'P1002-2', 'P1002-3', 'P1009','P1001-4','P1004-2'},subjects{s})) &&...
            strcmp(study_config.epochs.event,'EyesOpening') &&...
            strcmp(study_config.epochs.window, 'fixed')
        study_config.epochs.limits_wdw = [-5,15]; % in seconds
    end
    
    N = makeFolderFileNames(study_config, subjects{s});
    output_filepath = N.searchFolder_3arch_rej_ICcats;
    if ~exist(output_filepath, 'dir')
        mkdir(output_filepath);
    end
    
    if ~recompute && exist(fullfile(output_filepath,N.epochedFile),'file')...
            && (exist(fullfile(output_filepath,N.baselineEpochedFile),'file') || ~do_baseline)
        continue
    end
    
    EEG = pop_loadset('filename', N.postLabelingFile, 'filepath', N.searchFolder_2arch_rej_ICcats);
    TumbVisChan = strcmp({EEG.chanlocs.labels}, 'TumblerVisibility');
    
    if strcmp(task,'Tumbler')
        % Filter the dataset with the desired upper and lower frequencies
        % (definitive changes before epoching)
        lowcutoff = study_config.filterAnalysis.low_cut_off;
        highcutoff = study_config.filterAnalysis.high_cut_off;
        fprintf('Filtering between %.1f Hz and %.1f Hz...\n', lowcutoff, highcutoff)
        [EEG_filt] = custom_filter(EEG, lowcutoff, highcutoff);
        % No need for line Noise removal here (LP filtering)
        if any(TumbVisChan)
            % Revert the filter for TumbVisChan
            EEG_filt.data(TumbVisChan,:) = EEG.data(TumbVisChan,:);
        end
        %clear EEG
    else
        % data is already filtered in this case
        EEG_filt = EEG;
    end
    
    % Take the trials of interest for this subject only
    switch study_config.epochs.event
        case 'EyesOpening'
            trial_selection = strcmp(MergedDataBase.ID, subjects{s}) & ...
                MergedDataBase.completeEEG_full == 100;
            Trials_subj = MergedDataBase(trial_selection,:);
            ev_subj = MergedDataBase.urevent_seq(trial_selection,2);
            %Durations_ms_subj = MergedDataBase_withOutliers.duration_ms_seq(trial_selection,2);
            if do_baseline
                trial_selection_base = strcmp(MergedDataBase.ID, subjects{s}) & ...
                    MergedDataBase.completeEEG_full == 100;
                Base_subj = MergedDataBase(trial_selection_base,:);
                ev_base_subj = MergedDataBase.urevent_seq(trial_selection_base,1);
            end
        case 'TumblerVisible'
            Trials_subj = createTumbVisEvents(EEG, TumbVisChan, MergedDataBase(strcmp(MergedDataBase.ID, subjects{s}),:));
            for tr = 1:size(Trials_subj,1)
                EEG_filt = pop_editeventvals(EEG_filt, 'add', {1, 'TumblerVisible',...
                    Trials_subj.TumbVisEventLat(tr)/EEG_filt.srate, 1/EEG_filt.srate,...
                    Trials_subj.Block(tr), Trials_subj.Condition{tr}, Trials_subj.Trial(tr),...
                    Trials_subj.TrialType{tr}, [], [], []});
            end
        otherwise
            error('I don''t know this event');
    end
    
    %% Epoch Data:
    boundaryEvents_mask = strcmp({EEG_filt.event.type}, 'boundary');
    % delete boundary events
    EEG2 = pop_editeventvals(EEG_filt, 'delete', find(boundaryEvents_mask));
    clear EEG_filt
    
    % add a bunch of zeros at the end of the data to avoid reject the last
    % epochs because of insufficient length
    EEG2.pnts = EEG2.pnts + round(epoch_wdw(2)*EEG2.srate);
    EEG2.xmax = EEG2.xmax + epoch_wdw(2)*1000;
    EEG2.times = EEG2.xmin:(1000/EEG2.srate):EEG2.xmax;
    EEG2.data = cat(2, EEG2.data, zeros(EEG2.nbchan, round(epoch_wdw(2)*EEG2.srate)));
    EEG2.icaact = [];
    EEG2 = eeg_checkset(EEG2);
    
    if strcmp(task,'Tumbler')
        % Add a channel corresponding to bad samples
        EEG2.nbchan = EEG2.nbchan + 1;
        % Get rejected samples
        switch study_config.badSampsRejection
            case 'app'
                rejSamps = EEG2.etc.APP.rejectedSamples;
            case 'asr'
                rejSamps = EEG2.etc.ASR.rejectedSamples;
            case 'autoMoBI'
                rejSamps = EEG2.etc.autoMoBI.rejectedSamples;
        end
        EEG2.data = cat(1, EEG2.data, [rejSamps, ones(1,round(epoch_wdw(2)*EEG2.srate))]);
        EEG2.chanlocs(EEG2.nbchan).type = 'MASK';
        EEG2.chanlocs(EEG2.nbchan).labels = 'rejSamps';
        %     EEG2.chanlocs(EEG2.nbchan) = struct('type', 'MASK', 'labels', 'rejSamps', 'ref', '', 'urchan', [],...
        %         'X', [], 'Y', [], 'Z', [], 'unit', [], 'sph_theta', [], 'sph_phi', [], 'sph_radius', [], 'theta', [], 'radius', []);
    end
    
    % Epoch the dataset according to observation events:
    switch study_config.epochs.event
        case 'EyesOpening'
            epochEvents = find(logical(sum(ev_subj == [EEG2.event.urevent],1)));
            EEG_epoched = pop_epoch(EEG2, {}, epoch_wdw, 'eventindices', epochEvents, 'epochinfo', 'yes'); % time is indicated in seconds
        case 'TumblerVisible'
            epochEvents = find(strcmp({EEG2.event.type},'TumblerVisible'));
            EEG_epoched = pop_epoch(EEG2, {}, epoch_wdw, 'eventindices', epochEvents, 'epochinfo', 'yes'); % time is indicated in seconds
    end
    
    EEG_epoched.etc.epochInfo.Trials = Trials_subj;
    %EEG_epoched.etc.epochInfo.Durations_ms = Durations_ms_subj;
    if strcmp(task,'Tumbler')
        % copy bad samples information and remove the corresponding channel
        EEG_epoched.etc.epochInfo.rejectedSamples = logical(squeeze(EEG_epoched.data(end,:,:)))';
        EEG_epoched.nbchan = EEG_epoched.nbchan - 1;
        EEG_epoched.data = EEG_epoched.data(1:end-1,:,:);
        EEG_epoched.chanlocs = EEG_epoched.chanlocs(1:end-1);
    end
    EEG_epoched.icaact = [];
    EEG_epoched = eeg_checkset(EEG_epoched);
    
    pop_saveset(EEG_epoched, 'filename', N.epochedFile,'filepath', output_filepath);
    clear EEG_epoched
    
    if do_baseline
        % Epoch the dataset according to baseline events:
        epochEvents_base = find(logical(sum(ev_base_subj == [EEG2.event.urevent],1)));
        % Dark period is always timed to 1s in the experiment
        epoch_wdw_base = [-1,1].*study_config.epochs.bumper + study_config.epochs.limits_wdw_base;
        EEG_epoched_base = pop_epoch(EEG2, {}, epoch_wdw_base, 'eventindices', epochEvents_base, 'epochinfo', 'yes'); % time is indicated in seconds
        
        EEG_epoched_base.etc.epochInfo.Baselines = Base_subj;
        % copy bad samples information and remove the corresponding channel
        EEG_epoched_base.etc.epochInfo.rejectedSamples = logical(squeeze(EEG_epoched_base.data(end,:,:)))';
        EEG_epoched_base.nbchan = EEG_epoched_base.nbchan - 1;
        EEG_epoched_base.data = EEG_epoched_base.data(1:end-1,:,:);
        EEG_epoched_base.chanlocs = EEG_epoched_base.chanlocs(1:end-1);
        EEG_epoched_base.icaact = [];
        EEG_epoched_base = eeg_checkset(EEG_epoched_base);
        
        pop_saveset(EEG_epoched_base, 'filename', N.baselineEpochedFile, 'filepath', output_filepath);
        clear EEG_epoched_base
    end
    clear EEG2
end