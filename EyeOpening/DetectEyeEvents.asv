clear all;
config_PIONEER_Tumbler;

eyeDetectionChannel = {'R3Z'};
%eyeDetectionChannels = {'VEOGR','Z2Z','Z3Z','R3Z','L3Z'};
plot_general_level = false;

for subject_ind = subject_inds
    % Launch eeglab
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    % Overwrite subject for testing
    subject_ind = 4;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    EEG = pop_loadset('filename', N.preparedFile,'filepath', N.searchFolder_2);
    %EEG = pop_loadset('filename', N.IClabelledFile,'filepath', N.searchFolder_2arch_rej);
    %EEG.etc = rmfield(EEG.etc, 'filter');
    
    %% Prepare data (inspired from Sharma et al. 2020 method)
    EEG_sel = pop_select(EEG, 'channel', eyeDetectionChannel);
    if strcmp(eyeDetectionChannel{1},'VEOGR')
        % Invert polarity for EOG channel
        EEG_sel.data = -EEG_sel.data;
    end
    EEG_filt1 = custom_filter(EEG_sel, 0.75, 5);
    EEG_epoch1 = pop_epoch(EEG_filt1, {'OpenEyesSoundOn'}, [-7.5,27.5]);
    
    EEG_filt2 = custom_filter(EEG_sel, 0.05, 5);
    EEG_epoch2 = pop_epoch(EEG_filt2, {'OpenEyesSoundOn'}, [-7.5,27.5]);
    
    % General parameters
    opts_peaks.artifact_detection_th = 1000; % in microV
    opts_peaks.timeEO = 0; % in ms
    opts_peaks.timeEC = 20860; % in ms
    opts_peaks.accept_wdw = [-1,5]; % in seconds
    opts_peaks.percMax4Prom = 0;
    opts_peaks.wdw_diff1 = 0.2; % in seconds, for diff calculation
    opts_peaks.wdw_diff2 = 0.4; % in seconds, for diff calculation
    opts_peaks.diffQuant = 75; % for createPeaksTable2
    
    % Train set
    n_trials_train = 20;
    train_trials = randi([1 EEG_epoch1.trials],1,n_trials_train);
    while length(train_trials) ~= length(unique(train_trials))
        train_trials = randi([1 EEG_epoch1.trials],1,n_trials_train);
    end
    
    opts_peaks.quant = 90;
    opts_peaks.userLabel = false;
    opts_peaks.inspectedTrials = train_trials;
    Peaks = createPeaksTable2(EEG_epoch1, EEG_epoch2, opts_peaks);
    
    if plot_general_level
    % Plot
    opts_plot.subject = subject;
    opts_plot.step = 0.25; % in mV
    opts_plot.timeEO = opts_peaks.timeEO; % in ms
    opts_plot.timeEC = opts_peaks.timeEC; % in ms
    %plotEyeEventsOnTrials(EEG_epoch1, EEG_epoch2, Peaks_test, opts_plot)
    plotEyeEventsOnTrials(EEG_epoch1, EEG_epoch2, Peaks, opts_plot)
    end
    
    %% Add events to EEG data
    EEG_final = pop_loadset('filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej_ICcats);
    Trials = EEG_final.etc.TrialData;
    peaks_eo = strcmp(Peaks.AutoCategory,'EyeOpening');
    peaks_ec = strcmp(Peaks.AutoCategory,'EyeClosing');
    ref_lat = find(EEG_epoch1.times == 0);
    for tr = 1:size(Trials,1)
        bl = floor((tr-1)/10)+1;
        tr_bl = mod((tr-1),10)+1;
        ref_ev = intersect(intersect(find(strcmp({EEG_final.event.type},'OpenEyesSoundOn')),...
            findInStructWithEmpties(EEG_final.event, 'block',bl)),...
            findInStructWithEmpties(EEG_final.event, 'trial',tr_bl));
        
        peaks_tr = Peaks.Trial == tr;
        if any(peaks_tr & peaks_eo)
            lat_eo = EEG_final.event(ref_ev).latency + Peaks.Latency(peaks_tr & peaks_eo) - ref_lat;
            
            EEG_final = pop_editeventvals(EEG_final, 'add', {1, 'EyeOpeningDetection',...
                lat_eo/EEG_final.srate, 1/EEG_final.srate,...
                Trials.Block(tr), Trials.Condition{tr}, Trials.Trial(tr),...
                Trials.TrialType{tr}, [], [], []});
        end
        
        if any(peaks_tr & peaks_ec)
            lat_ec = EEG_final.event(ref_ev).latency + Peaks.Latency(peaks_tr & peaks_ec) - ref_lat;
            
            EEG_final = pop_editeventvals(EEG_final, 'add', {1, 'EyeClosingDetection',...
                lat_ec/EEG_final.srate, 1/EEG_final.srate,...
                Trials.Block(tr), Trials.Condition{tr}, Trials.Trial(tr),...
                Trials.TrialType{tr}, [], [], []});
        end
    end
end