%% ONLY CHANGE THESE PARTS!
user = 'sl';
task = 'Tumbler';

%% General foldernames and filenames
[study_config.study_folder, study_config.raw_data_folder, study_config.electrodes_folder,...
    study_config.raw_EEGLAB_data_folder, study_config.preprocessing_folder,...
    study_config.single_subject_analysis_folder, study_config.multi_subject_analysis_folder,...
    study_config.figures_folder, sourceFile, tab] = getMainFoldersNames_PIONEER(user, task);

%% Subjects definition:
study_config.subjects = table2struct(readtable(sourceFile, 'Sheet', tab));
for s = 1:length(study_config.subjects)
    if isempty(study_config.subjects(s).missingData) %|| isnan(study_config.subjects(s).missingData)
        study_config.subjects(s).missingData = {};
    else
        study_config.subjects(s).missingData = split(study_config.subjects(s).missingData,',');
    end
    
    if isnan(study_config.subjects(s).badElectrodes)
        study_config.subjects(s).badElectrodes = {};
    else
        study_config.subjects(s).badElectrodes = split(study_config.subjects(s).badElectrodes,',');
    end
end

included = find(strcmp({study_config.subjects.excluded},'No'));

%% Define subjects to include in the analysis
% Can be customized (specific indexes etc...)
subject_inds = included;
subject_inds = [9]; % WARNING: reads the corresponding excel file line. Or see:
% 1 = P1001
% 2 = P1001-2
% 3 = P1001-3
% 4 = P1004
% 5 = P1002
% 6 = P1002-2
% 7 = P1002-3
% 8 = P1009
% 9 = P1001-4
% 10 = P1004-2
subject_ind = subject_inds(1);
study_config.current_subject = subject_ind;

if strcmp(study_config.subjects(subject_ind).id, 'P1002-2')
    study_config.filenames = {'block002', 'block003', 'block004', 'block005', 'block006', 'block007', 'block008', 'block009'};
    % Enter the number of streams expected for each category in the .xdf file:
    % Categories: {EEG, ET, MOCAP, Events}
    study_config.stream_count = [1,0,0,1];
    % Enter the TYPES to put in each category
    % (used to identify the correct streams in each category)
    study_config.eeg_streams = {'EEG'};
    study_config.eye_tracker_streams = {};
    study_config.rigidbody_streams = {};
    study_config.event_streams = {'Markers'};
    
    % For Mocap streams, specify their order by names for consistency when
    % merging multiple .xdf files
    study_config.rb_streams_order = {};
else
    study_config.filenames = {'block001', 'block002', 'block003', 'block004', 'block005', 'block006', 'block007', 'block008'};
    % Enter the number of streams expected for each category in the .xdf file:
    % Categories: {EEG, ET, MOCAP, Events}
    if any(strcmp(study_config.subjects(subject_ind).id, {'P1002-3','P1009','P1001-4','P1004-2'}))
        study_config.stream_count = [1,0,2,1];
    else
        study_config.stream_count = [2,0,2,1];
    end
    % Enter the TYPES to put in each category
    % (used to identify the correct streams in each category)
    study_config.eeg_streams = {'EEG'};
    study_config.eye_tracker_streams = {};
    study_config.rigidbody_streams = {'MoCap'};
    study_config.event_streams = {'Markers'};
    
    % For Mocap streams, specify their order by names for consistency when
    % merging multiple .xdf files
    study_config.rb_streams_order = {'RigidBody1', 'RigidBody2'};
end

if any(strcmp(study_config.subjects(subject_ind).id, {'P1002-2','P1002-3','P1009','P1001-4','P1004-2'}))
    study_config.eeg_streams_order = {'EE225-000000-000430'};
    study_config.recording_unit = 'Volt';
    
    % Enter channels that you did not use at all:
    % For ANT data, the last 2 channels are 'trigger' and 'counter' and shoud be removed
    % Depends on the total number of channels per amplifier, specify for each amplifier
    study_config.channels_to_remove = {65:66};
    
    % Enter EOG channel names here:
    study_config.eog_channels  = {'EOG'};
    % Enter Bipolar channel names here:
    study_config.bip_channels  = {};
    
    %% Channel location files
    % standard locations for the cap used
    study_config.capName = 'CA-208';
    study_config.channel_locations_filename = sprintf('%s.elc', study_config.capName);
    % (optional) If you have an electrodes digitization file for each participant:
    % (the name of the participant will be added at the beginning of this string)
    study_config.indiv_channel_locations_filename = '';
else
    % Specific to ANT: 128ch creates 2 EEG streams (one per amplifier) and it
    % is important to know which one should be imported first to assign the
    % right channel labels (indicate stream NAMES here)
    study_config.eeg_streams_order = {'EE225-000000-000430', 'EE225-000000-000433'};
    study_config.recording_unit = 'Volt';
    
    % Enter channels that you did not use at all:
    % If none:
    %study_config.channels_to_remove = [];
    % By labels:
    %study_config.channels_to_remove = {'N29' 'N30' 'N31'};
    % For ANT data, the last 2 channels are 'trigger' and 'counter' and shoud be removed
    % Depends on the total number of channels per amplifier, specify for each amplifier
    study_config.channels_to_remove = {65:67,65:66};
    
    % Enter EOG channel names here:
    study_config.eog_channels  = {'VEOGR'};
    % If none:
    %study_config.eog_channels = {};
    % Enter Bipolar channel names here:
    study_config.bip_channels  = {};
    
    %% Channel location files
    % standard locations for the cap used
    study_config.capName = 'CW-04852';
    study_config.channel_locations_filename = sprintf('%s_wavegard.elc', study_config.capName);
    % (optional) If you have an electrodes digitization file for each participant:
    % (the name of the participant will be added at the beginning of this string)
    study_config.indiv_channel_locations_filename = '';
end

if strcmp(task,'Tumbler')
    study_config.channel_selection = {};
    % OR Only use a subset of channels for the whole analysis (simulate a sparser dataset)
    %     study_config.channel_selection = {'VEOGR','Z2Z','L3Z','R3Z','Z5Z','L6Z','R6Z','Z7Z',...
    %         'L8Z','R8Z','Z9Z','L10Z','R10Z','Z11Z','L12Z','R12Z','Z13Z','L14Z','R14Z','Z15Z',...
    %         'L16Z','R16Z','Z17Z','L18Z','R18Z','Z19Z','L20Z','R20Z','L2A','L4A','L6A','L2C','L4C','L6C','L8C',...
    %         'L2E','L4E','L6E','L8E','L10E','L2G','L4G','L5G','L7G','R2A','R4A','R6A','R2C','R4C','R6C','R8C',...
    %         'R2E','R4E','R6E','R8E','R10E','R2G','R4G','R5G','R7G'};
else
    study_config.channel_selection = {};
end

%% Definition of the global preprocessing architecture:
% 'bemobil' = architecture developped in Gramann lab
% 'simple' = simple pipeline for quick results
study_config.globalArchitecture = 'bemobil';

%% Preprocessing
% Electrodes that should be excluded right away (known recording issue)
% One line per subject concerned
study_config.moveElecInwards = 0; % in mm
%study_config.mocap_lowpass = 6;
%study_config.rigidbody_derivatives = 2;
% Time buffer to include before and after each trial selected for analysis
study_config.trialBuffer = 3; % in seconds
% Maximal length of NaN series to interpolate (fillNaNs function)
study_config.maxNans2Replace = 3; % in samples
% Resampling frequency
study_config.resample_freq = 250; % in Hz

% Filters depending on the step:
study_config.filterPreProc.low_cut_off = 1.5; % in Hz, put [] for no HP filter
study_config.filterPreProc.high_cut_off = []; % in Hz, put [] for no LP filter
study_config.filterICLabel.low_cut_off = 0.5; % in Hz, put [] for no HP filter
study_config.filterICLabel.high_cut_off = []; % in Hz, put [] for no LP filter
study_config.filterAnalysis.low_cut_off = 0.5; % in Hz, put [] for no HP filter
study_config.filterAnalysis.high_cut_off = 45; % in Hz, put [] for no LP filter

study_config.badSampsRejection = 'app'; % 'manual', 'app', 'asr'
study_config.do_second_tempRej = false; % Keep to false

%% Parameters for ASR pipeline:
study_config.ASR_use = 'reject'; % 'rewrite' or 'reject'
study_config.burst_crit = 25; % Burst criterion for ASR. Put [] for default

%% Parameters for the APP pipeline:
study_config.APP.censorBiweight = 6; % censor value, it corresponds to a certain number of std for a normal distribution
% c=6 is 4 std; c=7.5 is 5 std; c=9 is 6 std.
study_config.APP.z_criterion = 3.5; % rejection criterion according to literature
study_config.APP.inner_fence = 1.5; % inner_fence criterion for extreme outliers

% from M. Hubert, E. Vandervieren, An adjusted boxplot for skewed distributions,
% Computational Statistics & Data Analysis, Volume 52, Issue 12, 2008, Pages 5186-5201, ISSN 0167-9473,
% https://doi.org/10.1016/j.csda.2007.11.008.
%study_config.APP.skew_side_limit = 3;
%study_config.APP.opp_side_limit = 4;

% from Janir:
%study_config.APP.skew_side_limit = 4;
%study_config.APP.opp_side_limit = 3.5;

study_config.APP.skew_side_limit = 2;
study_config.APP.opp_side_limit = 4;

%% Parameters for bemobil pipeline:
%%%% Line noise removal technique:
% 'cleanLine' = use the Cleanline plugin (as in APP)
% 'cleanLinePREP' = use the PREP plugin, CleanLineNoise function
study_config.lineNoiseRemoval_method = 'cleanLinePREP';

%% Different ICA methods:
% 'runica'
% 'amica'
study_config.ICAmethod = 'amica';
study_config.max_threads = 4;
study_config.num_models = 1;
study_config.max_iter = 2000;

%% Dipole fitting:
study_config.dipfit.doDipoleFitting = true;
% Indicate a single transformation for all participants. Keep to [] if you have individual data.
% For this transform, check coregistration_CW04852.m file
study_config.dipfit.transform = [-0.0091,-36.3402,-55.9678,...
    0.23856,0.00766,-1.59112,...
    0.9298,1.1298,1.0526];
study_config.dipfit.use_fiducials = false;
study_config.dipfit.residualVariance_threshold = 100;
study_config.dipfit.do_remove_outside_head = 'off';
study_config.dipfit.number_of_dipoles = 1;

%% IC_label
% classes: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}
% The thresholds are extracted for the latest released paper about IClabel
% They are extracted from Table 3 to optimize accuracy of the classifier (Test)

% Different ICLabel options:
% 'lite' ; [0.53,0.17,0.06,0.10,0.42,0.15,0.29];
% 'default'; [0.35,0.30,0.04,0.03,0.84,0.05,0.26];
study_config.iclabel = 'default';
study_config.ICdetect_thresholds = [0.35,0.30,0.04,0.03,0.84,0.05,0.26];
study_config.ICcategorization = 'auto'; % 'auto' or 'manual'
study_config.ICselection = 'manual'; % 'auto' or 'manual'
study_config.cats2keep = {'Brain','BrainWithNoise'}; %'OccBrain','Brain','BrainWithNoise'

%% Epoching
study_config.epochs.event = 'EyesOpening';
study_config.epochs.window = 'fixed'; % 'full' = based on durations, 'fixed' = using limits_wdw
study_config.epochs.limits_wdw_base = [0,1]; % in seconds
study_config.epochs.limits_wdw = [-5,20]; % in seconds
study_config.epochs.bumper = 1; % in seconds (additional time to include before and after the trial
% to prevent artifacts in ERSP computation = padding)
% study_config.epochs.event = 'TumblerVisible';
% study_config.epochs.window = 'fixed'; % 'full' = based on durations, 'fixed' = using limits_wdw
% study_config.epochs.limits_wdw_base = [0,1]; % in seconds
% study_config.epochs.limits_wdw = [-0.5,1]; % in seconds
% study_config.epochs.bumper = 0.5; % in seconds (additional time to include before and after the trial
% % to prevent artifacts in ERSP computation = padding)

%% Amplitude features
study_config.amps.chunks = 0;
study_config.amps.normStyle = 'none';

%% Computing Power spectral densities
study_config.psd.method = 'pwelch';
study_config.psd.output = 'psd'; % For pwelch: 'psd' or 'power' (see pwelch spectrumtype option)
study_config.psd.window = study_config.resample_freq;
study_config.psd.overlap = 0.8*study_config.psd.window;
study_config.psd.chunks = 0; % Divide epoch data in chunks (multiple chunks for one trial)
% 0: no chunking
% any other value: chunk length in seconds
% Supported: 'pwelch'
study_config.psd.FoI = 1:40; % frequencies of interest, in Hz, make sure there is a sufficient margin with respect to filters
study_config.psd.normStyle = 'acrossTrials'; % How which dimensions to apply normalization
% Supported: 'none', 'acrossChans', 'acrossTrials', 'acrossChans&Trials'
study_config.psd.normTrialsGroup = {'perSubject', 'perCondition'};
%study_config.psd.normTrialsGroup = {'perSubject'};
% How trials should be grouped together for normalization (can be multiple options at the same time)
% Supported: perSubject, perCondition, perBlock, perTrialType, perAnswer
study_config.psd.normModel = 'gain';
% Supported: 'additive', 'gain'

%% Classifier specifications
study_config.class.model = 'linear'; % 'linear', 'diaglinear', 'svm' or 'logistic'
study_config.class.condAnalysis = 'separate'; % only 'separate' for now: analyze GogglesON and GogglesOFF separately
study_config.class.contrast = 'Answer'; % 'TrialType': build a classifier to distinguish between trial types
% 'Answer': build a classifier to distinguish between subject answers (Present vs Absent)
% 'Visibility_binary': build a classifier to distinguish between absence
% and presence on the camera FoV (with and without object)
% 'Visibility_binary2': build a classifier to distinguish between presence
% on the camera FoV and no object
study_config.class.thresh_vis_bin = 12; % in percent, to separate between visible and invisible tumbler trials
study_config.class.feat2test = [1:10]; % vector of best features to test
%study_config.class.feat2test = [1,5,10,15,20,25]; % vector of best features to test
study_config.class.nbfolds = 'loo'; % Number of folds for the cross-validation analysis.
% Set to 'loo' if you want a Leave-one-out cross validation.

%% Defining features for classifier
study_config.feat.type = 'psd'; % 'psd' or 'amp'
%%%%% Configurations for 'amp':
% Specific to Tumbler2020:
% study_config.feat.folderName = 'occChans-allTimes';
% study_config.feat.chans = {'O1','Oz','O2'}; % 'all' or cell of channels
% study_config.feat.times = 'all'; % 'all' or range [min, max] or vector

% study_config.feat.folderName = 'occChans-allTimes'; % Name for distinguishing analyses by folders
% study_config.feat.chans = {'Z17Z', 'Z18Z', 'Z19Z',...
%     'Z17L', 'Z18L', 'L18Z',...
%     'Z17R', 'Z18R', 'R18Z'}; % 'all' or cell of channels
% study_config.feat.times = 'all'; % 'all' or range [min, max] or vector

% study_config.feat.folderName = 'occChansExt-allTimes'; % Name for distinguishing analyses by folders
% study_config.feat.chans = {'Z16Z', 'Z17Z', 'Z18Z', 'Z19Z', 'Z20Z',...
%     'Z16L', 'Z17L', 'Z18L', 'Z19L',...
%     'L17Z', 'L18Z', 'L19Z', 'L20Z', 'L17L', 'L18L', 'L19L',...
%     'Z16R', 'Z17R', 'Z18R', 'Z19R',...
%     'R17Z', 'R18Z', 'R19Z', 'R20Z', 'R17R', 'R18R', 'R19R'}; % 'all' or cell of channels
% study_config.feat.times = 'all'; % 'all' or range [min, max] or vector

%%%%% Configurations for 'psd':
% study_config.feat.folderName = 'allChans-allFreqs'; % Name for distinguishing analyses by folders
% study_config.feat.chans = 'all'; % 'all' or cell of channels
% study_config.feat.freqs = 'all'; % 'all' or range [min, max] or vector

study_config.feat.folderName = 'allChans-alpha';
study_config.feat.chans = 'all'; % 'all' or cell of channels
study_config.feat.freqs = [8,14]; % 'all' or range [min, max] or vector

% Specific to Tumbler2020:
% study_config.feat.folderName = 'occChans-alpha';
% study_config.feat.chans = {'O1','Oz','O2'}; % 'all' or cell of channels
% study_config.feat.freqs = [8,14]; % 'all' or range [min, max] or vector

% study_config.feat.folderName = 'occChans-alpha';
% study_config.feat.chans = {'Z17Z', 'Z18Z', 'Z19Z',...
%     'Z17L', 'Z18L', 'L18Z',...
%     'Z17R', 'Z18R', 'R18Z'}; % 'all' or cell of channels
% study_config.feat.freqs = [8,14]; % 'all' or range [min, max] or vector

% study_config.feat.folderName = 'Frontal-allFreqs';
% study_config.feat.chans = {'L5Z','Z5Z','R5Z','L6Z','Z6Z','R6Z'}; % 'all' or cell of channels
% study_config.feat.freqs = 'all'; % 'all' or range [min, max] or vector

% study_config.feat.folderName = 'occChansExt-allFreqs';
% study_config.feat.chans = {'Z16Z', 'Z17Z', 'Z18Z', 'Z19Z', 'Z20Z',...
%     'Z16L', 'Z17L', 'Z18L', 'Z19L',...
%     'L17Z', 'L18Z', 'L19Z', 'L20Z', 'L17L', 'L18L', 'L19L',...
%     'Z16R', 'Z17R', 'Z18R', 'Z19R',...
%     'R17Z', 'R18Z', 'R19Z', 'R20Z', 'R17R', 'R18R', 'R19R'}; % 'all' or cell of channels
% study_config.feat.freqs = 'all'; % 'all' or range [min, max] or vector

% study_config.feat.folderName = 'occChansExt-alpha';
% study_config.feat.chans = {'Z16Z', 'Z17Z', 'Z18Z', 'Z19Z', 'Z20Z',...
%     'Z16L', 'Z17L', 'Z18L', 'Z19L',...
%     'L17Z', 'L18Z', 'L19Z', 'L20Z', 'L17L', 'L18L', 'L19L',...
%     'Z16R', 'Z17R', 'Z18R', 'Z19R',...
%     'R17Z', 'R18Z', 'R19Z', 'R20Z', 'R17R', 'R18R', 'R19R'}; % 'all' or cell of channels
% study_config.feat.freqs = [8,14]; % 'all' or range [min, max] or vector

%% Filenames suffix
study_config.merged_filename = 'allBlocks.set'; % In all pipelines
study_config.prepared_filename = 'prepared.set'; % In all pipelines
study_config.ASRin_filename = 'inputASR.set'; % In simple pipeline
study_config.ASRout_filename = 'outputASR.set'; % In simple pipeline
study_config.BadChansRemoved_filename = 'inter_avRef.set'; % In bemobil pipeline
study_config.beforeICA_filename = 'cleanedForICA.set'; % In all pipelines
study_config.icaOutput_filename = 'postICA.set'; % In all pipelines
study_config.dipolesFitted_filename = 'postICA_withdips.set'; % In all pipelines
study_config.icLabelled_filename = 'automaticIClabels.set'; % In all pipelines
study_config.icaSelect_filename = 'cleaned_with_ICA.set'; % In all pipelines
study_config.postICA_tempRej_filename = 'final.set'; % In all pipelines
study_config.epoched_filename = 'epoched.set'; % In all pipelines
study_config.base_epoched_filename = 'epoched_baselines.set'; % In all pipelines

study_config.interpolated_avRef_filename = 'interpolated_avRef.set'; %
study_config.filtered_filename = 'filtered.set';
study_config.without_bad_temp_segments = 'badTempSegmentsRemoved.set';
study_config.without_eyes = 'no_eyes.set';
study_config.ica_filename_output = 'postICA.set';
study_config.warped_dipfitted_filename = 'warped_dipfitted.set';
study_config.copy_weights_interpolate_avRef_filename = 'interp_avRef_ICA.set';
study_config.categorizedICs = 'ICLabel_computed.set';
study_config.single_subject_cleaned_ICA_filename = 'cleaned_with_ICA.set';
