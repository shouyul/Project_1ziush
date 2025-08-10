%% processing loop
clear all;
config_PIONEER_Tumbler;

skipImport = false; % Boolean to avoid running Mobilab step
overwriteImport = false; % Boolean to force Mobilab step to happen anyway
skipPrep = false; % Boolean to avoid preparation
overwritePrep = false | overwriteImport; % Boolean to force the preparation to happen anyway
skipPreproc = false; % Boolean to avoid preprocessing
overwritePreproc = false | overwritePrep; % Boolean to force the preprocessing to happen anyway
overwritePercVis = false | overwritePreproc; % Boolean to force the computation of tumbler visibility to happen anyway

for subject_ind = subject_inds
    % Launch eeglab
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    % Overwrite subject for testing
    %subject_ind = 2;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    % Special case for P1001-3 (non homogeneous type of recording)
    if strcmp(subject,'P1001-3')
        study_config.channel_locations_filename = sprintf('%s_wavegard_noEOG.elc', study_config.capName);
        study_config.eog_channels = {};
        study_config.channels_to_remove = {65:66,64:65};
        study_config.channel_selection = setdiff(study_config.channel_selection,{'VEOGR'});
    end
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    %% Importing XDF files
    if ~skipImport && (~exist([N.searchFolder_1 N.postimportFile],'file') || overwriteImport)
        % load xdf files and process them with mobilab, export to eeglab, split MoBI and merge all conditions for EEG
        [~, ~, ~] = xdf2set_PIONEER(ALLEEG, CURRENTSET, subject_ind, study_config, overwriteImport);
    end
    
    % Load anyway because EEG_merged output from xdf2set is only EEG data
    if ~skipPrep && (~exist([N.searchFolder_2 N.prepared_withMoCapFile],'file') || overwritePrep)
        %continue
        EEG_merged = pop_loadset('filename', N.postimport_withMoCapFile, 'filepath', N.searchFolder_1);
        [ALLEEG, EEG_merged, CURRENTSET] = eeg_store(ALLEEG, EEG_merged, CURRENTSET);
    end
    
    %continue
    %% Initial preparation of the data (definitive changes)
    if ~skipPrep && (~exist([N.searchFolder_2 N.prepared_withMoCapFile],'file') || overwritePrep)
        EEG_merged = changeUnit2MicroVolt(EEG_merged, study_config);
        
        % Fill NaNs (if comprised between valid samples:
        EEG_merged = fillNaNs(EEG_merged, study_config);
        
        % Check events and report for missing data
        EEG_merged = events_check_Tumbler(EEG_merged, study_config);
        EEG_merged = BehaviorReportTumbler(EEG_merged, study_config);
        
        % Resampling
        EEG_merged = pop_resample(EEG_merged, study_config.resample_freq);
        EEG_merged = eeg_checkset(EEG_merged);
        
        % Remove NaN regions
        eeg_chans = strcmp({EEG_merged.chanlocs.type}, 'EEG') | strcmp({EEG_merged.chanlocs.type}, 'EOG');
        nanTimePoints = sum(isnan(EEG_merged.data(eeg_chans,:)),1)>0;
        EEG_selected = pop_select(EEG_merged, 'nopoint', mask2intervals(nanTimePoints));
        
        % Select Complete Trials only (the dataset remains 'continuous' but remove data between trials)
        %trials_intervals = getIntervals(EEG_merged, 'fullTrial', study_config.trialBuffer, true);
        %EEG_selected = pop_select(EEG_merged,'time',trials_intervals);
        
        % Check there is no more NaN
        if ~sum(isnan(EEG_selected.data(eeg_chans,:)),'all')==0
            error('Still NaNs in the selected data set')
        end
        EEG_selected = eeg_checkset(EEG_selected);
        
        if ~isempty(study_config.channel_selection)
            EEG_selected = pop_select(EEG_selected,'channel',...
                union(study_config.channel_selection, ...
                {EEG_selected.chanlocs(contains({EEG_selected.chanlocs.type}, 'MOCAP')).labels}));
        end
        
        EEG_prepared = pop_saveset(EEG_selected, 'filename', N.prepared_withMoCapFile,'filepath', N.searchFolder_2);
        [ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
    elseif ~skipPreproc && (~exist([N.searchFolder_2 N.MoCapPreparedFile],'file') || overwritePreproc)
        EEG_prepared = pop_loadset('filename', N.prepared_withMoCapFile,'filepath', N.searchFolder_2);
        [ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
    else % No need to load the data
    end
    % Save RAM
    clear EEG_merged
    
    if ~skipPreproc && (~exist([N.searchFolder_2 N.MoCapPreparedFile],'file') || overwritePreproc)
        %% Extract MoCap data only
        mocap_chans = strcmp({EEG_prepared.chanlocs.type}, 'MOCAP');
        EEG_mocap = pop_select(EEG_prepared, 'channel', find(mocap_chans));
        
        EEG_mocap = dealWithMOCAPNaNs(EEG_mocap, 'RB1');
        EEG_mocap = dealWithMOCAPNaNs(EEG_mocap, 'RB2');
        
        % Check there is no more NaN
        if ~sum(isnan(EEG_mocap.data),'all')==0
            error('Still NaNs in the selected data set')
        end
        EEG_mocap = eeg_checkset(EEG_mocap);
        
        % Identify individual markers within rigid bodies
        tol_params = struct();
        tol_params.ori_good = 10; % in degrees
        tol_params.dist_good = 0.1; % in perc (0.1 = 10%)
        tol_params.ori_bad = 5; % in degrees
        tol_params.dist_bad = 0.05; % in perc (0.1 = 10%)
        tol_params.dims = 0.1; % Not used for RB1
        [EEG_mocap_RB1, bad_rb1] = findMarkersRole(EEG_mocap, 'RB1', 4, tol_params);
        % Interpolate bad data
        EEG_mocap_interp = EEG_mocap_RB1;
        RB1_chans = contains({EEG_mocap_RB1.chanlocs.labels},'RB1');
        EEG_mocap_interp.data(RB1_chans,:) = interp1(EEG_mocap_RB1.times(~bad_rb1),...
            EEG_mocap_RB1.data(RB1_chans,~bad_rb1)', EEG_mocap_RB1.times, 'spline')';
        
        
        tol_params.ori_good = 25; % in degrees
        tol_params.dist_good = 0.5; % Can be more tolerant since dimensions are checked later (avoid being prompted too often)
        tol_params.ori_bad = 5; % in degrees
        tol_params.dist_bad = 0.05; % in perc (0.1 = 10%)
        tol_params.dims = 0.1; % in perc (0.1 = 10%)
        [EEG_mocap_RB2, bad_rb2] = findMarkersRole(EEG_mocap, 'RB2', 4, tol_params);
        % Interpolate bad data
        RB2_chans = contains({EEG_mocap_RB2.chanlocs.labels},'RB2');
        EEG_mocap_interp.chanlocs(RB2_chans) = EEG_mocap_RB2.chanlocs(RB2_chans);
        EEG_mocap_interp.data(RB2_chans,:) = interp1(EEG_mocap_RB2.times(~bad_rb2),...
            EEG_mocap_RB2.data(RB2_chans,~bad_rb2)', EEG_mocap_RB2.times, 'spline')';
        
        % LP filter
        lowcutoff = [];
        highcutoff = 6;
        fprintf('Lowpass Filtering (%.1f Hz)...\n', highcutoff)
        [EEG_mocap_LP] = custom_filter(EEG_mocap_interp, lowcutoff, highcutoff);
        
        EEG_mocap_LP = pop_saveset(EEG_mocap_LP, 'filename', N.MoCapPreparedFile,'filepath', N.searchFolder_2);
        [ALLEEG, EEG_mocap_LP, CURRENTSET] = eeg_store(ALLEEG, EEG_mocap_LP, CURRENTSET);
    else
        EEG_mocap_LP = pop_loadset('filename', N.MoCapPreparedFile,'filepath', N.searchFolder_2);
        [ALLEEG, EEG_mocap_LP, CURRENTSET] = eeg_store(ALLEEG, EEG_mocap_LP, CURRENTSET);
    end
    
    % Rigid body properties
    RB1 = computeRigidBody(EEG_mocap_LP, 'RB1');
    RB2 = computeRigidBody(EEG_mocap_LP, 'RB2');
    marker_rad = 0.007; % in meters
    
    CC_events = find(strcmp({EEG_mocap_LP.event.type}, 'ConeCalibration'));
    %CC_events = find(strcmp({EEG_mocap_LP.event.type}, 'BlockStart'));
    for cc = CC_events
        % Ignore CC events in GogglesOFF blocks
        if strcmp(EEG_mocap_LP.event(cc).condition,'GogglesON')
            % Extract time around event
            latencies = ceil(EEG_mocap_LP.event(cc).latency-EEG_mocap_LP.srate/2):...
                floor(EEG_mocap_LP.event(cc).latency+EEG_mocap_LP.srate/2);
            
            % data from actual measurements
            switch subject
                case 'P1001'
                    if EEG_mocap_LP.event(cc).block == 1
                        tableLeftAngle = [-0.020,-0.055,-0.038]-marker_rad;% in meters
                    else
                        tableLeftAngle = [-0.043,-0.040,-0.041]-marker_rad;% in meters
                    end
                    tumblerWRTtable = [0.597,0.536,0];% in meters
                case 'P1001-2'
                    if EEG_mocap_LP.event(cc).trial ~= 2
                        continue
                    end
                    tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                    tumblerWRTtable = [0.521,0.436,0];% in meters
                case 'P1001-3'
                    if EEG_mocap_LP.event(cc).trial ~= 2
                        continue
                    end
                    tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                    tumblerWRTtable = [0.549,0.473,0];% in meters
                                    case 'P1002-3'
                    if EEG_mocap_LP.event(cc).trial ~= 2
                        continue
                    end
                    tableLeftAngle = [0,0,-0.042-marker_rad];% in meters
                    tumblerWRTtable = [0.58,0.559,0];% in meters
                case 'P1004'
                    tableLeftAngle = [0,0,-0.04-marker_rad];% in meters
                    tumblerWRTtable = [0.589,0.53,0];% in meters
                case 'P1009'
                    error('Enter tumbler coordinates (see Alex''s notes)');
                case 'P1001-4'
                    tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                    tumblerWRTtable = [0.555,0.426,0];% in meters
                case 'P1004-2'
                    warning('Please note that a different table was used for that recording');
                    tableLeftAngle = [0,0,-0.039-marker_rad];% in meters
                    tumblerWRTtable = [0.713,0.373,0];% in meters
            end
            
            figure;
            plotRigidBody(RB2, latencies, 'ave', true, tableLeftAngle, tumblerWRTtable);
            title(sprintf('%s - Cone Calibration - Block %d', subject, EEG_mocap_LP.event(cc).block))
            saveCurrentFig(fullfile(study_config.figures_folder, 'MocapAnimations', filesep),...
                sprintf('%s_ConeCalibration-B%d', subject, EEG_mocap_LP.event(cc).block), {'png'}, [])
        end
    end
    
    %% Compute percentage of tumbler visibility
    % Compute Rotation matrices for cone
    cone_RM = cat(2,reshape(RB2.Zvector,size(RB2.pos,1),1,size(RB2.pos,2)),...
        reshape(RB2.Xvector,size(RB2.pos,1),1,size(RB2.pos,2)),...
        reshape(RB2.Yvector,size(RB2.pos,1),1,size(RB2.pos,2)));
    
    events = EEG_mocap_LP.event;
    samp_step = 8;
    if ~isfile(fullfile(N.searchFolder_2, sprintf('%s_DataBase_withMocap.mat', subject))) || overwritePercVis
        Trials = EEG_mocap_LP.etc.TrialData;
        %perc_step = 5; edges = 0:perc_step:100;
        PercTumblerVis_EC = cell(size(Trials,1),1);
        %tumbVis_EC = nan(size(Trials,1),length(edges)-1);
        PercTumblerVis_EO = cell(size(Trials,1),1);
        %tumbVis_EO = nan(size(Trials,1),length(edges)-1);
        parfor t = 1:size(Trials,1)
            if strcmp(Trials.Condition{t}, 'GogglesON') && strcmp(Trials.TrialType{t}, 'WithObject')
                switch Trials.ID{t}
                    case 'P1001'
                        if Trials.Block(t) == 1
                            tableLeftAng = [-0.020,-0.055,-0.038]-marker_rad;% in meters
                        else
                            tableLeftAng = [-0.043,-0.040,-0.041]-marker_rad;% in meters
                        end
                        tumblerWRTtab = [0.597,0.536,0];% in meters
                    case 'P1001-2'
                        tableLeftAng = [0,0,-0.041-marker_rad];% in meters
                        tumblerWRTtab = [0.521,0.436,0];% in meters
                    case 'P1001-3'
                        tableLeftAng = [0,0,-0.041-marker_rad];% in meters
                        tumblerWRTtab = [0.549,0.473,0];% in meters
                                            case 'P1002-3'
                        tableLeftAng = [0,0,-0.042-marker_rad];% in meters
                        tumblerWRTtab = [0.58,0.559,0];% in meters
                    case 'P1004'
                        tableLeftAng = [0,0,-0.04-marker_rad];% in meters
                        tumblerWRTtab = [0.589,0.53,0];% in meters
                    case 'P1009'
                        error('Enter tumbler coordinates (see Alex''s notes)');
                    case 'P1001-4'
                        tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                        tumblerWRTtable = [0.555,0.426,0];% in meters
                    case 'P1004-2'
                        warning('Please note that a different table was used for that recording');
                        tableLeftAngle = [0,0,-0.039-marker_rad];% in meters
                        tumblerWRTtable = [0.713,0.373,0];% in meters
                end
                
                trialStart_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,1));
                obsStart_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,2));
                obsEnd_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,3));
                
                latencies_EC = ceil(events(trialStart_ev).latency):samp_step:floor(events(obsStart_ev).latency);
                PercTumblerVis_EC{t} = percTumblerInCone(RB2.pos(:,latencies_EC), cone_RM(:,:,latencies_EC),...
                    tableLeftAng, tumblerWRTtab);
                %[tumbVis_EC(t,:),~] = histcounts(tumbVolVis_EC{t}, edges);
                
                latencies_EO = ceil(events(obsStart_ev).latency):samp_step:floor(events(obsEnd_ev).latency);
                PercTumblerVis_EO{t} = percTumblerInCone(RB2.pos(:,latencies_EO), cone_RM(:,:,latencies_EO),...
                    tableLeftAng, tumblerWRTtab);
                %[tumbVis_EO(t,:),~] = histcounts(tumbVolVis_EO{t}, edges);
            end
        end
        
        Trials = addvars(Trials, PercTumblerVis_EC, PercTumblerVis_EO);
        save(fullfile(N.searchFolder_2, sprintf('%s_DataBase_withMocap', subject)), 'Trials');
    else
        load(fullfile(N.searchFolder_2, sprintf('%s_DataBase_withMocap.mat', subject))); % Load Trials variable
    end
    
    %% Plot Trials animation (time consuming)
        options = struct();
        options.record = true;
        options.animationSpeed = 1;
        options.srate = EEG_mocap_LP.srate;
        options.sampStepVisData = samp_step;
        for t = 1:size(Trials,1)
            if strcmp(Trials.Condition{t}, 'GogglesON')
                if strcmp(Trials.TrialType{t}, 'WithObject')
                    options.plot_tumbler = true;
                else
                    options.plot_tumbler = false;
                    continue
                end
    
                if isfile(fullfile(study_config.figures_folder, 'MocapAnimations',...
                        sprintf('%s_B%d-T%d_EO.mp4', Trials.ID{t},Trials.Block(t), Trials.Trial(t))))
                    continue
                end
    
                switch Trials.ID{t}
                    case'P1001'
                        if Trials.Block(t) == 1
                            options.tableLeftAngle = [-0.020,-0.055,-0.038]-marker_rad;% in meters
                        else
                            options.tableLeftAngle = [-0.043,-0.040,-0.041]-marker_rad;% in meters
                        end
                        options.tumblerWRTtable = [0.597,0.536,0];% in meters
                    case 'P1001-2'
                        options.tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                        options.tumblerWRTtable = [0.521,0.436,0];% in meters
                                            case 'P1001-3'
                        options.tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                        options.tumblerWRTtable = [0.549,0.473,0];% in meters
                                                                    case 'P1002-3'
                        options.tableLeftAngle = [0,0,-0.042-marker_rad];% in meters
                        options.tumblerWRTtable = [0.58,0.559,0];% in meters
                    case 'P1004'
                        options.tableLeftAngle = [0,0,-0.04-marker_rad];% in meters
                        options.tumblerWRTtable = [0.589,0.53,0];% in meters
                    case 'P1009'
                        error('Enter tumbler coordinates (see Alex''s notes)');
                    case 'P1001-4'
                        tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                        tumblerWRTtable = [0.555,0.426,0];% in meters
                    case 'P1004-2'
                        warning('Please note that a different table was used for that recording');
                        tableLeftAngle = [0,0,-0.039-marker_rad];% in meters
                        tumblerWRTtable = [0.713,0.373,0];% in meters
                end
    
                trialStart_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,1));
                obsStart_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,2));
                obsEnd_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,3));
    
                options.percTumbVis = Trials.PercTumblerVis_EC{t};
                options.titleFile = fullfile(study_config.figures_folder, 'MocapAnimations',...
                    sprintf('%s_B%d-T%d_EC', Trials.ID{t},Trials.Block(t), Trials.Trial(t)));
                figure;
                title({sprintf('%s - Block %d, Trial %d', Trials.ID{t},Trials.Block(t), Trials.Trial(t)),...
                    sprintf('Eyes closed - play speed x%d',options.animationSpeed)}, 'FontSize', 16);
                animateGogglesUse(RB2, ceil(events(trialStart_ev).latency):floor(events(obsStart_ev).latency), options);
    
                options.percTumbVis = Trials.PercTumblerVis_EO{t};
                options.titleFile = fullfile(study_config.figures_folder, 'MocapAnimations',...
                    sprintf('%s_B%d-T%d_EO', Trials.ID{t},Trials.Block(t), Trials.Trial(t)));
                figure;
                title({sprintf('%s - Block %d, Trial %d', Trials.ID{t},Trials.Block(t), Trials.Trial(t)),...
                    sprintf('Eyes open - play speed x%d',options.animationSpeed)}, 'FontSize', 16);
                animateGogglesUse(RB2, ceil(events(obsStart_ev).latency):floor(events(obsEnd_ev).latency), options);
            end
        end
    
    %% Add percentage of tumbler visibility to EEG processed file
    EEG_final = pop_loadset('filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej_ICcats);
    
    latencies = cell(size(Trials,1),1);
    PercTumblerVis = cell(size(Trials,1),1);
    parfor t = 1:size(Trials,1)
        if strcmp(Trials.TrialType{t}, 'WithObject')
            if strcmp(Trials.Condition{t}, 'GogglesON')
                switch Trials.ID{t}
                    case 'P1001'
                        if Trials.Block(t) == 1
                            tabLeftAng = [-0.020,-0.055,-0.038]-marker_rad;% in meters
                        else
                            tabLeftAng = [-0.043,-0.040,-0.041]-marker_rad;% in meters
                        end
                        tumbWRTtab = [0.597,0.536,0];% in meters
                    case 'P1001-2'
                        tabLeftAng = [0,0,-0.041-marker_rad];% in meters
                        tumbWRTtab = [0.521,0.436,0];% in meters
                    case 'P1001-3'
                        tabLeftAng = [0,0,-0.041-marker_rad];% in meters
                        tumbWRTtab = [0.549,0.473,0];% in meters
                        case 'P1002-3'
                        tabLeftAng = [0,0,-0.042-marker_rad];% in meters
                        tumbWRTtab = [0.58,0.559,0];% in meters
                    case 'P1004'
                        tabLeftAng = [0,0,-0.04-marker_rad];% in meters
                        tumbWRTtab = [0.589,0.53,0];% in meters
                    case 'P1009'
                        error('Enter tumbler coordinates (see Alex''s notes)');
                    case 'P1001-4'
                        tableLeftAngle = [0,0,-0.041-marker_rad];% in meters
                        tumblerWRTtable = [0.555,0.426,0];% in meters
                    case 'P1004-2'
                        warning('Please note that a different table was used for that recording');
                        tableLeftAngle = [0,0,-0.039-marker_rad];% in meters
                        tumblerWRTtable = [0.713,0.373,0];% in meters
                end
                
                trialStart_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,1));
                trialEnd_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,end));
                
                latencies{t} = ceil(events(trialStart_ev).latency):samp_step:floor(events(trialEnd_ev).latency);
                PercTumblerVis{t} = percTumblerInCone(RB2.pos(:,latencies{t}), cone_RM(:,:,latencies{t}),...
                    tabLeftAng, tumbWRTtab);
            end
        else
            % No object
            trialStart_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,1));
            trialEnd_ev = findInStructWithEmpties(events, 'urevent', Trials.urevent_seq(t,end));
            latencies{t} = ceil(events(trialStart_ev).latency):floor(events(trialEnd_ev).latency);
            PercTumblerVis{t} = zeros(1,length(latencies{t}));
        end
    end
    
    EEG_final_wMC = EEG_final;
    EEG_final_wMC.nbchan = EEG_final.nbchan+1;
    EEG_final_wMC.data = cat(1, EEG_final_wMC.data, nan(1,EEG_final.pnts));
    EEG_final_wMC.chanlocs(EEG_final_wMC.nbchan).type = 'MOCAP_comp';
    EEG_final_wMC.chanlocs(EEG_final_wMC.nbchan).unit = 'perc';
    EEG_final_wMC.chanlocs(EEG_final_wMC.nbchan).labels = 'TumblerVisibility';
    
    %     EEG_final_wMC.chanlocs(EEG_final_wMC.nbchan) = struct('type', 'MOCAP_comp', 'unit', 'perc', 'labels', 'TumblerVisibility',...
    %         'X', [], 'Y', [], 'Z', [], 'sph_phi', [], 'sph_radius', [], 'theta', [], 'radius', [], 'sph_theta', [], 'ref', '', 'urchan', []);
    
    for t = 1:size(Trials,1)
        if strcmp(Trials.TrialType{t}, 'WithObject')
            if strcmp(Trials.Condition{t}, 'GogglesON')
                all_latencies = latencies{t}(1):latencies{t}(end);
                EEG_final_wMC.data(end,all_latencies) = interp1(latencies{t},PercTumblerVis{t},all_latencies, 'pchip');
            end
        else
            % No object
            EEG_final_wMC.data(end,latencies{t}) = PercTumblerVis{t};
        end
    end
    
    pop_saveset(EEG_final_wMC, 'filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej_ICcats);
end