
DataBase = [];
clear all;
config_PIONEER_Tumbler;

%% Compile some general data for further analyses
% Options for this script
recompute = true;
%mocap_needed = false; %boolean to load MOCAP data or not
switch study_config.badSampsRejection
    case 'app'
        badSampsMethod = 'APP';
    case 'asr'
        badSampsMethod = 'ASR';
    case 'autoMoBI'
        badSampsMethod = 'autoMoBI';
end

%% Get the folders names right:
pipe_name = study_config.globalArchitecture;
fig_path = fullfile(study_config.figures_folder, 'Epochs');
if ~exist(fig_path, 'dir')
    mkdir(fig_path)
end

N = makeFolderFileNames(study_config, study_config.subjects(study_config.current_subject).id);
fname = 'MergedDataBase_forEEG.mat';
file2save = fullfile(N.searchFolder_2, fname);
if exist(file2save, 'file') && ~exist('MergedDataBase', 'var') && ~recompute
    load(file2save)
end

for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    
    % clear RAM
    %STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    %subject_ind = 41;
    subject = study_config.subjects(subject_ind).id;
    study_config.current_subject = subject_ind;
    
    % skip the first subject
    if strcmp(study_config.subjects(subject_ind).excluded, 'Yes')
        continue
    elseif exist('MergedDataBase', 'var') && ~recompute && sum(strcmp({MergedDataBase.ID}, subject))>0
        continue
    else
        N = makeFolderFileNames(study_config, subject);
        filepath = N.searchFolder_2;
        
        switch task
            case 'Tumbler'
                % load DataBase
                load(fullfile(filepath, sprintf('%s_DataBase.mat', subject)));
                
                % load prepared Dataset
                EEG_prepared = pop_loadset('filename', N.preparedFile, 'filepath', N.searchFolder_2, 'loadmode', 'info');
                %[ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
                events = EEG_prepared.event;
                %urevents = EEG_prepared.urevent;
                %         ObsStarts = events(strcmp({events.type},'ObsStart'));
                %         ObsEnds = events(strcmp({events.type},'ObsEnd'));
                %         BaseStarts = events(strcmp({events.type},'BaseStart'));
                %         BaseEnds = events(strcmp({events.type},'BaseEnd'));
                %Bounds = events(strcmp({events.type},'boundary'));
                noBounds = events(~strcmp({events.type},'boundary'));
                
                % load Preprocessed dataset
                EEG_preproc = pop_loadset('filename', N.preICAFile, 'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
                badSamples = EEG_preproc.etc.(badSampsMethod).rejectedSamples;
                clear EEG_preproc
                
                fields = {'ID', 'Gender', 'Block', 'Condition', 'Trial', 'TrialType', 'Answer', 'Delay_ms', 'ConfidenceLevel',...
                    'urevent_seq', 'duration_ms_seq', 'completeEEG_full', 'completeEEG_details', 'cleanEEG_full', 'cleanEEG_details'};
                if ~exist('MergedDataBase', 'var')
                    MergedDataBase = cell2struct(cell(numel(fields),1), fields', 1);
                    lastLine = 0;
                else
                    lastLine = size(MergedDataBase,2);
                end
                N_newTrials = size(DataBase,1);
                
                for i = 1:N_newTrials
                    blk = DataBase.Block(i);
                    trl = DataBase.Trial(i);
                    
                    %             if strcmp(subject,'P1004') && blk==1 && (trl==6 || trl==7)
                    %                 disp('stop')
                    %             end
                    if any(strcmp(fieldnames(DataBase), 'ConfidenceLevel'))
                        cfd_lvl = DataBase.ConfidenceLevel(i);
                    else
                        cfd_lvl = NaN;
                    end
                    
                    n_events = length(DataBase.urevent_seq(i,:));
                    dur_ms_seq = nan(1,n_events-1);
                    n_phases = length(DataBase.EEGComplete_details(i,:));
                    percCleanEEG_det = nan(1,n_phases-1);
                    
                    % Check if this trial was kept in the prepared dataset
                    % (rejected if not EEG complete or the buffer condition was not satisfied)
                    if (any([noBounds.block] == blk & [noBounds.trial] == trl)) && ...
                        ~(strcmp(subject,'P1001-4') && blk==1 && (trl == 6 || trl == 8)) && ... % adding exception for missing eeg data
                        ~(strcmp(subject,'P1001-4') && blk==7 && trl == 2) % adding exception for missing eeg data
                        % Compute total EEG cleanliness
                        start_lat = floor(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,1)).latency);
                        stop_lat = ceil(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,end)).latency);
                        dur_samps = 1+stop_lat-start_lat;
                        percCleanEEG = 100*(sum(~badSamples(start_lat:stop_lat))/dur_samps);
                        
                        for e = 2:n_events
                            start_lat = floor(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,e-1)).latency);
                            stop_lat = ceil(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,e)).latency);
                            dur_samps = 1+stop_lat-start_lat;
                            
                            % Compute durations
                            dur_ms_seq(e-1) = 1000*dur_samps/EEG_prepared.srate;
                        end
                        
                        for ph = 1:n_phases
                            switch ph
                                case 1
                                    % Eyes closed phase
                                    start_ev = 1;
                                    stop_ev = 2;
                                case 2
                                    % Eyes open phase
                                    start_ev = 2;
                                    stop_ev = 3;
                                case 3
                                    % Answer phase
                                    start_ev = 3;
                                    stop_ev = 6;
                            end
                            start_lat = floor(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,start_ev)).latency);
                            stop_lat = ceil(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,stop_ev)).latency);
                            dur_samps = 1+stop_lat-start_lat;
                            % Compute EEG cleanliness section by section
                            if DataBase.EEGComplete_details(i,ph) == 100
                                percCleanEEG_det(ph) = 100*(sum(~badSamples(start_lat:stop_lat))/dur_samps);
                            end
                        end
                    else
                        % If incomplete, this trial was already rejected in the prepared dataset
                        percCleanEEG = NaN;
                    end
                    
                    values = {subject, DataBase.Gender(i), blk, DataBase.Condition(i), trl, DataBase.TrialType(i),...
                        DataBase.Answer(i), DataBase.Delay_ms(i), cfd_lvl, DataBase.urevent_seq(i,:),...
                        dur_ms_seq, DataBase.EEGComplete_trial(i), DataBase.EEGComplete_details(i,:),...
                        percCleanEEG, percCleanEEG_det};
                    
                    MergedDataBase(lastLine + i) = cell2struct(values', fields', 1);
                end
                
                % Get rid of large structures:
                clear events noBounds
                clear EEG_prepared
            case 'Tumbler2020'
                EEG_cleaned = pop_loadset('filename', N.postLabelingFile, 'filepath', N.searchFolder_2arch_rej_ICcats, 'loadmode', 'info');
                events = EEG_cleaned.event;
                noBounds = events(~strcmp({events.type},'boundary'));
                fields = {'ID', 'Gender', 'Block', 'Condition', 'Trial', 'TrialType', 'Answer', 'Delay_ms', 'ConfidenceLevel',...
                    'urevent_seq', 'duration_ms_seq', 'completeEEG_full'};
                if ~exist('MergedDataBase', 'var')
                    MergedDataBase = cell2struct(cell(numel(fields),1), fields', 1);
                    lastLine = 0;
                else
                    lastLine = size(MergedDataBase,2);
                end
                trialInds = unique([noBounds.trialIndex]);
                trialInds = trialInds(~isnan(trialInds));
                N_newTrials = length(trialInds);
                
                for i = 1:N_newTrials
                    trEvents = find([noBounds.trialIndex] == trialInds(i));
                    
                    blk = ceil(trialInds(i)/10);
                    switch noBounds(trEvents(1)).condition
                        case 'Monoc-Glass-NoCR'
                            cnd = 'GogglesON';
                        case 'Monoc-NoGlass-NoCR'
                            cnd = 'GogglesOFF';
                        case 'Binoc-NoGlass-NoCR'
                            cnd = 'GogglesOFF-Binoc';
                    end
                    trl = mod(trialInds(i)-1,10)+1;
                    switch noBounds(trEvents(1)).trialType
                        case 'NoObject'
                            trltyp = 'WithoutObject';
                        case 'Object'
                            trltyp = 'WithObject';
                    end
                    
                    if any(contains({noBounds(trEvents).reason},'Object'))
                        aswEvent = find(contains({noBounds(trEvents).reason},'Object'));                        
                        if strcmp(noBounds(trEvents(aswEvent)).reason,'ObjectPresent')
                            asw = 'Present';
                        elseif strcmp(noBounds(trEvents(aswEvent)).reason,'ObjectAbsent')
                            asw = 'Absent';
                        end
                        delay_ms = 1000*(noBounds(trEvents(aswEvent)).latency - noBounds(trEvents(aswEvent-1)).latency)/EEG_cleaned.srate;
                    else
                        aswEvent = NaN;
                        asw = 'none';
                        delay_ms = NaN;
                    end
                    cfd_lvl = NaN;
                    
                    n_events = length(trEvents);       
                                        if isnan(aswEvent)
                        seq = [noBounds(trEvents).urevent];
                    else
                        seq = [noBounds(trEvents(setdiff(1:n_events,aswEvent))).urevent];
                    end
                    
                    dur_ms_seq = nan(1,length(seq)-1);
                    for e = 2:length(seq)
                        dur_ms_seq(e-1) = 1000*(noBounds([noBounds.urevent] == seq(e)).latency -...
                            noBounds([noBounds.urevent] == seq(e-1)).latency)/EEG_cleaned.srate;
                    end                    
                    
                    values = {subject, study_config.subjects(study_config.current_subject).gender,...
                        blk, cnd, trl, trltyp, asw, delay_ms, cfd_lvl, seq, dur_ms_seq,100};
                    
                    MergedDataBase(lastLine + i) = cell2struct(values', fields', 1);
                end
        end
    end
end

MergedDataBase = struct2table(MergedDataBase);

%% save
if ~exist(file2save) || recompute
    save(file2save,'MergedDataBase')
end

% Description of the data
figure
subplot(1,3,1)
[~, edges] = histcounts(0.001*MergedDataBase.duration_ms_seq(:,1));
hold on
h1 = histogram(0.001*MergedDataBase.duration_ms_seq(strcmp(MergedDataBase.ID, 'P1001'),1),edges);
h4 = histogram(0.001*MergedDataBase.duration_ms_seq(strcmp(MergedDataBase.ID, 'P1004'),1),edges);
xlabel('Duration (s)')
ylabel('Trial count')
title('Eyes Closed phase')
legend([h1,h4], {'P1001', 'P1004'}, 'Location', 'Best')

subplot(1,3,2)
[~, edges] = histcounts(0.001*MergedDataBase.duration_ms_seq(:,2));
hold on
h1 = histogram(0.001*MergedDataBase.duration_ms_seq(strcmp(MergedDataBase.ID, 'P1001'),2),edges);
h4 = histogram(0.001*MergedDataBase.duration_ms_seq(strcmp(MergedDataBase.ID, 'P1004'),2),edges);
xlabel('Duration (s)')
title('Eyes Open phase')
legend([h1,h4], {'P1001', 'P1004'}, 'Location', 'Best')

subplot(1,3,3)
[~, edges] = histcounts(0.001*sum(MergedDataBase.duration_ms_seq(:,3:end),2));
hold on
h1 = histogram(0.001*sum(MergedDataBase.duration_ms_seq(strcmp(MergedDataBase.ID, 'P1001'),3:end),2),edges);
h4 = histogram(0.001*sum(MergedDataBase.duration_ms_seq(strcmp(MergedDataBase.ID, 'P1004'),3:end),2),edges);
xlabel('Duration (s)')
title('Answer phase')
legend([h1,h4], {'P1001', 'P1004'}, 'Location', 'Best')
