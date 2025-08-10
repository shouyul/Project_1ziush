% xdf2set_PIONEER - Import function from raw .xdf to .set EEGLAB file format
% Incorporate steps specific to PIONEER experiment (mainly events processing)
% Data is being loaded with loadxdf_AD() and streams are being assigned accordingly to the config
% Optionally, if the .xdf dataset has been listed in the missingData field (in study_config.subjects),
% a recovery step is proposed to the user (interactive, see completeEEGfromEEGO function).
% Then, stream by stream all samples are assigned the closest time stamp in
% the common timeline (in case of irregular frame rates badly handled by load_xdf). Will take quite some time
% Events are preprocessed in this function (see export_events_Tumbler and read_trialtypes_Tumbler functions)
% Eventually, data is saved in the .set format, for all streams and by stream types
% Once all datasets for the subject have been processed, they are merged into a single file (AllStreams and EEG types only)
%
% The intermediate files are stored on the disk.
%
% Usage:
%   >>  [ALLEEG, EEG_merged, CURRENTSET] = xdf2set_PIONEER(ALLEEG, CURRENTSET, subject_ind, study_config, overwrite)
%
% Inputs:
%   subject_ind     - subject number of the current subject relative to study_config.subjects list (necessary for filepaths and storage)
%   study_config    - configuration struct with all necessary information.
%   ALLEEG          - complete EEGLAB data set structure
%   CURRENTSET      - index of current EEGLAB EEG structure within ALLEEG
%	overwrite		- boolean to force import even if the file output already exist
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_merged				  - merged EEGLAB EEG structure, contains EEG data (not AllStreams) of all datasets.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of EEGLAB EEG structures are stored on disk according to their names in the bemobil_config
%
% Authors: Alexandre Delaux, 2023

function [ALLEEG, EEG_merged, CURRENTSET] = xdf2set_PIONEER(ALLEEG, CURRENTSET, subject_ind, study_config, overwrite)

subject = study_config.subjects(subject_ind).id;
%disp(['Subject ' num2str(subject)]);

input_filepath = [study_config.study_folder study_config.raw_data_folder subject];
output_filepath = [study_config.study_folder study_config.raw_EEGLAB_data_folder subject];

file_prefix{1} = [subject '_AllStreams_'];

for i_filename = 1:length(study_config.filenames)
    %i_filename = 4;
    full_filename = [subject, '_', study_config.filenames{i_filename}];
    
    %continue
    
    %% import
    input_file = fullfile(input_filepath, [full_filename '.xdf']);
    if isfile(input_file)
        if isfile([output_filepath, filesep, file_prefix{1}, study_config.filenames{i_filename}, '.set']) && ~overwrite
            disp([file_prefix{1}, full_filename, '.set already exists : not importing XDF file']);
            continue
        else
            disp(['Importing file: "' full_filename '.xdf" ...']);
            [AllStreams, ~] = load_xdf_AD(input_file);
            %[AllStreams, ~] = load_xdf_AD(input_file, 'Verbose', true);
            %[AllStreams, ~] = load_xdf_AD(input_file, 'Verbose', true, 'HandleClockResets', false);
        end
    else
        continue
    end
    
    % Detect Stream types
    EEG_stream_inds = zeros(1,study_config.stream_count(1));
    ET_stream_inds = zeros(1,study_config.stream_count(2));
    MOCAP_stream_inds = zeros(1,study_config.stream_count(3));
    Event_stream_inds = zeros(1,study_config.stream_count(4));
    for s = 1:numel(AllStreams)
        stream_type = AllStreams{s}.info.type;
        if (sum(contains(study_config.eeg_streams,stream_type))>0) || strcmp(stream_type, 'impedance')
            ind2replace = find(contains(study_config.eeg_streams_order, AllStreams{s}.info.name));
            if EEG_stream_inds(ind2replace) == 0
                EEG_stream_inds(ind2replace)=s;
            else
                % Catch if the same stream is present more than once
                fprintf('Old stream length: %d sample(s).\n', size(AllStreams{EEG_stream_inds(ind2replace)}.time_series,2))
                fprintf('New stream length: %d sample(s).\n', size(AllStreams{s}.time_series,2))
                kept = input('Which one do I keep? [0=old, 1=new]: '); % enter 0 or 1
                if kept
                    EEG_stream_inds(ind2replace)=s;
                end
            end
        elseif (sum(contains(study_config.eye_tracker_streams,stream_type))>0)
            curr_ind = find(ET_stream_inds,1,'last');
            if isempty(curr_ind)
                ET_stream_inds(1)=s;
            else
                ET_stream_inds(curr_ind+1)=s;
            end
        elseif (sum(contains(study_config.rigidbody_streams,stream_type))>0)
            ind2replace = find(contains(study_config.rb_streams_order, AllStreams{s}.info.name));
            if MOCAP_stream_inds(ind2replace) == 0
                MOCAP_stream_inds(ind2replace)=s;
            else
                % Catch if the same stream is present more than once
                fprintf('Old stream length: %d sample(s).\n', size(AllStreams{MOCAP_stream_inds(ind2replace)}.time_series,2))
                fprintf('New stream length: %d sample(s).\n', size(AllStreams{s}.time_series,2))
                kept = input('Which one do I keep? [0=old, 1=new]: '); % enter 0 or 1
                if kept
                    MOCAP_stream_inds(ind2replace)=s;
                end
            end
        elseif (sum(contains(study_config.event_streams,stream_type))>0)
            if sum(contains(study_config.eeg_streams_order,AllStreams{s}.info.name))==0
                % To avoid loading EEG Marker streams (useless in ANT)
                curr_ind = find(Event_stream_inds,1,'last');
                if isempty(curr_ind)
                    Event_stream_inds(1)=s;
                else
                    Event_stream_inds(curr_ind+1)=s;
                end
            end
        else
            error(['Could not recognize the stream type, please make sure all stream types fall into',...
                ' a category defined in the config file']);
        end
    end
    
    %% Complete Missing data from eego file (if necessary & possible)
    % Check this function works well for PIONEER experiment before using it
    if sum(strcmp(study_config.filenames{i_filename}, study_config.subjects(subject_ind).missingData))>0
        AllStreams = completeEEGfromEEGO_PIONEER(AllStreams, EEG_stream_inds, study_config.filenames{i_filename},...
            study_config.subjects, input_filepath, study_config.figures_folder);
    end
    
    %% Read Streams
    Data_all = eeg_emptyset;
    
    %% Get info necessary to create the data:
    %% EEG:
    n_chans_eeg = 0;
    t_min = +inf;
    t_max = -inf;
    srate_eeg = -1;
    segments2keep_eeg = cell(length(EEG_stream_inds),1);
    for ee = 1:length(EEG_stream_inds)
        EEG_stream = AllStreams{EEG_stream_inds(ee)};
        n_chans_eeg = n_chans_eeg + size(EEG_stream.time_series,1);
        
        srate_eeg = str2num(EEG_stream.info.nominal_srate);
        
        % Look at Segments to understand if there is weird data
        for seg=1:numel(EEG_stream.segments)
            if EEG_stream.segments(seg).duration == 0
                % Invalid segment
            else
                % Valid segment
                last_valid_seg = segments2keep_eeg{ee}(find(segments2keep_eeg{ee},1,'last'));
                if ~isempty(last_valid_seg)
                    % A previous valid segment has been found
                    if EEG_stream.segments(seg).t_begin < EEG_stream.segments(last_valid_seg).t_end
                        % An overlap exists --> Remove the overlap period
                        if EEG_stream.segments(seg).t_begin <= EEG_stream.segments(last_valid_seg).t_begin
                            % Remove the whole last valid seg
                            segments2keep_eeg{ee} = [segments2keep_eeg{ee}(1:end-1), seg];
                            warning('Overlap: removed a full segment of %.2f seconds from EEG data',...
                                EEG_stream.segments(last_valid_seg).duration);
                            
                            t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                            t_max = max(t_max,EEG_stream.segments(seg).t_end);
                        elseif EEG_stream.segments(seg).t_end <= EEG_stream.segments(last_valid_seg).t_end
                            % Do not add the current seg
                            warning('Overlap: removed a full segment of %.2f seconds from EEG data',...
                                EEG_stream.segments(seg).duration);
                        else
                            t_discard_begin = EEG_stream.segments(seg).t_begin;
                            t_discard_stop = EEG_stream.segments(last_valid_seg).t_end;
                            
                            % Modify last_valid_seg:
                            EEG_stream.segments(last_valid_seg).t_end = t_discard_begin;
                            EEG_stream.segments(last_valid_seg).duration = t_discard_begin - EEG_stream.segments(last_valid_seg).t_begin;
                            EEG_stream.segments(last_valid_seg).index_range(2) = ...
                                find(EEG_stream.time_stamps(EEG_stream.segments(last_valid_seg).index_range(1):...
                                EEG_stream.segments(last_valid_seg).index_range(2)) <= t_discard_begin, 1, 'last')...
                                + EEG_stream.segments(last_valid_seg).index_range(1) - 1;
                            EEG_stream.segments(last_valid_seg).num_samples = diff(EEG_stream.segments(last_valid_seg).index_range)+1;
                            EEG_stream.segments(last_valid_seg).effective_srate = (EEG_stream.segments(last_valid_seg).num_samples-1)...
                                /EEG_stream.segments(last_valid_seg).duration;
                            
                            % Modify seg:
                            EEG_stream.segments(seg).t_begin = t_discard_stop;
                            EEG_stream.segments(seg).duration = EEG_stream.segments(seg).t_end - t_discard_stop;
                            EEG_stream.segments(seg).index_range(1) = ...
                                find(EEG_stream.time_stamps(EEG_stream.segments(seg).index_range(1):...
                                EEG_stream.segments(seg).index_range(2)) >= t_discard_stop, 1)...
                                + EEG_stream.segments(seg).index_range(1) - 1;
                            EEG_stream.segments(seg).num_samples = diff(EEG_stream.segments(seg).index_range)+1;
                            EEG_stream.segments(seg).effective_srate = (EEG_stream.segments(seg).num_samples-1)...
                                /EEG_stream.segments(seg).duration;
                            
                            segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                            warning('Overlap: removed a partial segment of %.2f seconds from EEG data',...
                                t_discard_stop - t_discard_begin);
                            
                            t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                            t_max = max(t_max,EEG_stream.segments(seg).t_end);
                        end
                    else
                        segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                        t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                        t_max = max(t_max,EEG_stream.segments(seg).t_end);
                    end
                else
                    segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                    t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                    t_max = max(t_max,EEG_stream.segments(seg).t_end);
                end
            end
        end
        
        AllStreams{EEG_stream_inds(ee)} = EEG_stream;
    end
    
    % Remove channels if indicated:
    for ee = 1:length(EEG_stream_inds)
        if ~isempty(study_config.channels_to_remove{ee})
            n_chans_eeg = n_chans_eeg - length(study_config.channels_to_remove{ee});
        end
    end
    
    %% ET:
    if ~isempty(ET_stream_inds)
        if length(ET_stream_inds) == 1
            % Only one stream
            t_min = min(t_min, AllStreams{ET_stream_inds(1)}.time_stamps(1));
            t_max = max(t_max, AllStreams{ET_stream_inds(1)}.time_stamps(end));
            
            % Two channels are redundant and contain time data
            n_chans_et = str2num(AllStreams{ET_stream_inds(1)}.info.channel_count) - 2;
            
            % Check that nominal srate and effective srate are not too far appart:
            nom_srate = str2num(AllStreams{ET_stream_inds(1)}.info.nominal_srate);
            eff_srate = AllStreams{ET_stream_inds(1)}.info.effective_srate;
            if (eff_srate < nom_srate*0.95 || eff_srate > nom_srate*1.05)
                warning('Eye tracking data: Effective srate (%.2f Hz) far from nominal srate (%.2f Hz)', eff_srate, nom_srate)
            end
            %srate_et = nom_srate;
            srate_et = round(eff_srate);
        else
            error('xdf2set not ready for dealing with 2 eye tracking streams')
        end
    else
        n_chans_et = 0;
    end
    
    %% MOCAP:
    n_chans_mocap = 0;
    srate_mocap = -1;
    if ~isempty(MOCAP_stream_inds)
        for moc = 1:length(MOCAP_stream_inds)
            if MOCAP_stream_inds(moc) > 0
                n_chans_mocap = n_chans_mocap + str2num(AllStreams{MOCAP_stream_inds(moc)}.info.channel_count);
                mocap_channels_description = AllStreams{MOCAP_stream_inds(moc)}.info.desc.channels;
                % First check that nominal srate and effective srate are not too far appart:
                nom_srate = str2num(AllStreams{MOCAP_stream_inds(moc)}.info.nominal_srate);
                if nom_srate ~= 0
                    eff_srate = AllStreams{MOCAP_stream_inds(moc)}.info.effective_srate;
                    if isnan(eff_srate)
                        warning('Error with MOCAP stream %d', moc)
                        % Problem with recording, missing data for this
                        % block (only 1 sample)
                        MOCAP_stream_inds(moc) = 0;
                    else
                        t_min = min(t_min, AllStreams{MOCAP_stream_inds(moc)}.time_stamps(1));
                        t_max = max(t_max, AllStreams{MOCAP_stream_inds(moc)}.time_stamps(end));
                        
                        if (eff_srate < nom_srate*0.95 || eff_srate > nom_srate*1.05)
                            error('Effective srate too far from nominal srate')
                        end
                    end
                    
                    % Then check that nominal srate are the same for all Mocap recordings
                    if (srate_mocap==-1)
                        srate_mocap = nom_srate;
                    elseif (srate_mocap ~= nom_srate)
                        error('FATAL ERROR: Not the same srate between Mocap recordings')
                    end
                else
                    srate_mocap = 0;
                    t_min = min(t_min, AllStreams{MOCAP_stream_inds(moc)}.time_stamps(1));
                    t_max = max(t_max, AllStreams{MOCAP_stream_inds(moc)}.time_stamps(end));
                end
            elseif MOCAP_stream_inds(moc) == 0 && moc == 2
                % In some cases, it is normal that the second mocap stream
                % is missing (goggles OFF condition)
                n_chans_mocap = n_chans_mocap*2;
            else
                error('Unexpected error');
            end
        end
    end
    
    interval = 1/srate_eeg;
    times = round(t_min,4):interval:round(t_max,4)+interval;
    
    for typ = 1:length(study_config.stream_count)
        switch typ
            case 1 %EEG
                disp('Importing EEG...')
                %% Data
                data_eeg = nan(n_chans_eeg,length(times));
                for ee = 1:length(EEG_stream_inds)
                    EEG_stream = AllStreams{EEG_stream_inds(ee)};
                    
                    % To compute channels span:
                    local_data = EEG_stream.time_series;
                    if ~isempty(study_config.channels_to_remove{ee})
                        chans2keep = setdiff(1:size(local_data,1),study_config.channels_to_remove{ee});
                        local_data = local_data(chans2keep,:);
                    end
                    
                    if ee == 1
                        channels_span = 1:size(local_data,1);
                    else
                        next_chan = channels_span(end)+1;
                        channels_span = next_chan:next_chan+size(local_data,1)-1;
                    end
                    
                    for seg = segments2keep_eeg{ee}
                        eff_srate = EEG_stream.segments(seg).effective_srate;
                        time_span = EEG_stream.segments(seg).index_range(1):EEG_stream.segments(seg).index_range(2);
                        local_data = EEG_stream.time_series(chans2keep, time_span);
                        local_times = EEG_stream.time_stamps(time_span);
                        
                        local2global = zeros(1,length(local_times));
                        ppm = ParforProgressbar(length(local_times), 'showWorkerProgress', true,...
                            'title', sprintf('Finding time samples for EEG%d;Seg%d', ee, seg));
                        parfor t = 1:length(local_times)
                            ind = find(times>=local_times(t),1);
                            if isempty(ind)
                                if t == length(local_times) && local_times(t)<times(end)+interval/2
                                    local2global(t) = length(times);
                                else
                                    error('Seg %d: Sample %d could not be placed',seg,t);
                                end
                            elseif ind==1
                                if t == 1 && local_times(t)>times(1)-interval/2
                                    local2global(t) = 1;
                                else
                                    error('Seg %d: Sample %d could not be placed',seg,t);
                                end
                            else
                                if round(local_times(t)-times(ind-1),4)<=interval/2
                                    local2global(t) = ind-1;
                                elseif round(times(ind)-local_times(t),4)<=interval/2
                                    local2global(t) = ind;
                                else
                                    error('Seg %d: Sample %d could not be placed',seg,t);
                                end
                            end
                            
                            ppm.increment();
                        end
                        delete(ppm);
                        
                        if (length(local2global) ~= length(unique(local2global)))
                            rep_list = find(local2global == [0,local2global(1:end-1)]);
                            for r = length(rep_list):-1:1
                                if rep_list(r) > 2 && (eff_srate >= srate_eeg*0.95 && eff_srate <= srate_eeg*1.05)...
                                        && local2global(rep_list(r)-2)==local2global(rep_list(r))-2
                                    % Correct local2global knowing samples
                                    % should be regularly spaced given the
                                    % srate
                                    local2global(rep_list(r)-1)=local2global(rep_list(r))-1;
                                else
                                    % Remove the duplicate sample and
                                    % average the values in the data
                                    local2global = [local2global(1:rep_list(r)-1),local2global(rep_list(r)+1:end)];
                                    local_times = [local_times(1:rep_list(r)-1),local_times(rep_list(r)+1:end)];
                                    local_data = [local_data(:,1:rep_list(r)-2),...
                                        mean(local_data(:,rep_list(r)-1:rep_list(r)),2),...
                                        local_data(:,rep_list(r)+1:end)];
                                end
                            end
                            
                            if (length(local2global) ~= length(unique(local2global)))
                                error('Still repetitions to deal with')
                            end
                        end
                        
                        % Check that eeg srate and effective srate are not too far appart:
                        srate_scale = round(srate_eeg/eff_srate);
                        if (eff_srate*srate_scale < srate_eeg*0.95 || eff_srate*srate_scale > srate_eeg*1.05)
                            warning('Impossible to combine the segment srate with global srate in a satisfying way')
                            
                            % portion to visualize what's happening
                            instant_srate = 1./diff(EEG_stream.time_stamps(EEG_stream.segments(seg).index_range(1):...
                                EEG_stream.segments(seg).index_range(2)));
                            histogram(instant_srate);
                            saveCurrentFig([output_filepath filesep],...
                                sprintf('%s_SrateHistogram_%s_EEG%dSeg%d',subject,full_filename, ee, seg), {'fig'},[]);
                            
                            data_eeg(channels_span,local2global)=local_data;
                        else
                            % We can consider the srate are multiples
                            fprintf('EEG%d;Seg%d: nomSR= %d, effSR= %.2f\n', ee, seg, srate_eeg, eff_srate)
                            for j = 1:srate_scale
                                data_eeg(channels_span,local2global(1:end-j+1)+j-1)=local_data(:,1:end-j+1);
                            end
                        end
                    end
                end
                
                %% chanlocs
                chanlocs_eeg = struct();
                if ~isempty(study_config.channel_locations_filename) &&...
                        isfile(fullfile(input_filepath, study_config.channel_locations_filename))
                    chanlocs_eeg = readlocs(fullfile(input_filepath, study_config.channel_locations_filename));
                else
                    error('xdf2set: Missing chanlocs file.')
                    % TODO
                    % Read from stream.info.desc.channels
                    % Get inspired from eeg_load_xdf
                end
                
                fields = fieldnames(chanlocs_eeg);
                command = '';
                for f = 1:(numel(fields)+2)
                    if f==1
                        command = [command,'''type'',','cell(1, n_chans_eeg+n_chans_et+n_chans_mocap),'];
                    elseif f==2
                        command = [command,'''unit'',[],'];
                    else
                        command = [command,'''',fields{f-2},''',[],'];
                    end
                end
                command = command(1:end-1);
                eval(['chanlocs = struct(',command,');']);
                for ch = 1:n_chans_eeg
                    if ~isempty(study_config.eog_channels) && sum(strcmp(study_config.eog_channels, chanlocs_eeg(ch).labels),1)>0
                        chanlocs(ch).type = 'EOG';
                    elseif ~isempty(study_config.bip_channels) && sum(strcmp(study_config.bip_channels, chanlocs_eeg(ch).labels),1)>0
                        chanlocs(ch).type = 'BIP';
                    else
                        chanlocs(ch).type = 'EEG';
                    end
                    chanlocs(ch).unit = 'microVolt';
                    for f = 1:numel(fields)
                        chanlocs(ch).(fields{f}) = chanlocs_eeg(ch).(fields{f});
                    end
                end
                
                %% etc
                for ee = 1:length(EEG_stream_inds)
                    Data_all.etc.(['descEEG', num2str(ee)]) = AllStreams{EEG_stream_inds(ee)}.info.desc;
                    Data_all.etc.(['infoEEG', num2str(ee)]) = rmfield(AllStreams{EEG_stream_inds(ee)}.info,'desc');
                    Data_all.etc.(['infoEEG', num2str(ee)]).('segments') = AllStreams{EEG_stream_inds(ee)}.segments;
                    Data_all.etc.(['infoEEG', num2str(ee)]).('keptSegments') = segments2keep_eeg{ee};
                end
                disp ('...done')
                
            case 2 %ET
                if ~isempty(ET_stream_inds)
                    % Only one stream
                    disp('Importing ET...')
                    %% Data
                    data_et = nan(n_chans_et,length(times));
                    local_data = AllStreams{ET_stream_inds(1)}.time_series(1:n_chans_et,:);
                    local_times = AllStreams{ET_stream_inds(1)}.time_stamps;
                    eff_srate = AllStreams{ET_stream_inds(1)}.info.effective_srate;
                    
                    local2global = zeros(1,length(local_times));
                    ppm = ParforProgressbar(length(local_times), 'showWorkerProgress', true,...
                        'title', 'Finding time samples for ET');
                    parfor t = 1:length(local_times)
                        ind = find(times>=local_times(t),1);
                        if isempty(ind)
                            if t == length(local_times) && local_times(t)<times(end)+interval/2
                                local2global(t) = length(times);
                            else
                                error('Sample %d could not be placed',t);
                            end
                        elseif ind==1
                            if t == 1 && local_times(t)>times(1)-interval/2
                                local2global(t) = 1;
                            else
                                error('Sample %d could not be placed',t);
                            end
                        else
                            if round(local_times(t)-times(ind-1),4)<=interval/2
                                local2global(t) = ind-1;
                            elseif round(times(ind)-local_times(t),4)<=interval/2
                                local2global(t) = ind;
                            else
                                error('Sample %d could not be placed',t);
                            end
                        end
                        
                        ppm.increment();
                    end
                    delete(ppm);
                    
                    if (length(local2global) ~= length(unique(local2global)))
                        rep_list = find(local2global == [0,local2global(1:end-1)]);
                        for r = length(rep_list):-1:1
                            if rep_list(r) > 2 && (eff_srate >= srate_et*0.95 && eff_srate <= srate_et*1.05)...
                                    && local2global(rep_list(r)-2)==local2global(rep_list(r))-2
                                % Correct local2global knowing samples
                                % should be regularly spaced given the
                                % srate
                                local2global(rep_list(r)-1)=local2global(rep_list(r))-1;
                            else
                                % Remove the duplicate sample and
                                % average the values in the data
                                local2global = [local2global(1:rep_list(r)-1),local2global(rep_list(r)+1:end)];
                                local_times = [local_times(1:rep_list(r)-1),local_times(rep_list(r)+1:end)];
                                local_data = [local_data(:,1:rep_list(r)-2),...
                                    mean(local_data(:,rep_list(r)-1:rep_list(r)),2),...
                                    local_data(:,rep_list(r)+1:end)];
                            end
                        end
                        
                        if (length(local2global) ~= length(unique(local2global)))
                            error('Still repetitions to deal with')
                        end
                    end
                    
                    srate_scale = srate_eeg/srate_et;
                    if floor(srate_scale) == srate_scale
                        % EEG srate is a multiple of ET srate
                        for j = 1:srate_scale
                            data_et(:,local2global(1:end-j+1)+j-1)=local_data(:,1:end-j+1);
                        end
                    else
                        % TO DO
                        error('xdf2set not coded for when the eeg sampling rate is not a multiple of theet sampling rate.')
                    end
                    
                    %% chanlocs
                    for ch = 1:n_chans_et
                        chanlocs(ch+n_chans_eeg).type = 'ET';
                        data_type = AllStreams{ET_stream_inds(1)}.info.desc.channels.channel{ch}.type;
                        if strcmp(data_type, 'position')
                            chanlocs(ch+n_chans_eeg).unit = 'pixel';
                        elseif strcmp(data_type, 'area')
                            chanlocs(ch+n_chans_eeg).unit = 'a.u.';
                        elseif strcmp(data_type, 'ppd')
                            chanlocs(ch+n_chans_eeg).unit = 'pixelsPerDegree';
                        end
                        chanlocs(ch+n_chans_eeg).labels = AllStreams{ET_stream_inds(1)}.info.desc.channels.channel{ch}.label;
                    end
                    
                    %% etc
                    Data_all.etc.('descET') = AllStreams{ET_stream_inds(1)}.info.desc;
                    Data_all.etc.('infoET') = rmfield(AllStreams{ET_stream_inds(1)}.info,'desc');
                    disp ('...done')
                end
                
            case 3 %MOCAP
                if ~isempty(MOCAP_stream_inds)
                    disp('Importing MOCAP...')
                    %% Data
                    data_mocap = nan(n_chans_mocap,length(times));
                    for moc = 1:length(MOCAP_stream_inds)
                        if MOCAP_stream_inds(moc) > 0
                            local_data = AllStreams{MOCAP_stream_inds(moc)}.time_series;
                            local_times = AllStreams{MOCAP_stream_inds(moc)}.time_stamps;
                            
                            if moc == 1
                                channels_span = 1:size(local_data,1);
                            else
                                next_chan = channels_span(end)+1;
                                channels_span = next_chan:next_chan+size(local_data,1)-1;
                            end
                            
                            local2global = zeros(1,length(local_times));
                            ppm = ParforProgressbar(length(local_times), 'showWorkerProgress', true,...
                                'title', 'Finding time samples for MOCAP');
                            parfor t = 1:length(local_times)
                                ind = find(times>=local_times(t),1);
                                if isempty(ind)
                                    if t == length(local_times) && local_times(t)<times(end)+interval/2
                                        local2global(t) = length(times);
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                elseif ind==1
                                    if t == 1 && local_times(t)>times(1)-interval/2
                                        local2global(t) = 1;
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                else
                                    if round(local_times(t)-times(ind-1),4)<=interval/2
                                        local2global(t) = ind-1;
                                    elseif round(times(ind)-local_times(t),4)<=interval/2
                                        local2global(t) = ind;
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                end
                                ppm.increment();
                            end
                            delete(ppm);
                            
                            if (length(local2global) ~= length(unique(local2global)))
                                error('Still repetitions to deal with')
                            end
                            
                            % Fill gaps (because of lower srate) by
                            % replicating the data between real samples
                            for t = 1:length(local_times)-1
                                reps = local2global(t+1)-local2global(t);
                                data_mocap(channels_span,local2global(t):local2global(t+1)-1)=repmat(local_data(:,t),1,reps);
                            end
                            data_mocap(channels_span,local2global(end))=local_data(:,end);
                            
                            %% chanlocs
                            for ch = 1:length(channels_span)
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).type = 'MOCAP';
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).unit = AllStreams{MOCAP_stream_inds(moc)}.info.desc.channels.channel{ch}.unit;
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).labels = sprintf('RB%d_%s',moc,...
                                    AllStreams{MOCAP_stream_inds(moc)}.info.desc.channels.channel{ch}.label);
                            end
                            
                            %% etc
                            Data_all.etc.(['descMOCAP', num2str(moc)]) = AllStreams{MOCAP_stream_inds(moc)}.info.desc;
                            Data_all.etc.(['infoMOCAP', num2str(moc)]) = rmfield(AllStreams{MOCAP_stream_inds(moc)}.info,'desc');
                            
                        elseif MOCAP_stream_inds(moc) == 0
                            % In some cases, it is normal that the second mocap stream
                            % is missing (goggles OFF condition)
                            if moc == 1
                                channels_span = 1:(n_chans_mocap/2);
                            else
                                next_chan = channels_span(end)+1;
                                channels_span = next_chan:(next_chan+(n_chans_mocap/2)-1);
                            end
                            %% chanlocs
                            for ch = 1:length(channels_span)
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).type = 'MOCAP';
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).unit = mocap_channels_description.channel{ch}.unit;
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).labels = sprintf('RB%d_%s',moc,...
                                    mocap_channels_description.channel{ch}.label);
                            end
                            %% etc
                            Data_all.etc.(['descMOCAP', num2str(moc)]) = 'For compatibility only';
                            Data_all.etc.(['infoMOCAP', num2str(moc)]) = 'For compatibility only';
                        else
                            error('Unexpected error');
                        end
                    end
                    disp ('...done')
                end
                
            case 4 %Events
                disp('Importing Events...')
                if contains(study_config.study_folder, 'Tumbler')
                    all_events = export_events_Tumbler(AllStreams{Event_stream_inds}, times);
                    if strcmp(subject,'P1002')
                        % Button roles were inverted by this participant
                        pres_asw = strcmp({all_events.answer}, 'Present');
                        abs_asw = strcmp({all_events.answer}, 'Absent');
                        if ~isempty(pres_asw)
                            for e = find(pres_asw)
                                all_events(e).answer = 'Absent';
                            end
                        end
                        if ~isempty(abs_asw)
                            for e = find(abs_asw)
                                all_events(e).answer = 'Present';
                            end
                        end
                    end
                    
                    trialtypes_file = fullfile(input_filepath, [full_filename '.txt']);
                    all_events = read_trialtypes_Tumbler(all_events, trialtypes_file);
                elseif contains(study_config.study_folder, 'VEP')
                    all_events = export_events_VEP(AllStreams{Event_stream_inds}, times);
                end
                disp ('...done')
        end
    end
    
    data = data_eeg;
    if ~isempty(ET_stream_inds)
        data = cat(1, data, data_et);
    end
    if ~isempty(MOCAP_stream_inds)
        data = cat(1, data, data_mocap);
    end
    
    Data_all.data = data;
    [Data_all.nbchan,Data_all.pnts,Data_all.trials] = size(Data_all.data);
    [Data_all.filepath,fname,fext] = fileparts(input_file);
    Data_all.filename = [fname fext];
    Data_all.srate = srate_eeg;
    Data_all.xmin = 0;
    Data_all.xmax = (Data_all.pnts-1)/Data_all.srate;
    
    Data_all.chanlocs = chanlocs;
    
    Data_all.event = all_events;
    %Data_all.urevent = all_events;
    
    % Save file
    if ~isfolder(output_filepath)
        mkdir(output_filepath); % make sure that folder exists
    end
    pop_saveset(Data_all, 'filename',[file_prefix{1} study_config.filenames{i_filename}],'filepath', output_filepath);
    
    clear AllStreams data
    
    %% Split into unique channel types
    % get rid of memory mapped object storage
    %pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);
    
    disp('Splitting dataset into unique channel types...')
    Data_all = pop_loadset('filename',[file_prefix{1} study_config.filenames{i_filename}, '.set'],'filepath', output_filepath);
    
    unique_types = unique({Data_all.chanlocs.type});
    for unique_type = 1:length(unique_types)
        if strcmp(unique_types(unique_type),'EEG') && any(strcmp(unique_types,'EOG'))
            indices = find(strcmp({Data_all.chanlocs.type}, 'EEG') | strcmp({Data_all.chanlocs.type}, 'EOG'));
        elseif strcmp(unique_types(unique_type),'EOG')
            continue
        else
            indices = find(strcmp({Data_all.chanlocs.type}, unique_types(unique_type)));
        end
        file_prefix{end+1} = [subject '_' unique_types{unique_type} '_'];
        fprintf('Type "%s": %d of %d channels. ',unique_types{unique_type}, length(indices), length(Data_all.chanlocs))
        
        Data_split(unique_type) = pop_select(Data_all,'channel',indices);
        Data_split(unique_type) = eeg_checkset(Data_split(unique_type));
        
        % new data set in EEGLAB
        [ALLEEG, ~, CURRENTSET] = pop_newset(ALLEEG, Data_split(unique_type), CURRENTSET, 'gui', 'off');
        % save set
        pop_saveset(Data_split(unique_type), 'filename', [file_prefix{end} study_config.filenames{i_filename}], 'filepath', output_filepath);
    end
    disp('...done!')
    
    % clear RAM
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    clear Data_all Data_split
    disp(['XDF exporting done for ' full_filename]);
end

% This merges all EEG data files into one big file
input_filepath = output_filepath;

%file_prefix{1} = [subject '_AllStreams_'];

file_prefix = unique(file_prefix);
for f = 1:numel(file_prefix)
    % Only merge EEG file
    if contains(file_prefix{f},'EEG') || contains(file_prefix{f},'AllStreams')
        % get rid of memory mapped object storage
        %pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);
        
        % make sure EEGLAB has no files other than the ones to be merged
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        
        for i_filename = 1:length(study_config.filenames)
            full_filename = [file_prefix{f} study_config.filenames{i_filename} '.set'];
            if isfile([input_filepath filesep full_filename])
                EEG = pop_loadset('filename', full_filename, 'filepath', input_filepath);
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
            else
                continue
            end
        end
        
        % merges all files currently loaded in EEGLAB into one file and stores
        % the original filenames in EEG.etc.appended_files
        [ALLEEG, EEG_merged, CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,1:length(ALLEEG),...
            [file_prefix{f} study_config.merged_filename], output_filepath);
    end
end