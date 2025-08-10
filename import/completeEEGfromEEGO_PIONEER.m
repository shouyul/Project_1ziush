function AllStreams = completeEEGfromEEGO_PIONEER(AllStreams, EEG_stream_inds, block, subjects_info, input_filepath, figs_folder)
Nstreams = length(EEG_stream_inds);

delims = strfind(input_filepath, filesep);
subj_name = input_filepath(delims(end)+1:end);
task_name = input_filepath(delims(end-3)+1:delims(end-2)-1);
sbj = strcmp({subjects_info(:).id}, subj_name);

if strcmp(task_name,'VEPTask')
    warning('The script should work for VEP task, given that the BIP photodiode channel was used on amp 1 = A');
elseif strcmp(task_name,'TumblerTask')
    error('The script should first be adapted to the Tumbler task (in terms of channels to look for)');
else
    error('Unexpected folder architecture');
end

% Load cnt file
list = dir(input_filepath);
cntfiles = contains({list.name},'.cnt');
if sum(cntfiles)>1
    fprintf('Doing %s\n', block);
    cntfiles = find(cntfiles);
    for f =1:length(cntfiles)
        fprintf('Option %d: %s\n', f, list(cntfiles(f)).name)
    end
    chosen_f = input('Which file to load? ');
    input_cntfile = list(cntfiles(chosen_f)).name;
else
    input_cntfile = list(cntfiles).name;
end
EEG_cnt = pop_loadeep_v4_custom(fullfile(input_filepath, input_cntfile));

% Analyse pattern of interruptions
cnt_events = EEG_cnt.event;
switch subj_name
%     case 'BAL03'
%         % Add a missing amplifier reconnected event
%         cnt_events(end+1) = struct('latency', 887008,...
%             'type', '9002, Amplifier reconnected' ,'duration',0);
%     case 'DRA22'
%         % Add a missing amplifier disconnected event
%         cnt_events(end+1) = struct('latency', 1810356,...
%             'type', '9001, Amplifier disconnected' ,'duration',0);
%     case 'ECU24'
%         % Add a missing amplifier reconnected event
%         cnt_events(end+1) = struct('latency', 2201464.9,...
%             'type', '9002, Amplifier reconnected' ,'duration',0);
end
% Sorting the events based on latency (in case some were added at the end)
[~,inds] = sort([cnt_events(:).latency]);
cnt_events = cnt_events(inds);
% Rounding if scalar latencies were used for organization purposes
for e = 1:length(cnt_events)
    cnt_events(e).latency = round([cnt_events(e).latency]);
end

b = str2num(block(end-2:end));
%a = str2num(block(end));
if strcmp(subjects_info(sbj).blockInterruption, 'No')
    inter = b-1;
% elseif strcmp(subjects_info(sbj).blockInterruption, 'B1')
%     if b == 1
%         inter = b-1+a-1;
%     else
%         inter = b;
%     end
% elseif strcmp(subjects_info(sbj).blockInterruption, 'B2')
%     if b == 1
%         inter = b-1;
%     elseif b == 2
%         inter = b-1+a-1;
%     else
%         inter = b;
%     end
% elseif strcmp(subjects_info(sbj).blockInterruption, 'B3')
%      if b < 3
%         inter = b-1;
%     else
%         inter = b-1+a-1;
%     end
else
    inter = NaN;
end

blockStart_lat = NaN;
blockEnd_lat = NaN;
segments = [];
seg2skip = 0;
userInterruptions = 0;
ampliDisconnected = 0;
paired = true;
e = 1;
while e <= length(cnt_events)
    if contains(cnt_events(e).type, '0,')
        % Impedance event
        if isnan(blockStart_lat) && paired
            if userInterruptions == inter
                blockStart_lat = cnt_events(e).latency;
                segments(ampliDisconnected+1,1) = blockStart_lat;
            end
            paired = false;
        elseif (e == length(cnt_events) || contains(cnt_events(e+1).type, '0,'))...
                && isnan(blockEnd_lat) && ~paired
            if userInterruptions == inter
                blockEnd_lat = cnt_events(e).latency;
                segments(ampliDisconnected+1,2) = blockEnd_lat;
            else
                seg2skip = seg2skip+1;
            end
            paired = true;
            userInterruptions = userInterruptions+1;
        end
    elseif contains(cnt_events(e).type, '9001,')
        % Amplifier disconnected event
        if userInterruptions == inter
            segments(ampliDisconnected+1,2) = cnt_events(e).latency;
            ampliDisconnected = ampliDisconnected+1;
        else
            seg2skip = seg2skip+1;
        end
    elseif contains(cnt_events(e).type, '9002,')
        % Amplifier reconnected event
        if userInterruptions == inter
            segments(ampliDisconnected+1,1) = cnt_events(e).latency;
        end
        e = e+2;
        continue
    end
    
    if ~isnan(blockStart_lat) && ~isnan(blockEnd_lat)
        break
    end
    
    e = e+1;
end

%% Plot the segments
figure
hold on
for seg = 1:size(segments,1)
    plot(EEG_cnt.times(segments(seg,:))./1000,[1,1], '-+r')
    text(EEG_cnt.times(segments(seg,1))./1000,1.1,num2str(seg));
    if seg > 1
        xline(EEG_cnt.times(segments(seg,1))./1000, ':r');
    end
end

center = (EEG_cnt.times(segments(end,2)) + EEG_cnt.times(segments(1,1)))/(1000*2);
switch Nstreams
    case 1
        center1 = (AllStreams{EEG_stream_inds(1)}.segments(1).t_begin +...
            AllStreams{EEG_stream_inds(1)}.segments(end).t_end)/2;
        
        for seg = 1:length(AllStreams{EEG_stream_inds(1)}.segments)
            plot([AllStreams{EEG_stream_inds(1)}.segments(seg).t_begin, AllStreams{EEG_stream_inds(1)}.segments(seg).t_end]...
                + center - center1,[2,2]-0.05*(seg-1), '-+b')
            text(AllStreams{EEG_stream_inds(1)}.segments(seg).t_begin + center - center1, 2.1-0.05*(seg-1), num2str(seg));
        end
        ylim([0,3])
        yticks(1:2)
        yticklabels({'Eego recording', 'LSL stream'})
    case 2
        center1 = (AllStreams{EEG_stream_inds(1)}.segments(1).t_begin +...
            AllStreams{EEG_stream_inds(1)}.segments(end).t_end)/2;
        center2 = (AllStreams{EEG_stream_inds(2)}.segments(1).t_begin +...
            AllStreams{EEG_stream_inds(2)}.segments(end).t_end)/2;
        
        for seg = 1:length(AllStreams{EEG_stream_inds(1)}.segments)
            plot([AllStreams{EEG_stream_inds(1)}.segments(seg).t_begin, AllStreams{EEG_stream_inds(1)}.segments(seg).t_end]...
                + center - center1,[3,3]-0.05*(seg-1), '-+b')
            text(AllStreams{EEG_stream_inds(1)}.segments(seg).t_begin + center - center1, 3.1-0.05*(seg-1), num2str(seg));
        end
        
        for seg = 1:length(AllStreams{EEG_stream_inds(2)}.segments)
            plot([AllStreams{EEG_stream_inds(2)}.segments(seg).t_begin, AllStreams{EEG_stream_inds(2)}.segments(seg).t_end]...
                + center - center2,[2,2]-0.05*(seg-1), '-+g')
            text(AllStreams{EEG_stream_inds(2)}.segments(seg).t_begin + center - center2, 2.1-0.05*(seg-1), num2str(seg));
        end
        ylim([0,4])
        yticks(1:3)
        yticklabels({'Eego recording', 'LSL stream 2', 'LSL stream 1'})
end
xlabel('Time (s)')
title(sprintf('%s - %s', subj_name, block))

saveCurrentFig([figs_folder, 'MissingData', filesep],...
    sprintf('SegmentsRepartition_%s_%s', subj_name, block),{'png'}, [800,600])
% set(gcf, 'Position', [200, 300, 800, 600]);
% saveas(gcf, fullfile('C:\Users\VRuser\Documents\MATLAB\DATA\figures\MissingData',...
%     sprintf('SegmentsRepartition_%s_%s.png', subj_name, block)), 'png');

if b == 1
    % Enable breakpoint to check the cnt_events struct for each subject
    cnt_events;
end

%% Doing the stuff
min_err_thresh = 0.0021; % in microvolt
step_seg_search = 1; % in seconds
%step_checkpoints = 60; % in seconds

for i = 1:Nstreams
    fprintf('Doing Stream %d\n',i);
    stream = AllStreams{EEG_stream_inds(i)};
    chans2match = 1:size(stream.time_series,1)-2; % from stream
    dur2match = 1; % in seconds
    
    if exist('newStream', 'var')
        clear newStream
    end
    newStream.info = stream.info;
    
    % Fill with copied segment that should start the stream (if exists)
    if i == 2 && exist('segs2copy', 'var')
        seg2 = find(floor(segs2copy)==0);
        if ~isempty(seg2)
            newStream.segments = newStream2.segments(seg2);
            newStream.segments(1).index_range = [1,newStream.segments(1).num_samples];
            newStream.time_series = newStream2.time_series(:,...
                newStream2.segments(seg2).index_range(1):newStream2.segments(seg2).index_range(2));
            newStream.time_stamps = newStream2.time_stamps(...
                newStream2.segments(seg2).index_range(1):newStream2.segments(seg2).index_range(2));
        end
    end
    
    if strcmp(task_name,'VEPTask')
        if i == 1
            chans2search = [1:64,129]; % from EEG_cnt, should be 64 channels + bip channel
            %segs2search = segs2search1; % from EEG_cnt
        else
            chans2search = 65:128; % from EEG_cnt, other EEG channels
            %segs2search = segs2search2; % from EEG_cnt
        end
    else
        error('update script');
    end
    
    for seg = 1:length(stream.segments)
        %% Asking user decision
        skipped = false;
        kept = false;
        while true
            instruction = input(sprintf('What should I do with segment %d (Stream %d)? ',seg,i));
            if ischar(instruction)
                switch instruction
                    case 'skip'
                        % Skip this segment
                        skipped = true;
                        fprintf('Skipping Seg %d\n', seg);
                        break
                    case 'keep'
                        kept = true;
                        fprintf('Keeping Seg %d\n', seg);
                        break
                    case 'exp'
                        % Expand the segment
                        seg2search = input('Segment correspondance with cnt? ');
                        seg2copy = [];
                        break
                    case 'copy'
                        % Expand the segment and then copy the corresponding channels
                        % to the other stream
                        seg2search = input('Segment correspondance with cnt? ');
                        seg2copy = input('Copy to which segment index (final numerotation)? ');
                        break
                    otherwise
                        fprintf('%s not understood\n', instruction)
                        continue
                end
            else
                disp('Instruction should be a char')
                continue
            end
        end
        
        if ~skipped
            if ~kept
                %% Check if there is a real need to search this segment
                % 1. Needs to be copied to the other stream
                % 2. The cnt data contains significantly more samples (>1%)
                % AND the number of samples should be more than the number of samples to match
                n_samps = stream.segments(seg).num_samples;
                n_samps2match = round(dur2match * stream.segments(seg).effective_srate);
                reuse = false;
                if ~isempty(seg2copy) || diff(segments(seg2search,:)) > n_samps * 1.01...
                        && (n_samps >= n_samps2match)
                    fprintf('Doing seg %d\n', seg)
                    
                    %% Match segment beginning
                    ind_start_seg = stream.segments(seg).index_range(1);
                    % Search for the most suitable start index
                    disp('Searching for the most suitable starting index for segment beginning:')
                    ind_start_search = 1;
                    wdw_range = [-30, 30];
                    found = false;
                    max_ind_bad = 0;
                    min_ind_good = inf;
                    
                    % Storage variables
                    ind_start = [];
                    min_errors_start = [];
                    start_align = [];
                    num_valid_start = [];
                    ind = 1;
                    while max_ind_bad ~= min_ind_good-1
                        fprintf('ind_start_search = %d\n',ind_start_search);
                        % Samples to match
                        samples2match = (ind_start_seg+ind_start_search-1)+(0:(n_samps2match-1));
                        % Search window
                        searchWindow = (stream.time_stamps(ind_start_seg+ind_start_search-1) + ...
                            - stream.segments(seg).t_begin) + wdw_range; % in seconds
                        
                        % Define data to match
                        data2match = stream.time_series(chans2match,samples2match).*10^6; % conversion to microVolt as in cnt file
                        % Perform the search
                        [min_errors_start(ind), start_align(ind), align_time, num_valid_start(ind)] =...
                            findStreamAlignmentFromCNT(EEG_cnt, data2match, chans2search,...
                            segments(seg2search,:), searchWindow, 'backward', ind==1);
                        
                        % Determine whether min_err is sufficiently low
                        ind_start(ind) = ind_start_search;
                        if min_errors_start(ind) >= min_err_thresh
                            if ind == 1
                                % Ask user if he prefers to extend the wdw range
                                % or change the starting index
                                userAnswer1 = input('Extend wdw_range? ');
                                if userAnswer1
                                    wdw_range = input('New wdw_range? ');
                                    continue
                                end
                            end
                            max_ind_bad = ind_start_search;
                            if found
                                ind_start_search = ceil((max_ind_bad+min_ind_good)/2);
                            else
                                max_step = n_samps-(ind_start_search+n_samps2match);
                                if max_step > step_seg_search*EEG_cnt.srate
                                    ind_start_search = ind_start_search + (step_seg_search*EEG_cnt.srate);
                                elseif max_step > 1
                                    ind_start_search = ind_start_search + floor(max_step/2);
                                elseif max_step == 1
                                    ind_start_search = ind_start_search + 1;
                                else
                                    break
                                end
                            end
                        else
                            if ~found
                                found = true;
                            end
                            min_ind_good = ind_start_search;
                            ind_start_search = floor((max_ind_bad+min_ind_good)/2);
                            
                            % Restrain wdw_range to gain computation time
                            % align_time is computed with respect to the center of
                            % the current wdw_range
                            align_time = align_time + sum(wdw_range)/2;
                            
                            
                            %wdw_range = align_time + [-5,5];
                            % replaces the above line
                            if ind_start_seg+ind_start_search-1 > 0 % because of some exception, didn't try to get to the bottom of it
                                if (stream.time_stamps(ind_start_seg+ind_start_search-1) ...
                                        - stream.segments(seg).t_begin) + align_time < 0
                                    wdw_range = - (stream.time_stamps(ind_start_seg+ind_start_search-1) ...
                                        - stream.segments(seg).t_begin) + [-5,5];
                                end
                            else
                                wdw_range = align_time + [-5,5];
                            end
                            
                        end
                        ind = ind+1;
                    end
                    
                    if max_ind_bad ~= min_ind_good-1
                        % The while loop didn't ended correctly
                        disp('Failed to find a good fit for segment beginning');
                        reuse = true;
                    end
                    
                    if ~reuse
                        good_ind = find(ind_start == min_ind_good);
                        num_samps_add_begin = num_valid_start(good_ind);
                        start_ind_seg_cnt = start_align(good_ind);
                        % For expansion
                        if num_samps_add_begin < 0
                            start_ind_seg_stream = ind_start_seg;
                            num_samps_add_begin = 0;
                        else
                            start_ind_seg_stream = ind_start_seg+ind_start(good_ind)-1;
                        end
                        % For copy
                        if ~isempty(seg2copy)
                            start_ind_seg_stream_copy = ind_start_seg+ind_start(good_ind)-1;
                        end
                        
                        if ind > 2
                            % Otherwise no search has been conducted
                            [ind_start, sort_order] = sort(ind_start);
                            min_errors_start = min_errors_start(sort_order);
                            start_align = start_align(sort_order);
                            num_valid_start = num_valid_start(sort_order);
                            
                            % Plot the result of the search
                            figure
                            semilogy(ind_start,min_errors_start, '-x');
                            xlabel('ind start search');
                            ylabel('min err (in microVolt)');
                            title(sprintf('Stream %d - Seg %d - beginning', i, seg));
                        end
                        
                        %% Match segment ending
                        ind_end_seg = stream.segments(seg).index_range(2);
                        % Search for the most suitable start index
                        disp('Searching for the most suitable ending index for segment ending:')
                        ind_end_search = 0;
                        if ind > 2
                            % reset wdw_range
                            wdw_range = [-30, 30];
                        end
                        found = false;
                        min_ind_bad = 1;
                        max_ind_good = -inf;
                        
                        % Storage variables
                        ind_end = [];
                        min_errors_end = [];
                        end_align = [];
                        num_valid_end = [];
                        ind = 1;
                        while min_ind_bad ~= max_ind_good+1
                            fprintf('ind_end_search = %d\n', ind_end_search);
                            % Samples to match
                            samples2match = (ind_end_seg+ind_end_search)+(-(n_samps2match-1):0);
                            % Search window
                            searchWindow = (stream.time_stamps(ind_end_seg+ind_end_search) + ...
                                - stream.segments(seg).t_begin) + wdw_range; % in seconds
                            
                            % Define data to match
                            data2match = stream.time_series(chans2match,samples2match).*10^6; % conversion to microVolt as in cnt file
                            % Perform the search
                            [min_errors_end(ind), end_align(ind), align_time, num_valid_end(ind)] =...
                                findStreamAlignmentFromCNT(EEG_cnt, data2match, chans2search,...
                                segments(seg2search,:), searchWindow, 'forward', ind==1);
                            
                            % Determine whether min_err is sufficiently low
                            ind_end(ind) = ind_end_search;
                            if min_errors_end(ind) >= min_err_thresh
                                if ind == 1
                                    % Ask user if he prefers to extend the wdw range
                                    % or change the starting index
                                    userAnswer1 = input('Extend wdw_range? ');
                                    if userAnswer1
                                        wdw_range = input('New wdw_range? ');
                                        continue
                                    end
                                end
                                min_ind_bad = ind_end_search;
                                if found
                                    ind_end_search = floor((min_ind_bad+max_ind_good)/2);
                                else
                                    max_step = n_samps + ind_end_search - n_samps2match;
                                    if max_step > step_seg_search*EEG_cnt.srate
                                        ind_end_search = ind_end_search - (step_seg_search*EEG_cnt.srate);
                                    elseif max_step > 1
                                        ind_end_search = ind_end_search - floor(max_step/2);
                                    elseif max_step == 1
                                        ind_end_search = ind_end_search - 1;
                                    else
                                        break
                                    end
                                end
                            else
                                if ~found
                                    found = true;
                                end
                                max_ind_good = ind_end_search;
                                ind_end_search = ceil((min_ind_bad+max_ind_good)/2);
                                % Restrain wdw_range to gain computation time
                                % align_time is computed with respect to the center of
                                % the current wdw_range
                                align_time = align_time + sum(wdw_range)/2;
                                wdw_range = align_time + [-5,5];
                            end
                            ind = ind+1;
                        end
                        
                        if min_ind_bad ~= max_ind_good+1
                            % The while loop didn't ended correctly
                            disp('Failed to find a good fit for segment ending');
                            reuse = true;
                        end
                        
                        if ~reuse
                            good_ind = find(ind_end == max_ind_good);
                            num_samps_add_end = num_valid_end(good_ind);
                            stop_ind_seg_cnt = end_align(good_ind);
                            % For expansion
                            if num_samps_add_end < 0
                                stop_ind_seg_stream = ind_end_seg;
                                num_samps_add_end = 0;
                            else
                                stop_ind_seg_stream = ind_end_seg+ind_end(good_ind);
                            end
                            %For copy
                            if ~isempty(seg2copy)
                                stop_ind_seg_stream_copy = ind_end_seg+ind_end(good_ind);
                            end
                            
                            if ind > 2
                                % Otherwise no search has been conducted
                                [ind_end, sort_order] = sort(ind_end);
                                min_errors_end = min_errors_end(sort_order);
                                end_align = end_align(sort_order);
                                num_valid_end = num_valid_end(sort_order);
                                
                                % Plot the result of the search
                                figure
                                semilogy(ind_end,min_errors_end, '-x');
                                xlabel('ind end search');
                                ylabel('min err (in microVolt)');
                                title(sprintf('Stream %d - Seg %d - ending', i, seg));
                            end
                            
                            %% Expand the segment and save in the newStream structure
                            n_samps_tot = num_samps_add_begin +...
                                (stop_ind_seg_stream - start_ind_seg_stream + 1)...
                                + num_samps_add_end;
                            t_b = stream.time_stamps(start_ind_seg_stream)-(num_samps_add_begin/EEG_cnt.srate);
                            t_e = stream.time_stamps(stop_ind_seg_stream)+(num_samps_add_end/EEG_cnt.srate);
                            
                            % Segments struct
                            Segment = struct('num_samples', n_samps_tot, 'index_range', [1,n_samps_tot],...
                                't_begin',t_b, 't_end', t_e, 'duration', t_e-t_b,...
                                'effective_srate', n_samps_tot/(t_e-t_b));
                            
                            % Time stamps
                            Time_stamps = [stream.time_stamps(start_ind_seg_stream)+(-num_samps_add_begin:-1)./EEG_cnt.srate,...
                                stream.time_stamps(start_ind_seg_stream:stop_ind_seg_stream),...
                                stream.time_stamps(stop_ind_seg_stream)+(1:num_samps_add_end)./EEG_cnt.srate];
                            
                            % Time series
                            Time_series = zeros(size(stream.time_series,1),n_samps_tot);
                            if num_samps_add_begin > 0
                                Time_series(chans2match,1:num_samps_add_begin) =...
                                    EEG_cnt.data(chans2search, start_ind_seg_cnt + (-num_samps_add_begin:-1)).*10^(-6);
                            end
                            Time_series(:,num_samps_add_begin+1+(0:(stop_ind_seg_stream-start_ind_seg_stream))) =...
                                stream.time_series(:,start_ind_seg_stream:stop_ind_seg_stream);
                            if num_samps_add_end > 0
                                Time_series(chans2match,(end-num_samps_add_end+1):end) =...
                                    EEG_cnt.data(chans2search, stop_ind_seg_cnt + (1:num_samps_add_end)).*10^(-6);
                            end
                            
                            if ~isfield(newStream,'segments')
                                newStream.segments = Segment;
                                newStream.time_series = Time_series;
                                newStream.time_stamps = Time_stamps;
                            else
                                num_existing_segments = numel(newStream.segments);
                                num_existing_samp = length(newStream.time_stamps);
                                Segment.index_range = num_existing_samp + Segment.index_range;
                                newStream.segments(num_existing_segments+1) = Segment;
                                newStream.time_series = cat(2, newStream.time_series, Time_series);
                                newStream.time_stamps = [newStream.time_stamps, Time_stamps];
                            end
                            
                            fprintf('Segment %d enlarged by %.0f%%\n', seg,...
                                (100*Segment.num_samples/stream.segments(seg).num_samples)-100);
                            
                            
                            %% Fill with copied segment (if exists)
                            if i == 2 && exist('segs2copy', 'var')
                                num_existing_segments = numel(newStream.segments);
                                seg2 = find(floor(segs2copy)==num_existing_segments);
                                if ~isempty(seg2)
                                    num_existing_samp = length(newStream.time_stamps);
                                    Segment = newStream2.segments(seg2);
                                    Segment.index_range = num_existing_samp + [1,Segment.num_samples];
                                    newStream.segments(num_existing_segments+1) = Segment;
                                    newStream.time_series = cat(2, newStream.time_series, newStream2.time_series(:,...
                                        newStream2.segments(seg2).index_range(1):newStream2.segments(seg2).index_range(2)));
                                    newStream.time_stamps = [newStream.time_stamps, newStream2.time_stamps(...
                                        newStream2.segments(seg2).index_range(1):newStream2.segments(seg2).index_range(2))];
                                end
                            end
                            
                            %% Copying data to the other stream
                            if ~isempty(seg2copy)
                                % Check if numsamples interpolated are consistent
                                numSampsinStream = (stop_ind_seg_stream_copy - start_ind_seg_stream_copy + 1);
                                numSampsinCNT = (stop_ind_seg_cnt - start_ind_seg_cnt + 1);
                                if numSampsinStream ~= numSampsinCNT
                                    warning('Number of corresponding samples:\n In original stream: %d\n In CNT: %d',...
                                        numSampsinStream, numSampsinCNT);
                                end
                                
                                start_ind_copy = start_ind_seg_cnt - num_samps_add_begin;
                                stop_ind_copy = stop_ind_seg_cnt + num_samps_add_end;
                                
                                % Reuse Segment and Time_stamps
                                Segment.index_range = [1 Segment.num_samples];
                                % Recreate Time_series by interpolation (in case there is
                                % no exact correspondance
                                x = (EEG_cnt.times(start_ind_copy:stop_ind_copy) - EEG_cnt.times(start_ind_seg_cnt))/1000;
                                v = EEG_cnt.data(setdiff(1:EEG_cnt.nbchan,chans2search),start_ind_copy:stop_ind_copy).*10^(-6);
                                xq = Time_stamps - Time_stamps((start_ind_seg_stream_copy - start_ind_seg_stream) + num_samps_add_begin + 1);
                                Time_series = zeros(size(v,1)+2,length(xq));
                                Time_series(1:end-2,:) = interp1(x, v', xq, 'pchip')';
                                
                                % Save data
                                if i == 1
                                    % newStream for stream 2 is yet to be created
                                    if ~exist('newStream2', 'var')
                                        newStream2.info = AllStreams{EEG_stream_inds(2)}.info;
                                        newStream2.segments = Segment;
                                        newStream2.time_series = Time_series;
                                        newStream2.time_stamps = Time_stamps;
                                        segs2copy = seg2copy;
                                    else
                                        num_existing_segments = numel(newStream2.segments);
                                        num_existing_samp = length(newStream2.time_stamps);
                                        Segment.index_range = num_existing_samp + Segment.index_range;
                                        newStream2.segments(num_existing_segments+1) = Segment;
                                        newStream2.time_series = cat(2, newStream2.time_series, Time_series);
                                        newStream2.time_stamps = [newStream2.time_stamps, Time_stamps];
                                        segs2copy = [segs2copy, seg2copy];
                                    end
                                else
                                    % newStream for stream 1 has already been computed
                                    stream2complete = AllStreams{EEG_stream_inds(1)};
                                    newStream1 = stream2complete;
                                    
                                    if isfield(stream2complete, 'segments')
                                        if floor(seg2copy) < numel(stream2complete.segments)
                                            % Delete data that should fall after insertion
                                            newStream1.segments = newStream1.segments(1:floor(seg2copy));
                                            newStream1.time_series = newStream1.time_series(:,...
                                                1:stream2complete.segments(floor(seg2copy)).index_range(2));
                                            newStream1.time_stamps = newStream1.time_stamps(...
                                                1:stream2complete.segments(floor(seg2copy)).index_range(2));
                                        end
                                        
                                        % Insert the copied data
                                        num_existing_samp = length(newStream1.time_stamps);
                                        Segment.index_range = num_existing_samp + Segment.index_range;
                                        newStream1.segments(ceil(seg2copy)) = Segment;
                                        newStream1.time_series = cat(2, newStream1.time_series, Time_series);
                                        newStream1.time_stamps = [newStream1.time_stamps, Time_stamps];
                                        
                                        if floor(seg2copy) < numel(stream2complete.segments)
                                            % Re-insert the data deleted
                                            num_existing_samp = length(newStream1.time_stamps);
                                            for seg2 = ceil(seg2copy):numel(stream2complete.segments)
                                                num_existing_segments = numel(newStream1.segments);
                                                Segment = stream2complete.segments(seg2);
                                                Segment.index_range = num_existing_samp + [1, Segment.num_samples];
                                                newStream1.segments(num_existing_segments+1) = Segment;
                                                num_existing_samp = num_existing_samp + Segment.num_samples;
                                            end
                                            newStream1.time_series = cat(2, newStream1.time_series,...
                                                stream2complete.time_series(:, stream2complete.segments(ceil(seg2copy)).index_range(1):end));
                                            newStream1.time_stamps = [newStream1.time_stamps,...
                                                stream2complete.time_stamps(stream2complete.segments(ceil(seg2copy)).index_range(1):end)];
                                        end
                                    else
                                        % rare case where no segments were kept in stream 1
                                        newStream1.segments = Segment;
                                        newStream1.time_series = Time_series;
                                        newStream1.time_stamps = Time_stamps;
                                    end
                                    
                                    AllStreams{EEG_stream_inds(1)} = newStream1;
                                end
                            end
                        end
                    end
                else
                    reuse = true;
                end
            else
                reuse = true;
            end
            
            if reuse
                %% Copy existing segment
                fprintf('Reusing seg %d as such\n', seg)
                if ~isfield(newStream,'segments')
                    newStream.segments = stream.segments(seg);
                    newStream.segments(1).index_range = [1, stream.segments(seg).num_samples];
                    newStream.time_series = stream.time_series(:,...
                        stream.segments(seg).index_range(1):stream.segments(seg).index_range(2));
                    newStream.time_stamps = stream.time_stamps(...
                        stream.segments(seg).index_range(1):stream.segments(seg).index_range(2));
                else
                    %% First fill with segment copied from other stream (if exists)
                    if i == 2 && exist('segs2copy', 'var')
                        num_existing_segments = numel(newStream.segments);
                        seg2 = find(floor(segs2copy)==num_existing_segments);
                        if ~isempty(seg2)
                            num_existing_samp = length(newStream.time_stamps);
                            Segment = newStream2.segments(seg2);
                            Segment.index_range = num_existing_samp + [1,Segment.num_samples];
                            newStream.segments(num_existing_segments+1) = Segment;
                            newStream.time_series = cat(2, newStream.time_series, newStream2.time_series(:,...
                                newStream2.segments(seg2).index_range(1):newStream2.segments(seg2).index_range(2)));
                            newStream.time_stamps = [newStream.time_stamps, newStream2.time_stamps(...
                                newStream2.segments(seg2).index_range(1):newStream2.segments(seg2).index_range(2))];
                        end
                    end
                    
                    num_existing_segments = numel(newStream.segments);
                    num_existing_samp = length(newStream.time_stamps);
                    newStream.segments(num_existing_segments+1) = stream.segments(seg);
                    newStream.segments(num_existing_segments+1).index_range = ...
                        num_existing_samp + [1, stream.segments(seg).num_samples];
                    newStream.time_series = cat(2, newStream.time_series, stream.time_series(:,...
                        stream.segments(seg).index_range(1):stream.segments(seg).index_range(2)));
                    newStream.time_stamps = [newStream.time_stamps, stream.time_stamps(...
                        stream.segments(seg).index_range(1):stream.segments(seg).index_range(2))];
                end
            end
        else
            if i == 2 && seg == length(stream.segments) && ~isfield(newStream,'segments')
            % rare case where no segments were kept in stream 2 but there are segments to copy from stream 1
            newStream = newStream2;
            end
        end
    end
    AllStreams{EEG_stream_inds(i)} = newStream;
end
close all
end