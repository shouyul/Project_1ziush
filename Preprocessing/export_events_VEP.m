function [events] = export_events_VEP(Event_stream, times)
% Interpret event XDF stream for PIONEER VEP experiment to put them in the EEGLAB set.
%
% Inputs:
%   - Event_stream :       Event stream loaded from the xdf.
%   - times :              Complete time points vector of the recording
%                               (to compute latencies)
%
% Outputs:
%   - events :              Structure containing all events loaded.

required_fields = {'type', 'latency', 'duration', 'block', 'trial', 'stimulus'};
required_types = {'str', 'num', 'num', 'num', 'num', 'str'};

events_count = length(Event_stream.time_stamps);
command = '';
for f = 1:numel(required_fields)
    if strcmp(required_fields{f}, 'duration')
        command = [command,'''',required_fields{f},''',num2cell(ones(1, events_count)),'];
    elseif strcmp(required_types{f}, 'str')
        command = [command,'''',required_fields{f},''','''','];
    elseif strcmp(required_types{f}, 'num')
        command = [command,'''',required_fields{f},''',[],'];
    end
end
% Remove the last ','
command = command(1:end-1);
eval(['events = struct(',command,');']);

%[~,events.latency] = min(abs(times' - Event_stream.time_stamps),[],1);
% events.type = cell(1,events_count);
% events.duration = ones(1,events_count);
% events.block = zeros(1,events_count);
% events.trial = zeros(1,events_count);
% events.stimulus = cell(1,events_count);

keys = {'block', 'trial', 'rep', 'event', 'type', 'stim', 'frame'};
copy_stim = false;
for e=1:events_count
    [~,events(e).latency] = min(abs(times - Event_stream.time_stamps(e)));
    
    pairs=strsplit(Event_stream.time_series{e}, ';');
    for k=1:length(keys)
        if (sum(contains(pairs, keys{k}))==1)
            pair = pairs(contains(pairs, keys{k}));
            splitPair = strsplit(pair{1}, ':');
            switch keys{k}
                case 'block'
                    events(e).block = str2double(splitPair{2});
                case 'trial'
                    events(e).trial = str2double(splitPair{2});
                case 'event'
                    events(e).type = splitPair{2};
                    if strcmp(events(e).type,'TrialStart')
                        copy_stim = true;
                    elseif strcmp(events(e).type,'TrialStop')
                        copy_stim = false;
                    end
                case 'type'
                    events(e).stimulus = splitPair{2};
                    if copy_stim
                        current_stim = events(e).stimulus;
                    end
                case 'frame'
                    fr = str2double(splitPair{2});
                    if fr == 1
                        events(e).type = [events(e).type, 'Show'];
                    elseif fr == 11
                        events(e).type = [events(e).type, 'Hide'];
                    end
            end
        end
    end
    
    if copy_stim
        events(e).stimulus = current_stim;
    end
end
end

