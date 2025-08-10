function [events] = export_events_Tumbler(Event_stream, times)
% Interpret event XDF stream for PIONEER Tumbler experiment to put them in the EEGLAB set.
%
% Inputs:
%   - Event_stream :       Event stream loaded from the xdf.
%   - times :              Complete time points vector of the recording
%                               (to compute latencies)
%
% Outputs:
%   - events :              Structure containing all events loaded.

required_fields = {'type', 'latency', 'duration', 'block', 'trial', 'answer', 'delay'};
required_types = {'str', 'num', 'num', 'num', 'num', 'str', 'num'};
channels_correspondance = {'EventName', '', '', 'Block', 'Trial', 'Answer', 'Delay'};

required_inds = nan(numel(required_fields),1);
channels = Event_stream.info.desc.channels.channel;
for ch = 1:numel(channels)
    ind = find(strcmp(channels_correspondance, channels{ch}.label));
    if ~isempty(ind)
        required_inds(ind) = ch;
    end
end

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

for t=1:events_count
    for i=1:numel(required_inds)
        if strcmp(required_fields{i}, 'duration')
            continue
        elseif strcmp(required_fields{i}, 'latency')
            [~,events(t).latency] = min(abs(times - Event_stream.time_stamps(t)));
        else
            switch required_types{i}
                case 'str'
                    events(t).(required_fields{i}) = Event_stream.time_series{required_inds(i),t};
                case 'num'
                    events(t).(required_fields{i}) = str2num(Event_stream.time_series{required_inds(i),t});
            end
        end
    end
end
end

