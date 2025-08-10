function new_events = read_trialtypes_Tumbler(events, trialtypes_file)
% For each block, import the block type (i.e. the condition) and the
% individual trial type from the specified text file.
%
% Inputs:
%   events              - event structure (output of export_events_Tumbler)
%   trialtypes_file     - [string] name of the trialtypes .txt file
%                           specifying trial types.
% Output:
%   new_events          - modified event structure with additionnal 
%                           'condition' and 'trialtype' columns.

fileID = fopen(trialtypes_file, 'r');
% First line is the condition
Cond = fgetl(fileID);
% Read next lines (trial types)
n_trials = 10;
TrialTypes = cell(1,n_trials);
for t = 1:n_trials
    TrialTypes{t} = fgetl(fileID);
end
fclose(fileID);

event_fields = fieldnames(events);
new_event_fields = [event_fields; {'condition'; 'trialtype'}];

command = '';
for f = 1:numel(new_event_fields)
    if any(strcmp({'type', 'answer', 'condition', 'trialtype'}, new_event_fields{f}))
        command = [command,'''',new_event_fields{f},''','''','];
    else
        command = [command,'''',new_event_fields{f},''',[],'];
    end
end
% Remove the last ','
command = command(1:end-1);
eval(['new_events = struct(',command,');']);

events_count = numel(events);
for e = 1:events_count
    for f = 1:numel(new_event_fields)
        if f <= numel(event_fields)
            if any(strcmp({'block','trial'},new_event_fields{f}))
                % Adapt numerotation
                new_events(e).(new_event_fields{f}) = events(e).(new_event_fields{f})+1;
            else
                new_events(e).(new_event_fields{f}) = events(e).(new_event_fields{f});
            end
        elseif strcmp('condition',new_event_fields{f})
            new_events(e).(new_event_fields{f}) = Cond;
        elseif strcmp('trialtype',new_event_fields{f})
            new_events(e).(new_event_fields{f}) = TrialTypes{new_events(e).trial};
        else
            error('Unknown field')
        end
    end
end

% Reorder fields:
new_events = orderfields(new_events, [1,2,3,4,8,5,9,6,7]);
end