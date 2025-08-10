function [EEG, bound_evts] = findBoundaries(EEG, study_config)
subj_ind = study_config.current_subject;
EEG = eeg_checkset(EEG, 'makeur');

events = EEG.event;
bound_evts = find(strcmp({events(:).type}, 'boundary'));

%% Clean around boundary if necessary (case by case)
if ~isempty(bound_evts)
    for b = length(bound_evts):-1:1
        events(bound_evts(b)).duration = 0;
        if strcmp(study_config.subjects(subj_ind).id, 'CIG01') && b == 1
            first2remove = find([events(:).TrialIndex] == 40,1);
            for e = bound_evts(b)-1:-1:first2remove
                events(e) = [];
            end
            bound_evts(b) = first2remove;
        end
    end
end

EEG.event = events;
EEG = eeg_checkset(EEG);
end