function EEG = BehaviorReportVEP(EEG, cfg)
% Should be called after fillNaNs and events_check but before resampling or selection
% if all the trials are present
%
% Behavioral reports specific to the experiment
%
% Output: Creates EEG.etc.TrialData and also saves it externally

events = EEG.event;
evts_noBounds = events(~(strcmp({events.type},'boundary') | strcmp({events.type},'ExpEnd')));

DataBase = struct('ID', [], 'Gender', [], 'Block', [], 'Trial', [], 'TrialType', [],...
    'urevent_seq', [], 'EEGComplete_trial',[]);

% Subject specific variables
subject = cfg.subjects(cfg.current_subject).id;
gender = cfg.subjects(cfg.current_subject).gender;

TrialInspection = EEG.etc.TrialsInspection;

i=1;
for tr = 1:size(TrialInspection,1)
    bl = TrialInspection.BlockInd(tr);
    tr_inBl = TrialInspection.TrialInd(tr);
    trialEvents = [evts_noBounds.block] == bl & [evts_noBounds.trial] == tr_inBl;
    
    %fprintf('Block %d, Trial %d.\n', bl, tr_inBl)
    DataBase(i) = struct('ID', subject, 'Gender', gender,...
        'Block', bl, 'Trial', tr_inBl, 'TrialType', TrialInspection.TrialType(tr),...
        'urevent_seq', TrialInspection.Trial_urevent_seq(tr,:),...
        'EEGComplete_trial', TrialInspection.Trial_perc(tr));
    i=i+1;
end

DataBase = struct2table(DataBase);

%% Save the output structure
EEG.etc.TrialData = DataBase;
dirname = [cfg.study_folder cfg.preprocessing_folder];
fname = [subject '_DataBase.mat'];
save(fullfile(dirname, fname), 'DataBase')
end