config_PIONEER_Tumbler;
dipfit_path = fileparts(which('pop_multifit'));

% Load empty dataset (to be able to load the proper chanlocs structure)
EEG = pop_loadset('filename', 'Empty_getChanlocs.set',...
    'filepath', fullfile(study_config.study_folder, study_config.raw_data_folder));
% Resize to 128 channels
EEG.nbchan = 128;
EEG.data = zeros(EEG.nbchan, EEG.pnts);
EEG.chanlocs = readlocs(fullfile(study_config.study_folder, study_config.raw_data_folder,...
    study_config.subjects(subject_ind).id, study_config.channel_locations_filename));
EEG = eeg_checkset(EEG);

EEG2 = EEG;
chanlocs = EEG2.chanlocs;
%matches = {'Z20Z','Iz';'Z2Z','Fpz';'Z10Z','Cz';'L12Z','CP1';'R12Z','CP2'};
matches = {'Z20Z','Iz';'Z2Z','Fpz';'L12Z','CP1';'R12Z','CP2'};

%% changing labels of the chanlocs struct for automatic warping
changed_chans = [];
% exact matches
for ce = 1:size(matches,1)
    chan_ind = find(strcmp({chanlocs.labels},matches{ce,1}));
    EEG2 = pop_chanedit(EEG2, 'changefield', {chan_ind, 'labels', matches{ce,2}});
    changed_chans = [changed_chans, chan_ind];
end
[~, trans] = coregister(EEG2.chanlocs,...
    fullfile(dipfit_path, 'standard_BEM', 'elec', 'standard_1005.elc'),...
    'warp', 'auto');


