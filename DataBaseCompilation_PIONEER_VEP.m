clear all;
config_PIONEER_VEP;

%% Compile some general data for further analyses
% Options for this script
recompute = false;
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
            case 'VEP'
                % load DataBase
                load(fullfile(filepath, sprintf('%s_DataBase.mat', subject)));
                
                % load prepared Dataset
                EEG_prepared = pop_loadset('filename', N.preparedFile, 'filepath', N.searchFolder_2, 'loadmode', 'info');
                %[ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
                events = EEG_prepared.event;
                noBounds = events(~(strcmp({events.type},'boundary') | strcmp({events.type},'ExpEnd')));
                
                % load Preprocessed dataset
                EEG_preproc = pop_loadset('filename', N.preICAFile, 'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
                badSamples = EEG_preproc.etc.(badSampsMethod).rejectedSamples;
                clear EEG_preproc
                
                fields = {'ID', 'Gender', 'Block', 'Trial', 'TrialType',...
                    'urevent_seq', 'completeEEG_full', 'cleanEEG_full'};
                if ~exist('MergedDataBase', 'var')
                    MergedDataBase = cell2struct(cell(numel(fields),1), fields', 1);
                    lastLine = 0;
                else
                    %lastLine = size(MergedDataBase,1);
                    lastLine = size(MergedDataBase,2); % JB: changed to 2, seemed a mistake (was OK in Tumbler version)
                end
                N_newTrials = size(DataBase,1);
                
                for i = 1:N_newTrials
                    blk = DataBase.Block(i);
                    trl = DataBase.Trial(i);                    
                    
                    % Check if this trial was kept in the prepared dataset
                    % (rejected if not EEG complete or the buffer condition was not satisfied)
                    if any([noBounds.block] == blk & [noBounds.trial] == trl) && ...
                        ~(strcmp(subject,'P1002') && blk==3 && trl == 1) && ... % adding exception for missing data
                        ~(strcmp(subject,'P1004-2') && blk==1 && trl == 1) % adding exception for missing data
                        % Compute total EEG cleanliness
                        start_lat = floor(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,1)).latency);
                        stop_lat = ceil(noBounds([noBounds.urevent] == DataBase.urevent_seq(i,end)).latency);
                        dur_samps = 1+stop_lat-start_lat;
                        percCleanEEG = 100*(sum(~badSamples(start_lat:stop_lat))/dur_samps);
                    else
                        % If incomplete, this trial was already rejected in the prepared dataset
                        percCleanEEG = NaN;
                    end
                    
                    values = {subject, DataBase.Gender(i), blk, trl, DataBase.TrialType(i),...
                        DataBase.urevent_seq(i,:), DataBase.EEGComplete_trial(i), percCleanEEG};
                    
                    %MergedDataBase(lastLine + i,:) = values;
                    MergedDataBase(lastLine + i) = cell2struct(values', fields', 1);
                end
                
                % Get rid of large structures:
                clear events noBounds
                clear EEG_prepared            
        end
    end
end

MergedDataBase = struct2table(MergedDataBase);

%% save
if ~exist(file2save) || recompute
    save(file2save,'MergedDataBase')
end