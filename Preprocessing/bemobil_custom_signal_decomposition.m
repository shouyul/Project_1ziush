% bemobil_signal_decomposition() - Computes a spatial filter for the EEG data set, which decomposes
% the data into components (e.g. statistically independent components using AMICA - by default).
% AMICA sometimes crashes bevore the first iteration with some combinations of data set length and
% number of threads. Therefore, if this happens and the according error message of WINDOWS is
% closed, AMICA is automatically restarted with one thread less. In case the number of threads are
% reduced to 0, the standard EEGLAB runica() is started.
%
% Usage:
%   >>  [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET, amica, numb_models, maxx_threads, data_rank, other_algorithm)
%   >>  [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET, amica, numb_models, maxx_threads, data_rank, other_algorithm, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   amica                   - Boolean value (1/0) to use amica or not
%   num_models              - number of models to learn, default is 1
%   max_threads             - maximum of CPU threads to be used for AMICA
%   data_rank               - rank of the data matrix (= number of channels minus number of interpolated
%       channels minus 1, if average referenced)
%   other_algorithm         - currently (20.6.2017) not yet implemented, hopefully SSD and JD will
%       be added here some day
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%    EEGLAB, runamica15
%
% Authors: Lukas Gehrke, Marius Klug, 2017

function [EEG] = bemobil_custom_signal_decomposition(EEG, config, data_rank)

subject = config.subjects(config.current_subject).id;
N = makeFolderFileNames(config, subject);
main_filepath = N.searchFolder_2arch_rej;
out_filepath = fullfile(main_filepath, [subject '_' upper(config.ICAmethod)]);

crashed = false;
if strcmpi(config.ICAmethod, 'amica')
    if isfield(EEG,'datfile') && length(EEG.datfile) > 0
        disp('Found datfile.');
        data = fullfile(EEG.filepath, EEG.datfile);
    else
        disp('No datfile field found in EEG structure. Will write temp file in current directory.');
        data = EEG.data(:,:);
    end
    
    % delete potentially preexistent folder since it will interfere in case AMICA crashes
    if ~isempty(dir(out_filepath));
        rmdir(out_filepath,'s');
    end
    
    disp('Starting AMICA...');
    maxx_threads = config.max_threads;
    while maxx_threads > 0
        % try/catch loop because AMICA can crash dependent on the data set and the number of threads
        try
            [w, s, mods] = runamica15(data,...
                'num_models', config.num_models,...
                'max_threads', maxx_threads,...
                'outdir', out_filepath,...
                'num_chans', EEG.nbchan,...
                'writestep', config.max_iter,...
                'pcakeep', data_rank);
            disp('AMICA successfull, storing weights and sphere.');
            EEG.etc.spatial_filter.algorithm = 'AMICA';
            EEG.etc.spatial_filter.AMICAmods = mods;
            % if successful, get out of the loop
            break
            
        catch
            % if error, reduce threads by one
            maxx_threads = maxx_threads - 1;
            warning(['AMICA crashed. Reducing maximum threads to ' num2str(maxx_threads)]);
        end
    end
    
    if maxx_threads == 0
        warning('AMICA crashed with all possible maximum thread options. Try increasing the maximum usable threads of your CPU. If the maximum number of threads has already been tried, you''re pretty much fucked. Ask Jason Palmer, the creator of AMICA.');
        disp('Continuing with default EEGLAB runica() ...');
        crashed = true;
    end
end

if strcmpi(config.ICAmethod, 'runica') || crashed
    disp('Starting RUNICA...');
    [w,s] = runica(EEG.data, 'pca', data_rank, 'extended', 1, 'maxsteps', config.max_iter);
    disp('runica successfull, storing weights and sphere.');
    EEG.etc.spatial_filter.algorithm = 'RUNICA';
end

% store the actual weights and sphere values in the EEG data set and calculate the rest of the
% spatial filter stuff
EEG.icaweights = w;
EEG.icasphere = s;
EEG = eeg_checkset(EEG);

EEG.etc.spatial_filter.original_data_path = out_filepath;