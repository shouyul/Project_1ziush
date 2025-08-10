% bemobil_dipfit() - Prepares data for dipole fitting and runs the dipole fitting procedure
%
% Inputs:
%   EEG                     - current EEGLAB EEG structure
%   cfg_dipfit              - configuration containing the following options:
%       transform                   - already computed transform matrix
%                                   (skips coregistration - relevant if you don't have individualized data for electrodes location)
%       use_fiducials               - boolean indicating if coregistration should be based on fiducials
%                                   (only relevant if you digitized electrodes)
%       residualVariance_threshold  - number percentage of residual variance accepted, default is '15'
%       do_remove_outside_head      - 'on' or 'off' to remove dipoles located outside the head
%       number_of_dipoles           - '1' or '2', 2 meaning bilateral dipole fitting
%
% Outputs:
%   EEG                     - current EEGLAB EEG structure
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, coregister, pop_dipfit_settings, pop_multifit

function [EEG] = dipfit(EEG, cfg_dipfit)

% Read parameters from the config file:
trans = cfg_dipfit.transform;
use_fids = cfg_dipfit.use_fiducials;
RV_th = cfg_dipfit.residualVariance_threshold;
rm_out_head = cfg_dipfit.do_remove_outside_head;
num_dips = cfg_dipfit.number_of_dipoles;

% load some standard data for dipfit
dipfitdefs;
dipfit_path = fileparts(which('pop_multifit'));

if isempty(trans)
    
    % use standard BEM headmodel
    % change relevant fiducial labels so that it matches standard 5/10 template
    if use_fids  
        chaninfo = EEG.chaninfo;
        for f = 1:numel(chaninfo.nodatchans)
            switch chaninfo.nodatchans(f).labels
                case 'nas'
                    chaninfo.nodatchans(f).labels = 'Nz';
                case 'lhj'
                    chaninfo.nodatchans(f).labels = 'LPA';
                case 'rhj'
                    chaninfo.nodatchans(f).labels = 'RPA';
                otherwise
                    error('Unknown fiducial name')
            end
        end
        disp('Coregistering electrodes to 10-5 template...')
        error('To test')
        [~, trans] = coregister(EEG.chanlocs,...
            fullfile(dipfit_path, 'standard_BEM', 'elec','standard_1005.elc'), 'chaninfo1', chaninfo,...
            'manual', 'off', 'alignfid', {'Nz', 'LPA','RPA'});
    else
        % Manual fitting
        [~, trans] = coregister(EEG.chanlocs,...
            fullfile(dipfit_path, 'standard_BEM', 'elec','standard_1005.elc'));
    end
end

% Exclude EOG channels
chans2select = 1:EEG.nbchan;
EEG_indices = strcmp({EEG.chanlocs.type}, 'EEG');
chans2select = chans2select(EEG_indices);

% settings
EEG = pop_dipfit_settings(EEG, 'hdmfile',fullfile(dipfit_path, 'standard_BEM', 'standard_vol.mat'),...
    'coordformat','MNI',...
    'mrifile',fullfile(dipfit_path, 'standard_BEM', 'standard_mri.mat'),...
    'chanfile',fullfile(dipfit_path, 'standard_BEM', 'elec', 'standard_1005.elc'),...
    'coord_transform', trans, 'chansel', chans2select);

% do the dipole fitting
EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)], 'threshold', RV_th,...
    'dipoles', num_dips, 'rmout', rm_out_head);

% save dipfit info in EEG.etc
EEG.etc.dipfit.transform = trans;
EEG.etc.dipfit.RV_threshold = RV_th;
EEG.etc.dipfit.remove_outside_head = rm_out_head;
EEG.etc.dipfit.number_of_dipoles = num_dips;

EEG = eeg_checkset(EEG);