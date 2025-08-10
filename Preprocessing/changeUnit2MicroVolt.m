function EEG = changeUnit2MicroVolt(EEG, cfg)
% Change EEG channels units

switch lower(cfg.recording_unit)
    case 'microvolt'
        coeff = 1;
    case 'millivolt'
        coeff = 10^3;
    case 'volt'
        coeff = 10^6;
    otherwise
        error('Unit not recognized')
end

% Select channels to change
chan_inds = contains({EEG.chanlocs.unit}, 'Volt');
EEG.data(chan_inds,:) = EEG.data(chan_inds,:).*coeff;

% Do the same for ICs (if exist)
% if ~isempty(EEG.icaact)
%     EEG.icaact = EEG.icaact.*coeff;
% end
end