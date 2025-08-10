function [Amps, TumbVis] = computeAmps_TumbVis(EEG, phase, cfg)

if ~strcmp(cfg.epochs.window, 'fixed')
    error('Not coded yet')
end

mocapChans = contains({EEG.chanlocs.type}, 'MOCAP');
nb_chans = EEG.nbchan - sum(mocapChans);
nb_trials = EEG.trials;

if cfg.amps.chunks == 0
    nb_chunks = 1;
else
    switch phase
        case 'EC'
            nb_chunks = abs(cfg.epochs.limits_wdw(1))/cfg.amps.chunks;
        case 'EO'
            nb_chunks = cfg.epochs.limits_wdw(2)/cfg.amps.chunks;
        otherwise
            error('Unknown phase')
    end
end

if any(mocapChans)
    TumbVisChan = strcmp({EEG.chanlocs.labels}, 'TumblerVisibility');
    TumbVis = zeros(nb_trials, nb_chunks);
else
    TumbVis = nan(nb_trials, nb_chunks);
end

if nb_chunks == 1
    switch phase
        case 'EC'
            EEG_sel = pop_select(EEG, 'time', [cfg.epochs.limits_wdw(1),0]);
        case 'EO'
            EEG_sel = pop_select(EEG, 'time', [0,cfg.epochs.limits_wdw(2)]);
        otherwise
            error('Unknown phase')
    end
    
    nb_times = EEG_sel.pnts;
    Amps = zeros(nb_chans, nb_times, nb_trials, nb_chunks);
    
    Amps(:,:,:,1) = EEG_sel.data(~mocapChans,:,:);
    if any(mocapChans)
        TumbVis(:,1) = mean(squeeze(EEG_sel.data(TumbVisChan,:,:)),1)';
    end
else
    for ch = 1:nb_chunks
        switch phase
            case 'EC'
                EEG_sel = pop_select(EEG, 'time', cfg.epochs.limits_wdw(1)+...
                    [ch-1,ch].*cfg.amps.chunks);
            case 'EO'
                EEG_sel = pop_select(EEG, 'time', [ch-1,ch].*cfg.amps.chunks);
            otherwise
                error('Unknown phase')
        end
        
        if ch == 1
            nb_times = EEG_sel.pnts;
            Amps = zeros(nb_chans, nb_times, nb_trials, nb_chunks);
        end
        
        Amps(:,:,:,ch) = EEG_sel.data(~mocapChans,:,:);
        if any(mocapChans)
            TumbVis(:,ch) = mean(squeeze(EEG_sel.data(TumbVisChan,:,:)),1)';
        end
    end
end
end