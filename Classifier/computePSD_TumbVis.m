function [PSD, TumbVis, nb_chunks] = computePSD_TumbVis(EEG, phase, cfg)

if ~strcmp(cfg.epochs.window, 'fixed')
    error('Not coded yet')
end

mocapChans = contains({EEG.chanlocs.type}, 'MOCAP');
nb_chans = EEG.nbchan - sum(mocapChans);
nb_freqs = length(cfg.psd.FoI);
nb_trials = EEG.trials;

if cfg.psd.chunks == 0
    nb_chunks = 1;
else
    switch phase
        case 'EC'
            nb_chunks = abs(cfg.epochs.limits_wdw(1))/cfg.psd.chunks;
        case 'EO'
            nb_chunks = cfg.epochs.limits_wdw(2)/cfg.psd.chunks;
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

PSD = zeros(nb_chans, nb_freqs, nb_trials, nb_chunks);
switch cfg.psd.method
    case 'pwelch'
        % Compute PSD with pwelch method
        if nb_chunks == 1
            switch phase
                case 'EC'
                    EEG_sel = pop_select(EEG, 'time', [cfg.epochs.limits_wdw(1),0]);
                case 'EO'
                    EEG_sel = pop_select(EEG, 'time', [0,cfg.epochs.limits_wdw(2)]);
                otherwise
                    error('Unknown phase')
            end
            
            for tr = 1:nb_trials                
                PSD(:,:,tr,1) = pwelch(squeeze(EEG_sel.data(~mocapChans,:,tr))',...
                    cfg.psd.window, cfg.psd.overlap, cfg.psd.FoI, EEG.srate, cfg.psd.output)';
                if any(mocapChans)
                    TumbVis(tr,1) = mean(squeeze(EEG_sel.data(TumbVisChan,:,tr)));
                end
            end
        else
            for ch = 1:nb_chunks
                switch phase
                    case 'EC'
                        EEG_sel = pop_select(EEG, 'time', cfg.epochs.limits_wdw(1)+...
                            [ch-1,ch].*cfg.psd.chunks);
                    case 'EO'
                        EEG_sel = pop_select(EEG, 'time', [ch-1,ch].*cfg.psd.chunks);
                    otherwise
                        error('Unknown phase')
                end
                
                for tr = 1:nb_trials
                    PSD(:,:,tr,ch) = pwelch(squeeze(EEG_sel.data(~mocapChans,:,tr))',...
                        cfg.psd.window, cfg.psd.overlap, cfg.psd.FoI, EEG.srate, cfg.psd.output)';
                    if any(mocapChans)
                        TumbVis(tr,ch) = mean(squeeze(EEG_sel.data(TumbVisChan,:,tr)));
                    end
                end
            end
        end
    otherwise
        error('Not coded yet')
end
end