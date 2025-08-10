function [PSD, nb_chunks] = computePSD(EEG, cfg)

if ~strcmp(cfg.epochs.window, 'fixed')
    error('Not coded yet')
end

nb_chans = EEG.nbchan;
nb_freqs = length(cfg.psd.FoI);
nb_trials = EEG.trials;

if cfg.psd.chunks == 0
    nb_chunks = 1;
else
    nb_chunks = (cfg.epochs.limits_wdw(2)-cfg.epochs.limits_wdw(1))/cfg.psd.chunks;
end

PSD = zeros(nb_chans, nb_freqs, nb_trials, nb_chunks);
switch cfg.psd.method
    case 'pwelch'
        % Compute PSD with pwelch method
        if nb_chunks == 1
            for tr = 1:nb_trials
                PSD(:,:,tr,1) = pwelch(squeeze(EEG.data(:,:,tr))',...
                    cfg.psd.window, cfg.psd.overlap, cfg.psd.FoI, EEG.srate, cfg.psd.output)';
            end
        else
            for ch = 1:nb_chunks
                EEG_sel = pop_select(EEG, 'time', cfg.epochs.limits_wdw(1)+[ch-1,ch].*cfg.psd.chunks);                
                for tr = 1:nb_trials
                    PSD(:,:,tr,ch) = pwelch(squeeze(EEG_sel.data(:,:,tr))',...
                        cfg.psd.window, cfg.psd.overlap, cfg.psd.FoI, EEG.srate, cfg.psd.output)';
                end
            end
        end
    otherwise
        error('Not coded yet')
end
end