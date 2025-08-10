function [filteredEEG_noLN, lineNoiseOut] = removeLineNoise_custom(filteredEEG, lineNoiseRemoval_method, plot)
% Removes Line Noise depending on the specified 'lineNoiseRemoval_method'
% The input EEG data should be filtered
% Outputs are the EEG data with performed line noise removal and
% lineNoiseOut struct containing information about the parameters used

lineNoiseIn = SetLineNoiseRemovalStruct(filteredEEG);
switch lineNoiseRemoval_method
    case 'cleanLine'
        [filteredEEG_noLN, ~, lineNoiseOut, Spect_orig, Spect_clean, eval_freqs] = pop_cleanline(filteredEEG,...
            'LineFrequencies',lineNoiseIn.lineFrequencies,...
            'ScanForLines',true,...
            'LineAlpha',lineNoiseIn.p,...
            'Bandwidth',lineNoiseIn.fScanBandWidth,...
            'SignalType','Channels',...
            'ChanCompIndices',1:filteredEEG.nbchan,...
            'SlidingWinLength',lineNoiseIn.taperWindowSize,...
            'SlidingWinStep',lineNoiseIn.taperWindowStep,...
            'SmoothingFactor',lineNoiseIn.tau,...
            'PaddingFactor',lineNoiseIn.pad,...
            'ComputeSpectralPower',true,...
            'NormalizeSpectrum',false, ...
            'PlotFigures',false,...
            'VerboseOutput',0);
        lineNoiseOut = rmfield(lineNoiseOut,'EEG');
        
        % Add information to the lineNoiseOut struct:
        exact_eval_freqs = zeros(filteredEEG_noLN.nbchan, length(lineNoiseOut.linefreqs)); % stored in Hz
        avg_noise_red = zeros(filteredEEG_noLN.nbchan, length(lineNoiseOut.linefreqs)); %stored in dB
        for chan = 1:filteredEEG_noLN.nbchan
            for i_freq=1:length(lineNoiseOut.linefreqs)
                [~, f_ind] = min(abs(eval_freqs-lineNoiseOut.linefreqs(i_freq)));
                exact_eval_freqs(chan, i_freq) = eval_freqs(f_ind);
                avg_noise_red(chan, i_freq) = Spect_orig(chan,f_ind)-Spect_clean(chan,f_ind);
            end
        end
        
        lineNoiseOut.exact_eval_freqs = exact_eval_freqs;
        lineNoiseOut.avg_noise_red = avg_noise_red;
    case 'cleanLinePREP'
        [filteredEEG_noLN, lineNoiseOut] = cleanLineNoise(filteredEEG, lineNoiseIn);
    otherwise
        error('Not a valid Line noise removal method');
end

if plot
    figure
    subplot(2,1,1)
    spectopo(filteredEEG.data, filteredEEG.pnts, filteredEEG.srate, 'freq', 50,...
        'chanlocs', filteredEEG.chanlocs, 'freqrange', [0.1,60], 'percent', 50);
    subplot(2,1,2)
    spectopo(filteredEEG_noLN.data, filteredEEG_noLN.pnts, filteredEEG_noLN.srate, 'freq', 50,...
        'chanlocs', filteredEEG_noLN.chanlocs, 'freqrange', [0.1,60], 'percent', 50);
end
end

