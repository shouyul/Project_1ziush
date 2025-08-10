function EEG = fillNaNs(EEG, cfg)
% Fill NaN periods when they are sufficiently short (default in sampling)

subject = cfg.subjects(cfg.current_subject).id;
[EEG, bound_evts] = findBoundaries(EEG, cfg);
events = EEG.event;

replaced_samples = [];
if ~isempty(bound_evts)
    for b = 1:length(bound_evts)+1
        if b == 1
            range = 1:events(bound_evts(b)).latency-1;
        elseif b == length(bound_evts)+1
            range = round(events(bound_evts(b-1)).latency:length(EEG.times));
        else
            range = round(events(bound_evts(b-1)).latency:events(bound_evts(b)).latency-1);
        end
        data2inspect = EEG.data(:, range);
        NanDetection = isnan(data2inspect);
        % Find time samples where at least one channel contains a NaN:
        NanDetection = sum(NanDetection,1)>0;
        
        s=1;
        last_valid_sample = 0;
        while s <= length(NanDetection)
            if NanDetection(s)
                s = s+1;
            else
                if last_valid_sample>0 && s-last_valid_sample>1 && s-last_valid_sample<=cfg.maxNans2Replace+1
                    NaNChannels = isnan(data2inspect(:,last_valid_sample:s));
                    % Find channels where at least one time sample contains a NaN:
                    NaNChannels = sum(NaNChannels,2)>0;
                    
                    newdata = interp1([last_valid_sample, s],...
                        [data2inspect(NaNChannels,last_valid_sample),data2inspect(NaNChannels,s)]', last_valid_sample:s);
                    data2inspect(NaNChannels,last_valid_sample:s) = newdata';
                    
                    replaced_samples = [replaced_samples, range(1)-1+(last_valid_sample+1:s-1)];
                end
                last_valid_sample = s;
                s = s+1;
            end
        end
        
        EEG.data(:, range) = data2inspect;
    end
else
    range = 1:EEG.pnts;
    data2inspect = EEG.data(:, range);
    NanDetection = isnan(data2inspect);
    % Find time samples where at least one channel contains a NaN:
    NanDetection = sum(NanDetection,1)>0;
    
    s=1;
    last_valid_sample = 0;
    while s <= length(NanDetection)
        if NanDetection(s)
            s = s+1;
        else
            if last_valid_sample>0 && s-last_valid_sample>1 && s-last_valid_sample<=cfg.maxNans2Replace+1
                NaNChannels = isnan(data2inspect(:,last_valid_sample:s));
                % Find channels where at least one time sample contains a NaN:
                NaNChannels = sum(NaNChannels,2)>0;
                
                newdata = interp1([last_valid_sample, s],...
                    [data2inspect(NaNChannels,last_valid_sample),data2inspect(NaNChannels,s)]', last_valid_sample:s);
                data2inspect(NaNChannels,last_valid_sample:s) = newdata';
                
                replaced_samples = [replaced_samples, range(1)-1+(last_valid_sample+1:s-1)];
            end
            last_valid_sample = s;
            s = s+1;
        end
    end
    
    EEG.data(:, range) = data2inspect;
end

fprintf('%s: Filled %d NaNs by interpolation (%.1f%%)\n', subject,...
    length(replaced_samples), 100*length(replaced_samples)/length(EEG.times))
EEG.etc.NaNsFilled = replaced_samples;
end