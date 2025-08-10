function [Features, chan_labels, dim2_labels] = select_features(Data, chanlocs, FoIs, cfg)
nb_trials = size(Data,3);
switch cfg.type
    case 'psd'
        if ischar(cfg.chans) && ischar(cfg.freqs)
            % 'all' channels
            nb_chans = size(Data,1);
            chan_labels = {chanlocs.labels};
            % 'all' frequencies
            nb_freqs = size(Data,2);
            dim2_labels = FoIs;
            Features = reshape(Data, [nb_chans*nb_freqs,nb_trials]);
        elseif ischar(cfg.chans)
            % 'all' channels
            nb_chans = size(Data,1);
            chan_labels = {chanlocs.labels};
            
            % frequencies are indicated with a range or a vector
            if length(cfg.freqs) == 2
                % input as limits
                dim2_labels = FoIs(FoIs>=cfg.freqs(1) & FoIs<=cfg.freqs(2));
                Data_sel = Data(:,FoIs>=cfg.freqs(1) & FoIs<=cfg.freqs(2),:);
            else
                % input as frequency vector
                dim2_labels = cfg.freqs;
                freq_inds = any(FoIs == dim2_labels',1);
                Data_sel = Data(:,freq_inds,:);
            end
            
            nb_freqs = numel(dim2_labels);
            Features = reshape(Data_sel, [nb_chans*nb_freqs,nb_trials]);
        else
            % cell of channels
            chan_labels = {};
            Data_sel = [];
            for ch = 1:numel(cfg.chans)
                ind = strcmp({chanlocs.labels},cfg.chans{ch});
                if any(ind)
                    Data_sel = cat(1, Data_sel, Data(ind,:,:));
                    chan_labels = [chan_labels, cfg.chans{ch}];
                end
            end
            nb_chans = numel(chan_labels);
            
            if ischar(cfg.freqs)
                % 'all' frequencies
                nb_freqs = size(Data,2);
                dim2_labels = FoIs;
                Features = reshape(Data_sel, [nb_chans*nb_freqs,nb_trials]);
            else
                % frequencies are indicated with a range or a vector
                if length(cfg.freqs) == 2
                    % input as limits
                    dim2_labels = FoIs(FoIs>=cfg.freqs(1) & FoIs<=cfg.freqs(2));
                    Data_sel = Data_sel(:,FoIs>=cfg.freqs(1) & FoIs<=cfg.freqs(2),:);
                else
                    % input as frequency vector
                    dim2_labels = cfg.freqs;
                    freq_inds = any(FoIs == dim2_labels',1);
                    Data_sel = Data_sel(:,freq_inds,:);
                end
                nb_freqs = numel(dim2_labels);
                Features = reshape(Data_sel, [nb_chans*nb_freqs,nb_trials]);
            end
        end
    case 'amp'
        if ischar(cfg.chans) && ischar(cfg.times)
            % 'all' channels
            nb_chans = size(Data,1);
            chan_labels = {chanlocs.labels};
            % 'all' times
            nb_times = size(Data,2);
            dim2_labels = 1:nb_times;
            Features = reshape(Data, [nb_chans*nb_times,nb_trials]);
        elseif ischar(cfg.chans)
            error('Not coded yet');
        else
            % cell of channels
            chan_labels = {};
            Data_sel = [];
            for ch = 1:numel(cfg.chans)
                ind = strcmp({chanlocs.labels},cfg.chans{ch});
                if any(ind)
                    Data_sel = cat(1, Data_sel, Data(ind,:,:));
                    chan_labels = [chan_labels, cfg.chans{ch}];
                end
            end
            nb_chans = numel(chan_labels);
            
            if ischar(cfg.times)
                % 'all' frequencies
                nb_times = size(Data,2);
                dim2_labels = 1:nb_times;
                Features = reshape(Data_sel, [nb_chans*nb_times,nb_trials]);
            else
                error('Not coded yet');
            end
        end
end
end