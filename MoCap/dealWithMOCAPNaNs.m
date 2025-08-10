function EEG_mocap = dealWithMOCAPNaNs(EEG_mocap, rbname)

rb_chans = contains({EEG_mocap.chanlocs.labels}, rbname);
NaNs_rb = sum(isnan(EEG_mocap.data(rb_chans,:)),1)>0;
NaN_intervals = mask2intervals(NaNs_rb);

bnd_evts = find(strcmp({EEG_mocap.event.type}, 'boundary'));
bnd_lats = [EEG_mocap.event(bnd_evts).latency];

if any(NaNs_rb)
    warning('Still NaNs in the MoCap data set')
    switch rbname
        case 'RB1'
            disp('Correcting RB1 NaNs...')
            
            for i = 1:size(NaN_intervals,1)
                EEG_mocap = interpNaNs(EEG_mocap, rb_chans, NaN_intervals(i,:), bnd_lats);
            end
        case 'RB2'
            % NaNs in RB2 are normal in the GogglesOFF condition
            disp('Correcting RB2 NaNs...')
            for i = 1:size(NaN_intervals,1)
                EEG_mocap = interpNaNs(EEG_mocap, rb_chans, NaN_intervals(i,:), bnd_lats);
                
                %                 b = find(NaN_intervals(i,1) >= bnd_lats, 1, 'last');
                %                 if isempty(b)
                %                     % Probably the beginning of the recording
                %                     error('Not coded yet')
                %                 else
                %                     cond = EEG_mocap.event(bnd_evts(b)+1).condition;
                %                 end
                %                 switch cond
                %                     case 'GogglesOFF'
                %                         % Simply replace by zeros
                %                         EEG_mocap.data(rb_chans,NaN_intervals(i,1):NaN_intervals(i,2)) = 0;
                %                     case 'GogglesON'
                %                         % Same method as for RB1
                %                         len_interp_segment = (diff(NaN_intervals(i,:))+1)*10;
                %
                %                         if any(NaN_intervals(i,2) == floor(bnd_lats))
                %                             lat_start = NaN_intervals(i,2) - len_interp_segment;
                %                             x = lat_start:NaN_intervals(i,1)-1;
                %                             xq = lat_start:NaN_intervals(i,2);
                %                         elseif any(NaN_intervals(i,1) == ceil(bnd_lats))
                %                             lat_stop = NaN_intervals(i,1) + len_interp_segment;
                %                             x = NaN_intervals(i,2)+1:lat_stop;
                %                             xq = NaN_intervals(i,1):lat_stop;
                %                         else
                %                             error('NaN not touching any boundary');
                %                         end
                %                         EEG_mocap.data(rb_chans,xq) = interp1(x, EEG_mocap.data(rb_chans,x)', xq, 'linear', 'extrap')';
                %                 end
            end
            
        otherwise
    end
end

    function EEG_mocap = interpNaNs(EEG_mocap, rb_chans, interval, bnd_lats)
        len_interp_segment = (diff(interval)+1)*10; % Larger than actual NaN segment, *10 only works if NaN segment is not too big
        
        if any(interval(1) == [0, ceil(bnd_lats)]) && any(interval(2) == [floor(bnd_lats), EEG_mocap.pnts])
            % Full block of NaNs, simply replace by zeros
            EEG_mocap.data(rb_chans,interval(1):interval(2)) = 0;
        else
            if any(interval(1) == ceil(bnd_lats))
                lat_stop = interval(1) + len_interp_segment;
                x = interval(2)+1:lat_stop;
                xq = interval(1):lat_stop;
                EEG_mocap.data(rb_chans,xq) = interp1(x, EEG_mocap.data(rb_chans,x)', xq, 'linear', 'extrap')';
            elseif any(interval(2) == floor(bnd_lats))
                lat_start = interval(2) - len_interp_segment;
                x = lat_start:interval(1)-1;
                xq = lat_start:interval(2);
                EEG_mocap.data(rb_chans,xq) = interp1(x, EEG_mocap.data(rb_chans,x)', xq, 'linear', 'extrap')';
            else
                enclosed_bnd = bnd_lats > interval(1) & bnd_lats < interval(2);
                if sum(enclosed_bnd) == 1
                    lat_stop = ceil(bnd_lats(enclosed_bnd)) + len_interp_segment;
                    x = interval(2)+1:lat_stop;
                    xq = ceil(bnd_lats(enclosed_bnd)):lat_stop;
                    EEG_mocap.data(rb_chans,xq) = interp1(x, EEG_mocap.data(rb_chans,x)', xq, 'linear', 'extrap')';
                    lat_start = floor(bnd_lats(enclosed_bnd)) - len_interp_segment;
                    x = lat_start:interval(1)-1;
                    xq = lat_start:floor(bnd_lats(enclosed_bnd));
                    EEG_mocap.data(rb_chans,xq) = interp1(x, EEG_mocap.data(rb_chans,x)', xq, 'linear', 'extrap')';
                else
                    error('New case to deal with')
                end
            end
        end
    end
end