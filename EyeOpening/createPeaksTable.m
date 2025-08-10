function Peaks = createPeaksTable(EEG1, EEG2, trials, options)
% General parameters
latsLeft1 = -round(options.wdw_diff1*EEG1.srate):-1;
latsRight1 = 1:round(options.wdw_diff1*EEG1.srate);
latsLeft2 = -round(options.wdw_diff2*EEG1.srate):-1;
latsRight2 = 1:round(options.wdw_diff2*EEG1.srate);

varNames = {'Trial','Quantile','Latency','Max1','Ave1',...
    'Height1','Width1','Prominence1','Max2','Ave2','Height2','Diff1','Diff2',...
    'EyeOpening','EyeClosing','Blink','Rejected'};
varTypes = {'int32','double','int32','double','double',...
    'double','double','double','double','double','double','double','double',...
    'logical','logical','logical','logical'};

Peaks = table('Size',[0,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
i = 0;
for tr = trials
    EEGave1 = mean(EEG1.data(:,:,tr), 'all');
    EEGmax1 = max(EEG1.data(:,:,tr), [], 'all');
    EEGquant = quantile(EEG1.data(:,:,tr), options.quant/100, 'all');
    [heights, locs, widths, proms] = findpeaks(squeeze(EEG1.data(:,:,tr)),...
        'MinPeakHeight', EEGquant, 'MinPeakProminence', EEGmax1*options.percMax4Prom);
    
    EEGave2 = mean(EEG2.data(:,:,tr), 'all');
    EEGmax2 = max(EEG2.data(:,:,tr), [], 'all');
    
    Npeaks = length(locs);
    for p = 1:Npeaks
        i = i+1;
        Peaks.Trial(i) = tr;
        Peaks.Quantile(i) = options.quant;
        Peaks.Latency(i) = locs(p);
        Peaks.Max1(i) = EEGmax1;
        Peaks.Ave1(i) = EEGave1;
        Peaks.Height1(i) = heights(p);
        Peaks.Width1(i) = widths(p);
        Peaks.Prominence1(i) = proms(p);
        Peaks.Max2(i) = EEGmax2;
        Peaks.Ave2(i) = EEGave2;
        
        % Default values
        Peaks.EyeOpening(i) = false;
        Peaks.EyeClosing(i) = false;
        Peaks.Blink(i) = false;
        Peaks.Rejected(i) = false;
        if locs(p) > options.wdw_diff2*EEG1.srate && locs(p) < EEG1.pnts - options.wdw_diff2*EEG1.srate
            Peaks.Height2(i) = mean(EEG2.data(:,locs(p),tr),'all');
            
            meanLeft1 = mean(EEG2.data(:,locs(p)+latsLeft1,tr),'all');
            meanRight1 = mean(EEG2.data(:,locs(p)+latsRight1,tr),'all');
            diffLR1 = meanLeft1-meanRight1;
            Peaks.Diff1(i) = diffLR1;
            
            meanLeft2 = mean(EEG2.data(:,locs(p)+latsLeft2,tr),'all');
            meanRight2 = mean(EEG2.data(:,locs(p)+latsRight2,tr),'all');
            diffLR2 = meanLeft2-meanRight2;
            Peaks.Diff2(i) = diffLR2;
        else
            Peaks.Height2(i) = NaN;
            Peaks.Diff1(i) = NaN;
            Peaks.Diff2(i) = NaN;
            Peaks.Rejected(i) = true;
        end
    end
    
    if options.userLabel
        f = figure;hold on;
        f.WindowState = 'maximized';
        plot(EEG1.times,mean(EEG1.data(:,:,tr),1),':k');
        plot(EEG2.times,mean(EEG2.data(:,:,tr),1));
        xline(options.timeEO, 'Label', 'Eye Opening event');
        xline(options.timeEC, 'Label', 'Eye Closing event');
        yline(0,'--');
        yline(EEGquant,':');
        xlabel('Time (ms)');
        ylabel('Amplitude (microV)');
        title(sprintf('Trial %d', tr));
        
        peaks_tr = find(Peaks.Trial == tr & ~Peaks.Rejected)';
        scatter(EEG1.times(Peaks.Latency(peaks_tr)), Peaks.Height1(peaks_tr),'k','*');
        
        [selection, labels] = selectPeaks(Peaks, peaks_tr, EEG2.times);
        Peaks.Rejected(setdiff(peaks_tr,selection)) = true;
        for s = 1:length(selection)
            switch labels{s}
                case 'EyeOpening'
                    Peaks.EyeOpening(selection(s)) = true;
                case 'EyeClosing'
                    Peaks.EyeClosing(selection(s)) = true;
                case 'Blink'
                    Peaks.Blink(selection(s)) = true;
            end
        end
        %         for p = peaks_tr
        %             if ~isnan(Peaks.Height2(p))
        %                 pt = scatter(EEG2.times(Peaks.Latency(p)), Peaks.Height2(p),'m','*');
        %                 accepted = false;
        %                 while ~accepted
        %                     dec = input('Which category? [0=Exclude; 1=EyeOpening; 2=EyeClosing; 3=Blink] ');
        %                     if ~isempty(dec)
        %                         switch dec
        %                             case 0
        %                                 accepted = true;
        %                                 Peaks.Rejected(p) = true;
        %                                 delete(pt);
        %                                 scatter(EEG2.times(Peaks.Latency(p)), Peaks.Height2(p),'r','+');
        %                             case 1
        %                                 accepted = true;
        %                                 Peaks.EyeOpening(p) = true;
        %                                 delete(pt);
        %                                 scatter(EEG2.times(Peaks.Latency(p)), Peaks.Height2(p),'g','o');
        %                             case 2
        %                                 accepted = true;
        %                                 Peaks.EyeClosing(p) = true;
        %                                 delete(pt);
        %                                 scatter(EEG2.times(Peaks.Latency(p)), Peaks.Height2(p),'g','x');
        %                             case 3
        %                                 accepted = true;
        %                                 Peaks.Blink(p) = true;
        %                                 delete(pt);
        %                                 scatter(EEG2.times(Peaks.Latency(p)), Peaks.Height2(p),'g','v');
        %                             otherwise
        %                                 fprintf('Decision not understood.\n')
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        close gcf
    end
end

    function [selection, labels] = selectPeaks(Peaks, peaks_tr, times)
        peak_coords = [times(Peaks.Latency(peaks_tr))',Peaks.Height1(peaks_tr)];
        
        % give the user instructions
        fprintf('Select Peaks corresponding to EyeOpening, EyeClosing and Blink events...\n')
        disp('First use the mouse to click on a peak');
        disp('Then press "o" to label it EyeOpening, "c" to label it EyeClosing or "b" to label it Blink.');
        disp('Alternatively, press "r" to remove the last peak selected.');
        disp('Eventually, press "q" to finish selection');
        
        marker = '*'; markersize = 36; markercolor = 'm';
        modelAxes = get(gcf, 'Children');
        set(modelAxes,'NextPlot','add');
        
        done = false;
        
        selection = zeros(0,1);
        labels = cell(0,1);
        
        listen=true; % Listen boolean for mouse clicks
        while ~done
            k = waitforbuttonpress;
            if listen
                pk = get(modelAxes, 'currentpoint');
                pk = pk(1,1:2);
            end
            if k == 1 %checks if waitforbuttonpress was a key
                key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
                if strcmp(key, 'q')
                    % finished selecting points
                    done = true;
                    fprintf('Peak selection finished.\n')
                    continue
                end
                
                if ~listen
                    if strcmp(key, 'o')
                        labels(end+1,1) = {'EyeOpening'};
                        last_plot = scatter(times(Peaks.Latency(selection(end))), Peaks.Height2(selection(end)),...
                            markersize, 'g', 'o');
                    elseif strcmp(key, 'c')
                        labels(end+1,1) = {'EyeClosing'};
                        last_plot = scatter(times(Peaks.Latency(selection(end))), Peaks.Height2(selection(end)),...
                            markersize, 'g', 'x');
                    elseif strcmp(key, 'b')
                        labels(end+1,1) = {'Blink'};
                        last_plot = scatter(times(Peaks.Latency(selection(end))), Peaks.Height2(selection(end)),...
                            markersize, 'g', 'v');
                    elseif strcmp(key, 'r')
                        % remove last point
                        delete(last_plot)
                        if length(selection) == 1
                            selection = zeros(0,1);
                        else
                            selection = selection(1:end-1);
                        end
                    else
                        % Not a valid command
                        continue
                    end
                    listen = true;
                end
            else
                % it was a mouse click
                if listen
                    if isempty(pk)
                        fprintf('Click not registered.\n')
                    else
                        % Find the peak closest to the selected point
                        [~, pk_sel] = min(vecnorm(peak_coords-pk,2,2));
                        % Plot the selected point
                        if any(selection == peaks_tr(pk_sel))
                            disp('Peak selected twice, consider removing it.')
                            last_plot = scatter(peak_coords(pk_sel,1), peak_coords(pk_sel,2),...
                                markersize, 'r', marker);
                        else
                            last_plot = scatter(peak_coords(pk_sel,1), peak_coords(pk_sel,2),...
                                markersize, markercolor, marker);
                        end
                        
                        selection(end+1,1) = peaks_tr(pk_sel);
                        listen=false; % stop listening
                    end
                end
            end
        end
    end
end