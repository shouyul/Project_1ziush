function Peaks = createPeaksTable2(EEG1, EEG2, options)
% General parameters
latsLeft1 = -round(options.wdw_diff1*EEG1.srate):-1;
latsRight1 = 1:round(options.wdw_diff1*EEG1.srate);
latsLeft2 = -round(options.wdw_diff2*EEG1.srate):-1;
latsRight2 = 1:round(options.wdw_diff2*EEG1.srate);

varNames = {'Trial','Quantile','Latency','Max1','Ave1',...
    'Height1','Width1','Prominence1','Max2','Ave2','Height2','Diff1','Diff2','DiffIndex',...
    'AutoCategory','ManualCategory'};
varTypes = {'int32','double','int32','double','double',...
    'double','double','double','double','double','double','double','double','double',...
    'string','string'};

Peaks = table('Size',[0,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
i = 0;
for tr = 1:EEG1.trials
    EEGave1 = mean(EEG1.data(:,:,tr), 'all');
    EEGmax1 = max(EEG1.data(:,:,tr), [], 'all');
    EEGquant = quantile(EEG1.data(:,:,tr), options.quant/100, 'all');
    [heights, locs, widths, proms] = findpeaks(squeeze(EEG1.data(:,:,tr)),...
        'MinPeakHeight', EEGquant, 'MinPeakProminence', EEGmax1*options.percMax4Prom);
    
    EEGave2 = mean(EEG2.data(:,:,tr), 'all');
    EEGmax2 = max(EEG2.data(:,:,tr), [], 'all');
    
    if EEGmax2 > options.artifact_detection_th
        % Skip this trial if suspicion of artifacts
        fprintf('Trial %d probably contains large artifacts, skipping peak inspection...\n',tr)
        continue
    end
    
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
        Peaks.AutoCategory{i} = 'None';
        Peaks.ManualCategory{i} = 'None';
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
            
            Peaks.DiffIndex(i) = (diffLR2 - diffLR1)/abs(diffLR1);
        else
            Peaks.Height2(i) = NaN;
            Peaks.Diff1(i) = NaN;
            Peaks.Diff2(i) = NaN;
            Peaks.DiffIndex(i) = NaN;
        end
    end
end

% Remove peaks where diff was not calculated
Peaks = Peaks(~isnan(Peaks.Diff1),:);

% Reject peaks that do not satisfy important diff properties:
% 1. Diff1 and Diff2 should have the same sign
Peaks = Peaks(Peaks.Diff1.*Peaks.Diff2 > 0,:);

% 2.a Select the most positive differences
pos_peaks = Peaks.Diff1 > 0;
pos_quant = quantile(Peaks.Diff1(pos_peaks), options.diffQuant/100); %threshold for positive difference
% 2.b Select the most negative differences
neg_peaks = Peaks.Diff1 < 0;
neg_quant = -quantile(abs(Peaks.Diff1(neg_peaks)), options.diffQuant/100); %threshold for negative difference
Peaks = Peaks(Peaks.Diff1 > pos_quant | Peaks.Diff1 < neg_quant,:);

% 3. DiffIndex should have the same sign as Diff1
%Peaks = Peaks(Peaks.DiffIndex.*Peaks.Diff1 > 0,:);

for tr = 1:EEG1.trials
    peaks_tr = find(Peaks.Trial == tr)';
    
    if ~isempty(peaks_tr)
        % Automatic assignment with max(abs(diff)) for each category
        [vals, order] = sort(Peaks.Diff1(peaks_tr),'ascend');
        
        for i = length(peaks_tr):-1:1
            t = EEG1.times(Peaks.Latency(peaks_tr(order(i))));
            if vals(i) > 0
                if t > options.timeEO + options.accept_wdw(1)*1000 &&...
                        t < options.timeEO + options.accept_wdw(2)*1000 % conversion to ms
                    bestEO = order(i);
                    Peaks.AutoCategory{peaks_tr(bestEO)} = 'EyeOpening';
                    break
                else
                    continue
                end
            else
                bestEO = [];
                break
            end
        end
        
        for i = 1:length(peaks_tr)
            t = EEG1.times(Peaks.Latency(peaks_tr(order(i))));
            if vals(i) < 0
                if t > options.timeEC + options.accept_wdw(1)*1000 &&...
                        t < options.timeEC + options.accept_wdw(2)*1000 % conversion to ms
                    bestEC = order(i);
                    Peaks.AutoCategory{peaks_tr(bestEC)} = 'EyeClosing';
                    break
                else
                    continue
                end
            else
                bestEC = [];
                break
            end
        end
        
%         [~, bestEO] = max(Peaks.Diff1(peaks_tr));
%         if Peaks.Diff1(peaks_tr(bestEO)) > 0
%             Peaks.AutoCategory{peaks_tr(bestEO)} = 'EyeOpening';
%         else
%             bestEO = [];
%         end
%         
%         [~, bestEC] = min(Peaks.Diff1(peaks_tr));
%         if Peaks.Diff1(peaks_tr(bestEC)) < 0
%             Peaks.AutoCategory{peaks_tr(bestEC)} = 'EyeClosing';
%         else
%             bestEC = [];
%         end
        
        Peaks.AutoCategory(peaks_tr(setdiff(1:length(peaks_tr),[bestEO,bestEC]))) = 'Rejected';
    end
    
    if options.userLabel && any(options.inspectedTrials == tr)
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
        
        if ~isempty(peaks_tr)
        if ~isempty(bestEO)
            scatter(EEG1.times(Peaks.Latency(peaks_tr(bestEO))), Peaks.Height2(peaks_tr(bestEO)),'m','o');
        end
        if ~isempty(bestEC)
            scatter(EEG1.times(Peaks.Latency(peaks_tr(bestEC))), Peaks.Height2(peaks_tr(bestEC)),'m','x');
        end
        scatter(EEG1.times(Peaks.Latency(peaks_tr)), Peaks.Height1(peaks_tr),'k','*');
        end
        
        [selection, labels] = selectPeaks(Peaks, peaks_tr, EEG2.times);
        if isempty(selection)
            % Automatic selection was accepted
            Peaks.ManualCategory(peaks_tr) = Peaks.AutoCategory(peaks_tr);
        else
            % Manual selection by user
            Peaks.ManualCategory(setdiff(peaks_tr,selection)) = 'Rejected';
            for s = 1:length(selection)
                Peaks.ManualCategory(selection(s)) = labels{s};
            end
        end
        close gcf        
    end
end

    function [selection, labels] = selectPeaks(Peaks, peaks_tr, times)
        peak_coords = [times(Peaks.Latency(peaks_tr))',Peaks.Height1(peaks_tr)];
        
        % give the user instructions
        fprintf('Select Peaks corresponding to EyeOpening and EyeClosing events...\n')
        disp('If automatic selection is good, press "a"');
        disp('First use the mouse to click on a peak');
        disp('Then press "o" to label it EyeOpening or "c" to label it EyeClosing.');
        disp('Alternatively, press "r" to remove the last peak selected.');
        disp('Eventually, press "s" to finish selection');
        
        marker = '*'; markersize = 36; markercolor = 'm';
        modelAxes = get(gcf, 'Children');
        set(modelAxes,'NextPlot','add');
        
        done = false;
        
        selection = zeros(0,1);
        labels = cell(0,1);
        % supAxes = axes('pos',[0 0.95 1 1],'visible','off');
        % supText = text(supAxes,.5,0,['Select ' char(chanLabels(1))],...
        %     'FontSize',get(gcf,'defaultaxesfontsize')+4,...
        %     'horizontalalignment','center');
        
        listen=true; % Listen boolean for mouse clicks
        while ~done
            k = waitforbuttonpress;
            if listen
                pk = get(modelAxes, 'currentpoint');
                pk = pk(1,1:2);
            end
            if k == 1 %checks if waitforbuttonpress was a key
                key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
                if strcmp(key, 'a') && isempty(selection)
                    done = true;
                    fprintf('Automatic Peak selection accepted.\n')
                    continue
                elseif strcmp(key, 's') && ~isempty(selection)
                    % finished selecting points
                    done = true;
                    fprintf('Manual Peak selection finished.\n')
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
                    elseif strcmp(key, 'r')
                        % remove last point
                        delete(last_plot)
                        if length(selection) == 1
                            selection = zeros(0,1);
                        else
                            selection = selection(1:end-1);
                        end
                        
                        %             if ~isempty(selected)
                        %                 if ~isempty(marker)
                        %                     delete(findobj('marker', '*'));
                        %                     hs = plot3(selected(1:end-1,1), selected(1:end-1,2), selected(1:end-1,3), [markercolor marker]);
                        %                     set(hs, 'MarkerSize', markersize);
                        %                 end
                        %                 fprintf('Removed %s at [%9.4f %9.4f %9.4f]\n', char(chanLabels(size(selected,1))),...
                        %                     selected(end,1), selected(end,2), selected(end,3));
                        %                 selected = selected(1:end-1,:);
                        %                 set(supText,'String', ['Select ', char(chanLabels(size(selected,1)+1))]);
                        %             end
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
                        %set(supText,'String', ['Select ', char(chanLabels(size(selected,1)+1))]);
                        %fprintf('Selected %s at [%9.4f %9.4f %9.4f] \n', char(chanLabels(size(selected,1))),...
                        %    selected(end,1), selected(end,2), selected(end,3));
                        listen=false; % stop listening
                    end
                end
            end
        end
    end
end