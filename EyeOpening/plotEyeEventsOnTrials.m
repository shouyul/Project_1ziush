function plotEyeEventsOnTrials(EEG1, EEG2, Peaks, options)
% All trials at once
trials = double(unique(Peaks.Trial)');

meanData1 = squeeze(mean(EEG1.data(:,:,:),1))/1000; % convert to mV
meanData2 = squeeze(mean(EEG2.data(:,:,:),1))/1000; % convert to mV
decMeanData1 = zeros(size(meanData1));
decMeanData2 = zeros(size(meanData2));

for tr = trials
    decMeanData1(:,tr) = meanData1(:,tr) - options.step*(tr-1);
    decMeanData2(:,tr) = meanData2(:,tr) - options.step*(tr-1);
end

peaks4plot = nan(size(Peaks,1),1);
for p = 1:size(Peaks,1)
    tr = Peaks.Trial(p);
    peaks4plot(p) = decMeanData2(Peaks.Latency(p),tr);
end

% eyeOp_peaks = Peaks.EyeOpening;
% eyeCl_peaks = Peaks.EyeClosing;
% blink_peaks = Peaks.Blink;
% rej_peaks = Peaks.Rejected;
eyeOp_peaks = strcmp(Peaks.AutoCategory, 'EyeOpening');
eyeCl_peaks = strcmp(Peaks.AutoCategory, 'EyeClosing');
rej_peaks = strcmp(Peaks.AutoCategory, 'Rejected');

f = figure;hold on;
f.WindowState = 'maximized';
plot(EEG1.times,decMeanData1,':k');
plot(EEG2.times,decMeanData2);
scatter(EEG2.times(Peaks.Latency(eyeOp_peaks)), peaks4plot(eyeOp_peaks),'g','o')
scatter(EEG2.times(Peaks.Latency(eyeCl_peaks)), peaks4plot(eyeCl_peaks),'g','x')
%scatter(EEG2.times(Peaks.Latency(blink_peaks)), peaks4plot(blink_peaks),'g','v')
scatter(EEG2.times(Peaks.Latency(rej_peaks)), peaks4plot(rej_peaks),'r','*')

xline(options.timeEO, 'Label', 'Eye Opening event');
xline(options.timeEC, 'Label', 'Eye Closing event');
for tr = trials
    yl = yline(-options.step*(tr-1), 'Label',sprintf('Trial%d',tr));
    yl.LabelHorizontalAlignment = 'center';
    yl.LabelVerticalAlignment = 'bottom';
end
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
title(sprintf('%s - %s', options.subject, EEG2.chanlocs(:).labels))
end