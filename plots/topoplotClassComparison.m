function topoplotClassComparison(PSD, trials_sel, trials_cat, chanlocs, plot_params)

%Clims = 10*log10([min(PSD(:,plot_params.freqs2plot,trials_sel), [], 'all'),...
%     max(PSD(:,plot_params.freqs2plot,trials_sel), [], 'all')]/2);
%Clims = [min(mean(10*log10(PSD(:,plot_params.freqs2plot,trials_sel)),[2,3]), [], 'all'),...
%     max(mean(10*log10(PSD(:,plot_params.freqs2plot,trials_sel)),[2,3]), [], 'all')];
Clims = [-20,5];

for ch = 1:numel(chanlocs)
    if ~any(strcmp(plot_params.chans2plot, chanlocs(ch).labels))
    chanlocs(ch).labels = '.';
    end
end

%% Custom color map
N_colors_standard = 512;
myCmap = asymColorMapWhiteZero(Clims, N_colors_standard);
%set(0,'DefaultFigureColormap',myCmap)

% Without conv option
figure
subplot(1,2,1)
%topoplot(10*log10(mean(PSD(:,plot_params.freqs2plot,trials_sel&(~trials_cat)),[2,3])),...
%    chanlocs, 'maplimits', Clims, 'colormap', myCmap, 'electrodes', 'labels', 'conv', 'off');
topoplot(mean(10*log10(PSD(:,plot_params.freqs2plot,trials_sel&(~trials_cat))),[2,3]),...
    chanlocs, 'maplimits', Clims, 'colormap', myCmap, 'electrodes', 'labels', 'conv', 'off');
title(sprintf('%s', plot_params.classes{1}))
subplot(1,2,2)
%topoplot(10*log10(mean(PSD(:,plot_params.freqs2plot,trials_sel&(trials_cat)),[2,3])),...
%    chanlocs,'maplimits', Clims, 'colormap', myCmap, 'electrodes', 'labels', 'conv', 'off');
topoplot(mean(10*log10(PSD(:,plot_params.freqs2plot,trials_sel&(trials_cat))),[2,3]),...
    chanlocs,'maplimits', Clims, 'colormap', myCmap, 'electrodes', 'labels', 'conv', 'off');
title(sprintf('%s', plot_params.classes{2}))
cb = colorbar('southoutside');
cb.Position = [0.1,0.15,0.8,0.025];
cb.Label.String = 'Mean power amplitude (dB)';
set(cb, 'FontSize', 12)
if length(plot_params.freqs2plot)==1
    freq_str = sprintf('at %d Hz', plot_params.freqs2plot(1));
else
    freq_str = sprintf('between %d and %d Hz', plot_params.freqs2plot(1),plot_params.freqs2plot(end));
end
suptitle(sprintf('%s-%s %s', plot_params.name, plot_params.phase, freq_str))

%% Difference figure
Clims = [-2,2];
N_colors_standard = 512;
myCmap = asymColorMapWhiteZero(Clims, N_colors_standard);

figure
topoplot(mean(10*log10(PSD(:,plot_params.freqs2plot,trials_sel & trials_cat)),[2,3])-...
    mean(10*log10(PSD(:,plot_params.freqs2plot,trials_sel & ~trials_cat)),[2,3]),...
    chanlocs, 'maplimits', Clims, 'colormap', myCmap, 'electrodes', 'labels', 'conv', 'off');
title(sprintf('%s - %s', plot_params.classes{2}, plot_params.classes{1}))
cb = colorbar('southoutside');
cb.Position = [0.1,0.15,0.8,0.025];
cb.Label.String = 'Mean power amplitude (dB)';
set(cb, 'FontSize', 12)
if length(plot_params.freqs2plot)==1
    freq_str = sprintf('at %d Hz', plot_params.freqs2plot(1));
else
    freq_str = sprintf('between %d and %d Hz', plot_params.freqs2plot(1),plot_params.freqs2plot(end));
end
suptitle(sprintf('%s-%s %s', plot_params.name, plot_params.phase, freq_str))


% % With conv option
% figure
% subplot(1,2,1)
% p1 = topoplot(mean(Power_all_EO(:,all_f,Labels_all_EO==3),3), chanlocs_custom,...
%     'maplimits', Clims, 'electrodes', 'labels', 'conv', 'on');
% title('Stimulation - No Object')
% subplot(1,2,2)
% topoplot(mean(Power_all_EO(:,all_f,Labels_all_EO==4),3), EEG.chanlocs,...
%     'maplimits', Clims, 'electrodes', 'labels', 'conv', 'on');
% title('Stimulation - Object')
% cb =colorbar('AxisLocation', 'out', 'Position', [0.91,0.15,0.02,0.7]);
% cb.Label.String = 'Mean power amplitude (dB)';
% set(cb, 'FontSize', 12)
% suptitle('Mean (across trials) of power amplitude at 14 Hz averaged over the Eyes Open period')

%set(0,'DefaultFigureColormap',parula)
end