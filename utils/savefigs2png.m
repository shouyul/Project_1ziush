% path = ['C:\Users\Alexandre\Desktop\data\figures\IC_clustering\',...
% 'mainSTUDY_3_automatic\50clust_sig3_RSC4\Allclusters\user\Brain\PureClusters\'];
path = ['D:\Data_FourStreetsCity\figures\ERSPs\bemobil\APP\Brain_BrainWithNoise\depth95\'];
overwrite = false;
target_extension = 'png';

list_files = dir(path);
list_files = {list_files.name};
list_figures = list_files(contains(list_files, '.fig'));

for f = 1:length(list_figures)
    if isfile(sprintf('%s%s.%s', path, list_figures{f}(1:end-4), target_extension)) && ~overwrite
        continue
    else
%         if ~contains(list_figures{f}(1:end-4), 'OPA')
%             continue
%         end
        open([path, list_figures{f}])
        fig = gcf;
        fig.InvertHardcopy = 'off';
        set(fig, 'Position', [50, 50, 1250, 600]);
        ext = sprintf('-d%s', target_extension);
        print([path, list_figures{f}(1:end-4)], ext, '-r300');
        close(fig)
    end
end

% 
% ERSP_files_list = list_files(contains(list_files, 'ERSP'));
% 
% %% ave files
% ERSP_ave_files_list = ERSP_files_list(contains(ERSP_files_list, 'ave_'));
% %ERSP_ave_files_list = ERSP_ave_files_list(~contains(ERSP_ave_files_list, '_Alloave_'));
% ERSP_ave_files_list = ERSP_ave_files_list(~contains(ERSP_ave_files_list, '_3conds_'));
% ERSP_figs_list = ERSP_ave_files_list(contains(ERSP_ave_files_list, '.fig'));
% 
% for f = 1:length(ERSP_figs_list)
%     open([path, ERSP_figs_list{f}])
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     %saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'png');
%     saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'svg');
%     close(gcf)
% end
% 
% %% ave_3conds files
% ERSP_ave_files_list = ERSP_files_list(contains(ERSP_files_list, '_ave_3conds_'));
% ERSP_figs_list = ERSP_ave_files_list(contains(ERSP_ave_files_list, '.fig'));
% 
% for f = 1:length(ERSP_figs_list)
%     open([path, ERSP_figs_list{f}])
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'png');
%     saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'svg');
%     close(gcf)
% end
% 
% %% diff files
% ERSP_diff_files_list = ERSP_files_list(contains(ERSP_files_list, '_diff_'));
% ERSP_figs_list = ERSP_diff_files_list(contains(ERSP_diff_files_list, '.fig'));
% 
% for f = 1:length(ERSP_figs_list)
%     open([path, ERSP_figs_list{f}])
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     sdf('transparencyEnabled')
%     saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'png');
%     saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'svg');
%     close(gcf)
% end
% 
% %% allSubjs files
% % ERSP_allSubjs_files_list = ERSP_files_list(contains(ERSP_files_list, '_allSubjects_'));
% % ERSP_figs_list = ERSP_allSubjs_files_list(contains(ERSP_allSubjs_files_list, '.fig'));
% % 
% % for f = 1:length(ERSP_figs_list)
% %     open([path, ERSP_figs_list{f}])
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1])
% %     saveas(gcf, [path, ERSP_figs_list{f}(1:end-4)], 'png');
% %     close(gcf)
% % end