function [maxstat_val, cluster_vals, clusters] = clusteredStats(data, options_stats)
% Gives the clustered statistics for N conditions according to the options.
% 1. For every sample (ex: freq x time for ERSPs), the statistical test compare the N
%       conditions (ex: t-test if 2 conditions, ANOVA if 3 conditions).
%       !!!! N>3 has not been implemented here
% 2. Select all samples whose statistical value is larger than some
%       threshold (based on the test distribution).
% 3. Cluster the selected samples in connected sets on the basis of adjacency.
%       The notion of adjacency is straightforward for time and frequency
%       but needs some preparation when dealin with spatialized input. (see defineNeighbours function)
%       Possibility to exclude clusters that are too small
%       with respect to the size of the dataset (see options).
% 4. Calculate cluster-level statistics by taking the sum of the statistical values within a cluster.
% Follows the procedure detailed in:
% Maris, E. & Oostenveld, R. Nonparametric statistical testing of EEG- and MEG-data. J. Neurosci. Methods 164, 177?190 (2007).
%
% Inputs:
% data              - Data. Should be a 1 dimensional cell whose length corresponds to
%                   the number of conditions. For each 1 condition,
%                   the last dimension should correspond to trials (on
%                   which permutations will be computed). The first
%                   dimensions should be as indicated by the 'style' field in options.
% options_stats     - Options for this script. struct with fields:
%                   'fields' cell containing the names of the conditions.
%                   'model' string specificying the statistical model used.
%                       Supported: 'classic' or 'mixedEffects'
%                   'pairing' string indicating whether to compute paired
%                       statistical tests. Supported: 'on' or 'off'.
%                   'style' Specify the type of data with a string.
%                       Ex: 'timeXfreq' for ERSPs of a single channel/IC
%                       Supported: 'freq', 'timeXfreq', 'chanXfreq'
%                   'Chans' (optional) Only required if style is 'chanXfreq'.
%                       Channel indices corresponding to the input
%                       (following numerotation of the EEG.chanlocs field)
%                   'ElecFile' (optional) Only required if style is 'chanXfreq'.
%                       Path to the manufacturer electrode file (.elc)
%                   'MaxDeg' (optional) Only required if style is 'chanXfreq'.
%                       Solid angle in degree: Limit to define neighbouring channels.
%                   'removeSmallestClusters'. boolean to remove smallest
%                       clusters (see step 3 above).
%%%%% fields specific to 'mixedEffects' model:
%                   For each field specified in 'fields', a
%                       'trialsInfo_(fieldName)' table containing all the
%                       variables that are involved in the formula (see fitlme function)
%                   'formula' formula for the lme (see fitlme function)
%                   'MEterm' term of the formula that corresponds to the
%                       contrast between conditions.
%                   'original' boolean specifying whether the current
%                       data at test is the original distribution of trials (true)
%                       or a permuation (false)
%                   %%% The last three options are for saving purposes only
%                   (custom to 4SC experiment)
%                   'saveFolder' path to the saving folder
%                   'modelData' string specifying the data used for
%                       statistics (ex: 'test' for test trials only)
%                   'ROI' name of the ROI from which the data is extracted.
%                       ex: 'RSC', 'PPA', 'OPA'
%
% Outputs:
%   maxstat_val     - Largest of the cluster-level statistics (absolute value).
%   cluster_vals    - Computed statistical values for all clusters.
%   clusters        - vector or matrix containing the information
%                       about the clusters (each sample is
%                       associated with a cluster index - in correspondence
%                       with the cluster_vals vector - or 0 if it does not
%                       belong to any cluster). The size of this output
%                       will be the same as a cell of data, except for the
%                       last dimension

n_conds = length(options_stats.fields);

if n_conds < 2
    error('Cannot perform this analysis with only one condition')
end

if n_conds == 2
    alpha_crit = 0.975;
else
    alpha_crit = 0.95;
end

%% Perform the test comparing the samples:
switch options_stats.model
    case 'classic'
        stats = statcond(data, 'paired', options_stats.pairing,...
            'method', 'param', 'structoutput', 'on','verbose', 'off');
        values = stats.stat;
        
        % Define threshold
        if n_conds == 2
            if strcmp(options_stats.pairing, 'on')
                thresh = tinv(alpha_crit, stats.df);
            else
                thresh = tinv(alpha_crit, mean(stats.df));
            end
        else
            thresh = finv(alpha_crit, stats.df(1), stats.df(2));
        end
        
        switch options_stats.style
            case 'time'
                [clusters, cluster_count] = formClustersUniDim(values, thresh, n_conds);
            case 'timeXfreq'
                [clusters, cluster_count] = formClustersTimeFreq(values, thresh, n_conds);
            case 'freq'
                [clusters, cluster_count] = formClustersUniDim(values, thresh, n_conds);
            case 'chanXfreq'
                [clusters, cluster_count] = formClustersChanFreq(values, thresh, n_conds,...
                    options_stats.Chans, options_stats.ElecFile, options_stats.MaxDeg);
            otherwise
                error('Unknown style')
        end
        [maxstat_val, cluster_vals, clusters] = computeClusterStats(values, clusters, cluster_count,...
            options_stats.removeSmallestClusters);
        
        %         [maxstat_val, cluster_vals, clusters] = computeClusterStatsTimeFreq(values, thresh, n_conds,...
        %             options_stats.removeSmallestClusters);
        
    case 'mixedEffects'
        % Mixed-effects model
        [data_me, trInfo] = buildMEmodel(data, options_stats);
        
        %         options_stats2 = options_stats;
        %     for p = 1:3
        %         inds2keep = setdiff(1:3,p);
        %         options_stats2.fields = options_stats.fields(inds2keep);
        %         pairwiseStats.(['Pair' num2str(p)]).Conditions = options_stats.fields(inds2keep);
        %         [~, all_statVal_pair, all_clusters_pair, ~] = clusteredStats(data4stats(inds2keep), options_stats2);
        %
        %         pairwiseStats.(['Pair' num2str(p)]).StatValsData = all_statVal_pair;
        %         pairwiseStats.(['Pair' num2str(p)]).ClustersData = all_clusters_pair;
        %     end
        
        MEterm = options_stats.MEterm;
        formula = options_stats.formula;
        fields = options_stats.fields;
        main_values = nan(size(data_me,1),size(data_me,2));
        siz = size(main_values);
        main_thresholds = nan(siz);
        
        % One iteration to get info about the model
        % Create appopriate data structure
        tbl = [trInfo, table(squeeze(data_me(1,1,:)), 'VariableNames', {'Data'})];
        % Fit lme
        lme = fitlme(tbl, formula, 'DummyVarCoding', 'effects');
        
        if n_conds > 2
            n_pairs = size(options_stats.pairs,1);
            pair = cell(1,n_pairs);
            if n_conds == 3
                % Specific treatment for backward compatibility
                pair{1} = [3,2];
                %values_pair{1} = nan(siz);
                %thresholds_pair1 = nan(siz);
                pair{2} = [3,1];
                %values_pair{2} = nan(siz);
                %thresholds_pair2 = nan(siz);
                pair{3} = [2,1];
                %values_pair{3} = nan(siz);
                %thresholds_pair3 = nan(siz);
            else
                for p = 1:n_pairs
                    pair{p} = options_stats.pairs(p,:);
                    %values_pair{p} = nan(siz);
                end
            end
            values_pair = nan([n_pairs,siz]);
            thresholds_pairs = nan(siz);
            
            % Prepare for pairwise statistics
            vect_rep = cell(1,n_conds);
            if contains(MEterm,':')
                % Interaction contrast
                MEterms = strsplit(MEterm,':');
                all_coeffs_main1 =  contains(lme.CoefficientNames, MEterms{1}) & ~contains(lme.CoefficientNames, ':');
                all_coeffs_main2 =  contains(lme.CoefficientNames, MEterms{2}) & ~contains(lme.CoefficientNames, ':');
                all_coeffs_inter =  contains(lme.CoefficientNames, MEterms{1}) & contains(lme.CoefficientNames, MEterms{2});
                for c = 1:n_conds
                    vect_rep{c} = zeros(1,lme.NumCoefficients);
                    fields_byTerm = strsplit(fields{c},'_');
                    all_coeffs_field1 = contains(lme.CoefficientNames, fields_byTerm{1});
                    all_coeffs_field2 = contains(lme.CoefficientNames, fields_byTerm{2});
                    
                    % Fill coefficients for first main contrast
                    if any(all_coeffs_field1)
                        vect_rep{c}(all_coeffs_main1 & all_coeffs_field1) = 1;
                    else
                        vect_rep{c}(all_coeffs_main1) = -1;
                    end
                    
                    % Fill coefficients for second main contrast
                    if any(all_coeffs_field2)
                        vect_rep{c}(all_coeffs_main2 & all_coeffs_field2) = 1;
                    else
                        vect_rep{c}(all_coeffs_main2) = -1;
                    end
                    
                    % Fill coefficients for interaction
                    if any(all_coeffs_field1) && any(all_coeffs_field2)
                        vect_rep{c}(all_coeffs_inter & all_coeffs_field1 & all_coeffs_field2) = 1;
                    elseif any(all_coeffs_field1)
                        vect_rep{c}(all_coeffs_inter & all_coeffs_field1) = -1;
                    elseif any(all_coeffs_field2)
                        vect_rep{c}(all_coeffs_inter & all_coeffs_field2) = -1;
                    else
                        vect_rep{c}(all_coeffs_inter) = 1;
                    end
                end
            else
                % Main contrast
                all_coeffs_MEterm = contains(lme.CoefficientNames, MEterm) & ~contains(lme.CoefficientNames, ':');
                for c = 1:n_conds
                    vect_rep{c} = zeros(1,lme.NumCoefficients);
                    all_coeffs_field = contains(lme.CoefficientNames, fields{c});
                    if any(all_coeffs_field)
                        vect_rep{c}(all_coeffs_MEterm & all_coeffs_field) = 1;
                    else
                        vect_rep{c}(all_coeffs_MEterm) = -1;
                    end
                end
            end
        else
            % Necessary to define some variables that will be streamed to the parfor loop
            % even if they are not useful in this case
            n_pairs = 1;
            pair = {};
            vect_rep = {};
        end
        
        if options_stats.original
            % Informations about the model fit:
            fitted_values = permute(nan(size(data_me)),[3,1,2]);
            residuals_values = permute(nan(size(data_me)),[3,1,2]);
            Rsquared_values = nan(siz);
            
            D = designMatrix(lme);
            VIFs = diag(inv(corrcoef(D(:,2:end))));
            CoeffNames = lme.CoefficientNames(2:end);
        end
        
        ppm = ParforProgressbar(numel(main_values), 'showWorkerProgress', true,...
           'title', 'Computing LME for each data point');
        % Loop over each data point
        parfor xy = 1:numel(main_values)
            [x,y] = ind2sub(siz,xy);
            % Create appopriate data structure
            data_xy = squeeze(data_me(x,y,:));
            tbl = [trInfo, table(data_xy, 'VariableNames', {'Data'})];
            
            % Fit lme
            lme = fitlme(tbl, formula, 'DummyVarCoding', 'effects');
            if options_stats.original
                fitted_values(:,xy) = fitted(lme);
                residuals_values(:,xy) = residuals(lme);
                %plotResiduals(lme, 'fitted', 'ResidualType', 'Standardized')
                Rsquared_values(xy) = lme.Rsquared.Adjusted;
            end
            
            if n_conds == 2
                line = find(contains(lme.CoefficientNames, MEterm)...
                    & ~contains(lme.CoefficientNames,':'));
                main_value_temp = table2array(lme.Coefficients(line,4));
                main_threshold_temp = tinv(alpha_crit, table2array(lme.Coefficients(line,5)));
            else
                stats = anova(lme);
                % store values
                line = find(strcmp(stats.Term, MEterm));
                main_value_temp = stats.FStat(line);
                main_threshold_temp = finv(alpha_crit, stats.DF1(line), stats.DF2(line));
                
                %% Pairwise statistics                
                for p=1:n_pairs
                    if p == 1
                        [~,value_pair_temp, DF1, DF2] = coefTest(lme,...
                            vect_rep{pair{p}(1)}-vect_rep{pair{p}(2)});
                    else
                        [~,value_pair_temp, ~, ~] = coefTest(lme,...
                            vect_rep{pair{p}(1)}-vect_rep{pair{p}(2)});
                    end
                    
                    values_pair(p,xy) = value_pair_temp;
                end
                % coefTest returns an F-statistic
                thresholds_pairs(xy) = finv(alpha_crit, DF1, DF2);
            end
            main_values(xy) = main_value_temp;
            main_thresholds(xy) = main_threshold_temp;
            
            %             if n_conds == 3
            %                 values_pair1(xy) = value_pair1_temp;
            %                 values_pair2(xy) = value_pair2_temp;
            %                 values_pair3(xy) = value_pair3_temp;
            %                 thresholds_pairs(xy) = threshold_pairs_temp;
            %             end
            ppm.increment();
        end
        delete(ppm);
        
        if options_stats.original
            saveMEmodelMetaData(options_stats, CoeffNames, VIFs, fitted_values, residuals_values, Rsquared_values);
        end
        %{
                ppm = ParforProgressbar(size(data_me,2), 'showWorkerProgress', true,...
            'title', 'Computing LME for each data point');
        parfor y = 1:size(data_me,2)
            main_values_y = nan(size(data_me,1),1);
            main_thresholds_y = nan(size(data_me,1),1);
            
            if n_conds == 3
                values_pair1_y = nan(size(data_me,1),1);
                %thresholds_pair1_y = nan(size(data_me,1),1);
                values_pair2_y = nan(size(data_me,1),1);
                %thresholds_pair2_y = nan(size(data_me,1),1);
                values_pair3_y = nan(size(data_me,1),1);
                %thresholds_pair3_y = nan(size(data_me,1),1);
                thresholds_pairs_y = nan(size(data_me,1),1);
            end
            
            for x = 1:size(data_me,1)
                % Create appopriate data structure
                data_xy = squeeze(data_me(x,y,:));
                tbl = [trInfo, table(data_xy, 'VariableNames', {'Data'})];
                
                % Fit lme
                lme = fitlme(tbl, formula, 'DummyVarCoding', 'effects');
                if n_conds == 2
                    line = find(contains(lme.CoefficientNames, MEterm)...
                        & ~contains(lme.CoefficientNames,':'));
                    main_values_y(x) = table2array(lme.Coefficients(line,4));
                    main_thresholds_y(x) = tinv(alpha_crit, table2array(lme.Coefficients(line,5)));
                else
                    stats = anova(lme);
                    % store values
                    line = find(strcmp(stats.Term, MEterm));
                    main_values_y(x) = stats.FStat(line);
                    main_thresholds_y(x) = finv(alpha_crit, stats.DF1(line), stats.DF2(line));
                    
                    % Pairwise statistics
                    cond1 = find(contains(lme.CoefficientNames, MEterm)...
                        & contains(lme.CoefficientNames, fields{1})...
                        & ~contains(lme.CoefficientNames,':'));
                    vect1 = zeros(1,lme.NumCoefficients);
                    vect1(1) = 1; % Due to the 'effects' var coding
                    vect1(cond1) = 1;
                    cond2 = find(contains(lme.CoefficientNames, MEterm)...
                        & contains(lme.CoefficientNames, fields{2})...
                        & ~contains(lme.CoefficientNames,':'));
                    vect2 = zeros(1,lme.NumCoefficients);
                    vect2(1) = 1; % Due to the 'effects' var coding
                    vect2(cond2) = 1;
                    %                     cond3 = find(contains(lme.CoefficientNames, MEterm)...
                    %                         & contains(lme.CoefficientNames, options_stats.fields{3})...
                    %                         & ~contains(lme.CoefficientNames,':'));
                    vect3 = zeros(1,lme.NumCoefficients);
                    vect3(1) = 1; % Due to the 'effects' var coding
                    vect3(cond1) = -1;
                    vect3(cond2) = -1;
                    
                    [~,values_pair1_y(x), ~] = coefTest(lme,vect3-vect2);
                    %thresholds_pair1_y(x) = tinv(alpha_crit, df1);
                    [~,values_pair2_y(x), ~] = coefTest(lme,vect3-vect1);
                    %thresholds_pair2_y(x) = tinv(alpha_crit, df2);
                    [~,values_pair3_y(x), ~] = coefTest(lme,vect2-vect1);
                    %thresholds_pair3_y(x) = tinv(alpha_crit, df3);
                    thresholds_pairs_y(x) = tinv(alpha_crit, lme.DFE);
                end
            end
            main_values(:,y) = main_values_y;
            main_thresholds(:,y) = main_thresholds_y;
            
            if n_conds == 3
                values_pair1(:,y) = values_pair1_y;
                %thresholds_pair1(:,y) = thresholds_pair1_y;
                values_pair2(:,y) = values_pair2_y;
                %thresholds_pair2(:,y) = thresholds_pair2_y;
                values_pair3(:,y) = values_pair3_y;
                %thresholds_pair3(:,y) = thresholds_pair3_y;
                thresholds_pairs(:,y) = thresholds_pairs_y;
            end
            ppm.increment();
        end
        delete(ppm);
        %}
        main_thresh = mean(main_thresholds,'all');
        if n_conds == 2
            switch options_stats.style
                case 'timeXfreq'
                    [clusters, cluster_count] = formClustersTimeFreq(main_values, main_thresh, n_conds);
                otherwise
                    error('Unknown style')
            end
            
            [maxstat_val, cluster_vals, clusters] = computeClusterStats(main_values, clusters, cluster_count,...
                options_stats.removeSmallestClusters);
            %             [maxstat_val, cluster_vals, clusters] = computeClusterStatsTimeFreq(main_values, main_thresh, n_conds,...
            %                 options_stats.removeSmallestClusters);
        else
            switch options_stats.style
                case 'timeXfreq'
                    [clusters_main, cluster_count] = formClustersTimeFreq(main_values, main_thresh, n_conds);
                otherwise
                    error('Unknown style')
            end
            [maxstat_val.Main, cluster_vals.Main, clusters.Main] = computeClusterStats(main_values, clusters_main, cluster_count,...
                options_stats.removeSmallestClusters);
            %             [maxstat_val.Main, cluster_vals.Main, clusters.Main] = computeClusterStatsTimeFreq(main_values, main_thresh, n_conds,...
            %                 options_stats.removeSmallestClusters);
            
            thresh_pairs = mean(thresholds_pairs,'all');
            for p = 1:n_pairs
                switch options_stats.style
                    case 'timeXfreq'
                        [clusters_pair, cluster_count] = formClustersTimeFreq(squeeze(values_pair(p,:,:)), thresh_pairs, n_conds);
                    otherwise
                        error('Unknown style')
                end
                [maxstat_val.(sprintf('Pair%d',p)), cluster_vals.(sprintf('Pair%d',p)), clusters.(sprintf('Pair%d',p))] =...
                    computeClusterStats(squeeze(values_pair(p,:,:)), clusters_pair, cluster_count, options_stats.removeSmallestClusters);
                %             [maxstat_val.Pair1, cluster_vals.Pair1, clusters.Pair1] = computeClusterStatsTimeFreq(values_pair1, thresh_pairs, n_conds,...
                %                 options_stats.removeSmallestClusters);
            end
        end
end

% % Select stats above threshold
% if n_conds == 2
%     aboveTh = abs(values)>thresh;
% else
%     aboveTh = values>thresh;
% end

% %% Form clusters based on temporal/frequency adjacency
% clusters = zeros(size(values));
% cluster_count = 0;
% eq_clusters = [];
% for i = 1:size(values,1)
%     for j = 1:size(values,2)
%         if aboveTh(i,j)
%             % Check if one of the neighbours already belongs to a cluster
%             if j>1 && clusters(i,j-1)~=0 && i>1 && clusters(i-1,j)~=0
%                 % Conflicting clusters?
%                 if clusters(i,j-1) == clusters(i-1,j)
%                     % No
%                     clusters(i,j) = clusters(i-1,j);
%                 else
%                     % Yes, put them in the equivalent clusters list
%                     eq_clusters = [eq_clusters;[clusters(i-1,j),clusters(i,j-1)]];
%                     clusters(i,j) = clusters(i-1,j);
%                 end
%             elseif j>1 && clusters(i,j-1)~=0
%                 clusters(i,j) = clusters(i,j-1);
%             elseif i>1 && clusters(i-1,j)~=0
%                 clusters(i,j) = clusters(i-1,j);
%             else
%                 % Create a new cluster
%                 cluster_count = cluster_count+1;
%                 clusters(i,j) = cluster_count;
%             end
%         end
%     end
% end
%
% % Merge equivalent clusters for the neighbours of this channel
% if ~isempty(eq_clusters)
%     while ~isempty(intersect(unique(eq_clusters(:,1)), unique(eq_clusters(:,2))))
%         % Sort the equivalent clusters such as the lowest id is always in the
%         % first column
%         eq_clusters = sort(eq_clusters,2);
%
%         clusters2merge = unique(eq_clusters(:));
%
%         for cl = clusters2merge(2:end)' % Force the vector to be a line
%             if sum(eq_clusters(:) == cl) > 1
%                 lowestEq = find(eq_clusters(:,2) == cl,1);
%                 if ~isempty(lowestEq)
%                     duplicates = find(eq_clusters(:) == cl);
%                     IndtoAvoid = size(eq_clusters,1) + lowestEq;
%
%                     IndsToReplace = setdiff(duplicates,IndtoAvoid);
%                     eq_clusters(IndsToReplace) = eq_clusters(lowestEq,1);
%                 end
%             end
%         end
%
%         % Remove lines where both clusters are equal
%         eq_clusters(eq_clusters(:,1)==eq_clusters(:,2),:) = [];
%     end
%
%     for eq = 1:size(eq_clusters,1)
%         clusters(clusters == eq_clusters(eq,2)) = eq_clusters(eq,1);
%     end
% end

% if cluster_count == 0
%     maxstat_val = NaN;
%     cluster_vals = [];
% else
%     cluster_count = length(unique(clusters))-1; % There still is 0 values in the matrix
%
%     clust_inds = sort(unique(clusters));
%     clust_inds = clust_inds(2:end);
%
%     if options_stats.removeSmallestClusters
%         removedClust_count = 0;
%         %% Remove clusters that are too small (less than 0.1% of the samples)
%         for cl = 1:cluster_count
%             if numel(find(clusters == clust_inds(cl))) <= critSize
%                 % Remove the cluster
%                 clusters(clusters == clust_inds(cl)) = 0;
%                 removedClust_count = removedClust_count + 1;
%             else
%                 % Renumber the cluster from 1 to cluster_count
%                 clusters(clusters == clust_inds(cl)) = cl - removedClust_count;
%             end
%         end
%         cluster_count = cluster_count - removedClust_count;
%     end
%
%     %% Give the final clusters an index between 1 and cluster_count
%     for cl = 1:cluster_count
%         clusters(clusters == clust_inds(cl)) = cl;
%     end
%
%     %% Calculate cluster levels stats
%     cluster_vals = zeros(cluster_count,1);
%     for cl = 1:cluster_count
%         cluster_vals(cl) = sum(values(clusters == cl));
%     end
%
%     %% Take the largest
%     [~,c_max] = max(abs(cluster_vals));
%
%     % if length(c_max)>1 % multiple maximums
%     %     %Take the largest cluster
%     %     clusters_size = sum(clusters == c_max);
%     %     [~,c_max] = max(clusters_size);
%     %     if length(c_max)>1
%     %         warning('Still multiple maxima')
%     %     end
%     % end
%
%     maxstat_val = cluster_vals(c_max);
%     maxstat_val = abs(maxstat_val(1)); %In case there are multiple maxima
% end

    function [data_me, trInfo] = buildMEmodel(data, options_stats)
        data_me = [];
        for f = 1:numel(options_stats.fields)
            data_me = cat(3,data_me,data{f});
            trInfo_cond = options_stats.(sprintf('trialsInfo_%s', options_stats.fields{f}));
            if f==1
                if strcmp(options_stats.MEterm,'Baseline')
                    if contains(options_stats.fields{f},'base')
                        trInfo = [table(true(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    else
                        trInfo = [table(false(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    end
                else
                    trInfo = trInfo_cond;
                end
            else
                if strcmp(options_stats.MEterm,'Baseline')
                    if contains(options_stats.fields{f},'base')
                        trInfo_temp = [table(true(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    else
                        trInfo_temp = [table(false(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    end
                    trInfo = [trInfo; trInfo_temp];
                else
                    trInfo = [trInfo; trInfo_cond];
                end
            end
        end
    end

    function saveMEmodelMetaData(options, CoeffNames, VIFs, Fitted, Residuals, Rsquared)
        LMEsummary = struct;
        LMEsummary.Formula = options.formula;
        LMEsummary.CoeffNames = CoeffNames;
        LMEsummary.VIFs = VIFs;
        LMEsummary.fitted = Fitted;
        LMEsummary.residuals = Residuals;
        LMEsummary.Rsquare = Rsquared;
        
        if strcmp(options.MEterm, 'Baseline')
            save(fullfile(options.saveFolder, sprintf('lmeSummary_%s_%s_%svsBaseline',...
                options.ROI, options.modelData, options.fields{1})), 'LMEsummary');
        else
            if contains(options.MEterm, ':')
                MEterm_list = strsplit(options.MEterm,':');
                save(fullfile(options.saveFolder, sprintf('lmeSummary_%s_%s_%sX%sinteraction',...
                    options.ROI, options.modelData, MEterm_list{1}, MEterm_list{2})), 'LMEsummary');
            else
                save(fullfile(options.saveFolder, sprintf('lmeSummary_%s_%s_%scontrast',...
                    options.ROI, options.modelData, options.MEterm)), 'LMEsummary');
            end
        end
    end


    function [clusters, cluster_count] = formClustersUniDim(values, thresh, n_conds)
        % Select stats above threshold
        if n_conds == 2
            aboveTh = abs(values)>thresh;
        else
            aboveTh = values>thresh;
        end
        
        %% Form clusters based on temporal adjacency
        clusters = zeros(size(aboveTh));
        cluster_count = 0;
        for i = 1:length(aboveTh)
            if aboveTh(i)
                if i>1 && clusters(i-1)~=0
                    % Connect to existing cluster
                    clusters(i) = clusters(i-1);
                else
                    % Create a new cluster
                    cluster_count = cluster_count+1;
                    clusters(i) = cluster_count;
                end
            end
        end
        
    end

function [clusters, cluster_count] = formClustersTimeFreq(values, thresh, n_conds)
        % Select stats above threshold
        if n_conds == 2
            aboveTh = abs(values)>thresh;
        else
            aboveTh = values>thresh;
        end
        
        %% Form clusters based on temporal/frequency adjacency
        clusters = zeros(size(aboveTh));
        cluster_count = 0;
        eq_clusters = [];
        for i = 1:size(aboveTh,1)
            for j = 1:size(aboveTh,2)
                if aboveTh(i,j)
                    % Check if one of the neighbours already belongs to a cluster
                    if j>1 && clusters(i,j-1)~=0 && i>1 && clusters(i-1,j)~=0
                        % Conflicting clusters?
                        if clusters(i,j-1) == clusters(i-1,j)
                            % No
                            clusters(i,j) = clusters(i-1,j);
                        else
                            % Yes, put them in the equivalent clusters list
                            eq_clusters = [eq_clusters;[clusters(i-1,j),clusters(i,j-1)]];
                            clusters(i,j) = clusters(i-1,j);
                        end
                    elseif j>1 && clusters(i,j-1)~=0
                        clusters(i,j) = clusters(i,j-1);
                    elseif i>1 && clusters(i-1,j)~=0
                        clusters(i,j) = clusters(i-1,j);
                    else
                        % Create a new cluster
                        cluster_count = cluster_count+1;
                        clusters(i,j) = cluster_count;
                    end
                end
            end
        end
        % Merge equivalent clusters for the neighbours of this channel
        if ~isempty(eq_clusters)
            while ~isempty(intersect(unique(eq_clusters(:,1)), unique(eq_clusters(:,2))))
                % Sort the equivalent clusters such as the lowest id is always in the
                % first column
                eq_clusters = sort(eq_clusters,2);
                
                clusters2merge = unique(eq_clusters(:));
                
                for cl = clusters2merge(2:end)' % Force the vector to be a line
                    if sum(eq_clusters(:) == cl) > 1
                        lowestEq = find(eq_clusters(:,2) == cl,1);
                        if ~isempty(lowestEq)
                            duplicates = find(eq_clusters(:) == cl);
                            IndtoAvoid = size(eq_clusters,1) + lowestEq;
                            
                            IndsToReplace = setdiff(duplicates,IndtoAvoid);
                            eq_clusters(IndsToReplace) = eq_clusters(lowestEq,1);
                        end
                    end
                end
                
                % Remove lines where both clusters are equal
                eq_clusters(eq_clusters(:,1)==eq_clusters(:,2),:) = [];
            end
            
            for eq = 1:size(eq_clusters,1)
                clusters(clusters == eq_clusters(eq,2)) = eq_clusters(eq,1);
            end
        end
    end

    function [clusters, cluster_count] = formClustersChanFreq(values, thresh, n_conds, Channels, ElecFile, MaxDeg)
        % Select stats above threshold
        if n_conds == 2
            aboveTh = abs(values)>thresh;
        else
            aboveTh = values>thresh;
        end
        
        n_chans = size(values,1);
        n_freqs = size(values,2);
        %% Form clusters based on frequency adjacency
        clusters = zeros(n_chans,n_freqs);
        cluster_count = 0;
        for ch = 1:n_chans
            for f = 1:n_freqs
                if aboveTh(ch,f)
                    % Only spectral clustering for now
                    if f>1 && clusters(ch,f-1)~=0
                        % Check if previous time is already part of a cluster
                        clusters(ch,f) = clusters(ch,f-1);
                    else
                        % Create a new cluster
                        cluster_count = cluster_count+1;
                        clusters(ch,f) = cluster_count;
                    end
                end
            end
        end
        
        %% Group clusters by neighbouring channels
        % Read File where electrode coordinates are encoded
        % Here the spherical coordinates given by the manufacturer are used
        % (Not the precise individual channel location files)
        cap_coords = readeetraklocs(ElecFile);
        
        eq_clusters = [];
        for ch = 1:n_chans
            if iscell(Channels)
                chan = Channels{ch};
            else
                chan = Channels(ch);
            end
            neighbours = defineNeighbours(cap_coords, chan, MaxDeg);
            % Select neighbours that are in the subset
            neighbours = intersect(neighbours, Channels);
            
            if ~isempty(neighbours)
                for n = neighbours(:)' % Force the vector to be a line
                    
                    if iscell(Channels)
                        n_ind = find(strcmp(Channels,n));
                    else
                        n_ind = find(Channels == n);
                    end
                    
                    if n_ind > ch % otherwise the pair has already been inspected
                        % Keep track of the current clusters under inspection
                        curr_cluster_ch = 0;
                        curr_cluster_n = 0;
                        for f = 1:n_freqs
                            if clusters(ch,f) ~= 0 && clusters(n_ind,f) ~= 0 && clusters(n_ind,f)~=clusters(ch,f)
                                if curr_cluster_ch ~= clusters(ch,f) || curr_cluster_n ~= clusters(n_ind,f)
                                    % This pair was not yet encountered
                                    eq_clusters = [eq_clusters;[clusters(ch,f),clusters(n_ind,f)]];
                                end
                            end
                            curr_cluster_ch = clusters(ch,f);
                            curr_cluster_n = clusters(n_ind,f);
                        end
                    end
                end
            end
        end
        % Merge equivalent clusters
        if ~isempty(eq_clusters)
            while ~isempty(intersect(unique(eq_clusters(:,1)), unique(eq_clusters(:,2))))
                % Sort the equivalent clusters such as the lowest id is always in the
                % first column
                eq_clusters = sort(eq_clusters,2);
                
                clusters2merge = unique(eq_clusters(:));
                
                for cl = clusters2merge(2:end)' % Force the vector to be a line
                    if sum(eq_clusters(:) == cl) > 1
                        lowestEq = find(eq_clusters(:,2) == cl,1);
                        if ~isempty(lowestEq)
                            duplicates = find(eq_clusters(:) == cl);
                            IndtoAvoid = size(eq_clusters,1) + lowestEq;
                            
                            IndsToReplace = setdiff(duplicates,IndtoAvoid);
                            eq_clusters(IndsToReplace) = eq_clusters(lowestEq,1);
                        end
                    end
                end
                
                % Remove lines where both clusters are equal
                eq_clusters(eq_clusters(:,1)==eq_clusters(:,2),:) = [];
            end
            
            for eq = 1:size(eq_clusters,1)
                clusters(clusters == eq_clusters(eq,2)) = eq_clusters(eq,1);
            end
        end
    end

    function [maxstat_val, cluster_vals, clusters] = computeClusterStats(values, clusters, cluster_count, rmSmallClusters)
        if cluster_count == 0
            maxstat_val = NaN;
            cluster_vals = [];
        else
            clust_inds = sort(unique(clusters));
            if any(clusters==0, 'all')
                cluster_count = length(unique(clusters))-1; % There still is 0 values in the data
                clust_inds = clust_inds(2:end);
            end

            if rmSmallClusters
                %% Remove clusters that are too small (less than 0.1% of the samples)
                removedClust_count = 0;
                critSize = 0.001*numel(clusters);
                for cl = 1:cluster_count
                    if numel(find(clusters == clust_inds(cl))) <= critSize
                        % Remove the cluster
                        clusters(clusters == clust_inds(cl)) = 0;
                        removedClust_count = removedClust_count + 1;
                    else
                        % Renumber the cluster from 1 to cluster_count
                        clusters(clusters == clust_inds(cl)) = cl - removedClust_count;
                    end
                end
                cluster_count = cluster_count - removedClust_count;
            end
            
            %% Give the final clusters an index between 1 and cluster_count
            for cl = 1:cluster_count
                clusters(clusters == clust_inds(cl)) = cl;
            end
            
            %% Calculate cluster levels stats
            cluster_vals = zeros(cluster_count,1);
            for cl = 1:cluster_count
                cluster_vals(cl) = sum(values(clusters == cl));
            end
            
            %% Take the largest
            [~,c_max] = max(abs(cluster_vals));
            
            % if length(c_max)>1 % multiple maximums
            %     %Take the largest cluster
            %     clusters_size = sum(clusters == c_max);
            %     [~,c_max] = max(clusters_size);
            %     if length(c_max)>1
            %         warning('Still multiple maxima')
            %     end
            % end
            
            maxstat_val = cluster_vals(c_max);
            maxstat_val = abs(maxstat_val(1)); %In case there are multiple maxima
        end
    end

    function neighbours = defineNeighbours(cap_coords, channel, max_deg_dev)
        %% Locate the channel
        if ischar(channel)
            ch_ref = find(strcmp({cap_coords(:).labels},channel));
        else
            ch_ref = channel;
        end
        vect_ref = [cap_coords(ch_ref).X, cap_coords(ch_ref).Y, cap_coords(ch_ref).Z];
        vect_ref = vect_ref./norm(vect_ref); % Normalise vector
        
        %% Apply the neighbouring rule
        dot_crit = cos(max_deg_dev*pi()/180);
        neighbours = [];
        for ch = 1:numel(cap_coords)
            if ch ~= ch_ref
                vect = [cap_coords(ch).X, cap_coords(ch).Y, cap_coords(ch).Z];
                vect = vect./norm(vect); % Normalise vector
                if dot(vect, vect_ref) > dot_crit
                    neighbours = [neighbours,ch];
                end
            end
        end
        
        if ischar(channel)
            %Give the result as a cell of chars
            neighbours = {cap_coords(neighbours).labels};
        end
    end
end
