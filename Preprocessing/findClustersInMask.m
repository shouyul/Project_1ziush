function [all_inds, clusters] = findClustersInMask(mask)
% Helper function for fincHighCorrChannels and get_elecdists
% mask should be a square matrix (symmetrical)

[inds1,inds2] = find(mask);
all_inds = unique([inds1,inds2]);

clusters = {};
pairedWith = cell(length(all_inds),1);
for i = 1:length(all_inds)
    pairedWith{i} = union(inds2(inds1 == all_inds(i)), inds1(inds2 == all_inds(i)));
    % make sure this is a row vector:
    if size(pairedWith{i},1) > 1
        pairedWith{i} = pairedWith{i}';
    end
    
    % Form clusters
    numExistingCl = 0;
    for cl = 1:numel(clusters)
        if any(clusters{cl} == all_inds(i))
            numExistingCl = numExistingCl + 1;
            clusters{cl} = union(clusters{cl}, [all_inds(i),pairedWith{i}]);            
            continue
        end
        
        if any(clusters{cl} == pairedWith{i}', 'all')
            numExistingCl = numExistingCl + 1;
            clusters{cl} = union(clusters{cl}, [all_inds(i),pairedWith{i}]);
        end
    end
    
    if numExistingCl == 0
        if numel(clusters) == 0
            clusters{1} = [all_inds(i),pairedWith{i}];
        else
            clusters{end+1} = [all_inds(i),pairedWith{i}];
        end
    elseif numExistingCl > 1
        % Some clusters are overlapping: merge them
        eq_clusters = [];
        merge = [];
        for cl = 1:numel(clusters)
            if any(clusters{cl} == all_inds(i))
                eq_clusters = union(eq_clusters, cl);
                if isempty(merge)
                    merge = clusters{cl};
                else
                    merge = union(merge, clusters{cl});
                end
            end
        end       
        
        clusters{eq_clusters(1)} = merge;        
        clusters = clusters(setdiff(1:numel(clusters),eq_clusters(2:end)));
    end
end


end