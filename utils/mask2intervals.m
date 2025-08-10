function intervals = mask2intervals(mask)
% Converts a continuous unidirectional mask (boolean vector) to 
% the corresponding intervals mapping the start and stop of continous
% regions marked by the mask (ones in the the mask).

bounds = diff(mask);
starts = find(bounds == 1)+1;
stops = find(bounds == -1);
if mask(1) == 1
    starts(end+1) = 1;
    starts = sort(starts);
end

if mask(end) == 1
    stops(end+1) = length(mask);
end
intervals = [starts; stops]';
end