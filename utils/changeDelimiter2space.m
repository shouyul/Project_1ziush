function [newcell] = changeDelimiter2space(cell, delimiter)
% Convert the content of a string cell to a new string cell where the specified delimiter has been replaced by space
% Useful in labels for plots to avoid the wrong interpretation of _
newcell = cell(1,numel(cell));
for c=1:numel(cell)
    if ischar(cell{c})
        strArray = split(cell{c},delimiter);
        newstr = '';
        for i = 1:length(strArray)
            newstr = [newstr strArray{i} ' '];
        end
        newstr = newstr(1:end-1); % remove the last space
        newcell{c} = newstr;
    end
end
end

