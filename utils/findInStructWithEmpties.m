function res = findInStructWithEmpties(S, fieldName, value)
try
    idx = false(size(S));
    SS = {S.(fieldName)};
    inds = ~cellfun('isempty', SS);
    idx(inds) = [SS{inds}]==value;
    res = find(idx);
catch
    disp('Problem with findInStructWithEmpties function')
end
end