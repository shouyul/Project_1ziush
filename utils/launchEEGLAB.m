%% Needed to start eeglab on linux
% Make sure graph2d\axis.p function is first in the matlab path
path= fullfile(matlabroot, 'toolbox', 'matlab', 'graph2d');
addpath(genpath(path), '-begin')
which axis % Should return something like 'C:\Program Files\MATLAB\R2019a\toolbox\matlab\graph2d\axis.p

%% Options to start eeglab:
eeglab;
pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1,...
    'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0,...
    'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1,...
    'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 0);