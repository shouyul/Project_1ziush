%% Bemobil-APP pipeline
app_KC = struct('ID', [], 'Brain', [], 'BrainWithNoise', [], 'NoiseWithBrain',[], 'Doubts', []);

app_KC(1) = struct('ID', 'P1001',...
    'Brain', [1:4,8,12,20,22,43,47],...
    'BrainWithNoise', [5,7,9,14:15,21,34,41:42,49,52],...
    'NoiseWithBrain', [16,24,26,28,37:38,46,54:56],...
    'Doubts', [1,21,49]);

% 7&9 are complementary
% Components that seem influenced by goggles wearing:
% 3,7,14,21
% PCA reduction to 60 ICs

app_KC(end+1) = struct('ID', 'P1004',...
    'Brain', [1,3,12,17],...
    'BrainWithNoise', [30,35,40:41],...
    'NoiseWithBrain', [5,21:22,32,44,50,56],...
    'Doubts', [1,3,30]);

% 4&7 are complementary
% 12&17: 'triple wave'
% Components that seem influenced by goggles wearing:
% 1,3,4,7,12,17,35,40
% PCA reduction to 60 ICs

%% Bemobil-ASR pipeline
asr_KC = struct('ID', [], 'Brain', [], 'BrainWithNoise', [], 'NoiseWithBrain',[], 'Doubts', []);



