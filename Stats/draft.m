
% MUST BE MERGED WITH EXISTING STAT SCRIPTS.
% quick and dirty behavioral stats for PIONEER-EEG presentation on March 20 2024

P1001_off = [0 20; 0 20];
P1001_on = [20 0; 1 19];
P1002_off = [0 20; 0 20];
P1002_on = [20 0; 10 10];
P1004_off = [15 5; 14 6];
P1004_on = [19 1; 7 13];
P1009_off = [10 10; 2 18];
P1009_on = [0 20; 0 20];


n = 40; % Total number of trials per condition

disp('P1001 – Binomial');
[~,p_value,~] = binotest(P1001_off(1,1) + P1001_off(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);
[~,p_value,~] = binotest(P1001_on(1,1) + P1001_on(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);

disp('P1002 – Binomial');
[~,p_value,~] = binotest(P1002_off(1,1) + P1002_off(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);
[~,p_value,~] = binotest(P1002_on(1,1) + P1002_on(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);

disp('P1004 – Binomial');
[~,p_value,~] = binotest(P1004_off(1,1) + P1004_off(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);
[~,p_value,~] = binotest(P1004_on(1,1) + P1004_on(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);

disp('P1009 – Binomial');
[~,p_value,~] = binotest(P1009_off(1,1) + P1009_off(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);
[~,p_value,~] = binotest(P1009_on(1,1) + P1009_on(2,2), n, 0.5); % chance level is 0.5 here
disp(p_value);

disp('P1001 – Fisher');
[~,p,~] = fishertest(P1001_off);
disp(p);
[~,p,~] = fishertest(P1001_on);
disp(p);

disp('P1002 – Fisher');
[~,p,~] = fishertest(P1002_off);
disp(p);
[~,p,~] = fishertest(P1002_on);
disp(p);

disp('P1004 – Fisher');
[~,p,~] = fishertest(P1004_off);
disp(p);
[~,p,~] = fishertest(P1004_on);
disp(p);

disp('P1009 – Fisher');
[~,p,~] = fishertest(P1009_off);
disp(p);
[~,p,~] = fishertest(P1009_on);
disp(p);







