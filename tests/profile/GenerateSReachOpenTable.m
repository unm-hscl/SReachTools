%% Run the performance tests.
profileResults = runperf('tests/profile/ProfileSReachOpen.m');

%% Display the results. 

fullTable = vertcat(profileResults.Samples);

results = varfun(@mean, fullTable, ...
    'InputVariables', 'MeasuredTime', ...
    'GroupingVariables', 'Name');

disp(results);