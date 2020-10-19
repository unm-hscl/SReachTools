%% Run the performance tests.
profileResults = runperf('tests/profile/ProfileKernelEmbeddings.m');

%% Display the results.

fullTable = vertcat(profileResults.Samples);

results = varfun(@mean, fullTable, ...
    'InputVariables', 'MeasuredTime', ...
    'GroupingVariables', 'Name');

disp(results);

%% Run the performance tests.
profileResultsRFF = runperf('tests/profile/ProfileKernelEmbeddingsRFF.m');

%% Display the results.

fullTable = vertcat(profileResultsRFF.Samples);

resultsRFF = varfun(@mean, fullTable, ...
    'InputVariables', 'MeasuredTime', ...
    'GroupingVariables', 'Name');

disp(resultsRFF);

%% Plot the results.

% The following code assumes there are 25 values of M and 50 values of D.

M = 100:100:2500;
KernelEmbeddingsTime = zeros(size(M, 2), 1);

for p = 1:size(M, 2)
    KernelEmbeddingsTime(p) = table2array(results(p, 3));
end

D = 100:100:1000;
KernelEmbeddingsRFFTime = zeros(size(M, 2), size(D, 2));

for p = 1:size(D, 2)
    KernelEmbeddingsRFFTime(:, p) = ...
        table2array(resultsRFF(p:size(D, 2):size(resultsRFF, 1), 3));
end

%%

figure
hold on
plot(M, KernelEmbeddingsTime);
hold off

%%

figure 
hold on
for p = 1:size(D, 2)
    plot(M, KernelEmbeddingsRFFTime(:, p));
end
hold off