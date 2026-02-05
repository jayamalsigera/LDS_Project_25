%% Test Connectivity Percentage Distribution
% Run multiple iterations to analyze the distribution of connectivity
% percentages across randomly generated spatial networks

numIterations = 10000;
numNodes = 100;
numSensors = 20;
maxLength = 5000;

percentages = zeros(numIterations, 1);

h = waitbar(0, 'Running iterations...');
for i = 1:numIterations
  G = createSpatialNetwork(numNodes, numSensors, maxLength);
  percentages(i) = getConnectivityPercentage(G);
  waitbar(i/numIterations, h, sprintf('Iteration %d/%d', i, numIterations));
end
close(h);

stats = table( ...
  mean(percentages), ...
  median(percentages), ...
  std(percentages), ...
  min(percentages), ...
  max(percentages), ...
  'VariableNames', {'Mean', 'Median', 'StdDev', 'Min', 'Max'}, ...
  'RowNames', {'Connectivity %'});
disp(stats);

% Plot histogram
figure;
histogram(percentages, 30);
xlabel('Connectivity Percentage (%)');
ylabel('Frequency');
title(sprintf('Distribution of Connectivity Percentages (N=%d, %d iterations)', numNodes, numIterations));
grid on;
