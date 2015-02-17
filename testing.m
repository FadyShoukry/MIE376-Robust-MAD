% This script tests the robust and non-robust formulations of the MAD model

clear
clc

%% Fetching Required Data

connect = yahoo; %the source
assets = {'GOOGL', 'AAPL', 'ALTR', 'FB', 'YHOO', 'GS', 'TXN', 'IBM', 'SSNLF', 'MSIQX'}; 

% Get the data from Yahoo
% The modelling data
for i = 1:length(assets);
    tmp = fetch(connect, assets{i}, 'adj close', 'Aug 3 2013', 'Feb 3 2014', 'd');
    model_data(:, i) = tmp(:, 2);
end
% Testing data
for i = 1:length(assets);
    tmp = fetch(connect, assets{i}, 'adj close', 'Feb 4 2014', 'April 4 2014', 'd');
    testing_data(:, i) = tmp(:, 2);
end

% Re-order the data from oldest to newest
model_data = flipud(model_data);
testing_data = flipud(testing_data);

% Calculate the per-period return rates for both data sets
model_r_it = (model_data(2:end,:)./model_data(1:end-1,:))-1;
testing_r_it = (testing_data(2:end,:)./testing_data(1:end-1,:))-1;

% Get the size of the data
% N is the number of assets
% T is the number of periods
[T, N] = size(model_r_it);

% Compute the range of R values for comparison
R = linspace(min(model_r_it(:)), max(model_r_it(:)), 100);


%% Get the optimum solution for both models for a range of R values
% For a confidence level of 50%
% initialize the result vectors for speed
x_RMAD = zeros(N, length(R));
x_MAD = zeros(N, length(R));
fval_RMAD = zeros(1, length(R));
fval_MAD = zeros(1, length(R));

% For a confidence level of 95%
for i = 1:length(R)
    [x_RMAD(:,i), fval_RMAD(i)] = RCMAD(model_r_it, 0.95, R(i));
    [x_MAD(:,i), fval_MAD(i)] = MAD(model_r_it, R(i));
end

%% Plot the 2 sets of data
figure
plot(R, fval_MAD, 'r', ...
     R, fval_RMAD, 'b');
% Make the graph pretty
title('Optimal Objective Value for Different Returns (R)', ...
       'interpreter', 'latex');
ylabel('Optimal Objective Value', ...
       'interpreter', 'latex');
xlabel('Rate of Return (R)', ...
       'interpreter', 'latex');
legend('Non-Robust', 'Robust with 95% Confidence', ...
        'Location', 'NorthWest');

%% Find the optimal portfolios for different values of alpha for a given R
R = mean2(model_r_it); % Use the mean value of all returns as your R
confidence = [0.25, 0.5, 0.75, 0.95];

for a = 1:length(confidence)
    [x_RMAD_test(:,a), fval_RMAD_test(a)] = RCMAD(model_r_it, confidence(a), R);
end

[x_MAD_test, fval_MAD_test] = MAD(model_r_it, R); % find the optimal portfolio with MAD

%% Compare performance over the last 2 months of data
% Compute the return of RMAD portfolios for all values of confidence
% over all the time periods in the last 2 months
RMAD_returns = (testing_r_it*x_RMAD_test)';
MAD_returns = (testing_r_it*x_MAD_test)';

% Compute the difference between the RMAD and MAD returns
RMAD_performance = RMAD_returns - repmat(MAD_returns, 4, 1);

%% Plot the difference for the 4 values of confidence
x = 1:42;
for i = 1:4
    figure
    plot(x, RMAD_performance(i, :), '*r', 'MarkerSize', 10);
    title(['Increase in Performance of RC MAD, $1- \alpha$ = ' num2str(confidence(i))], ...
           'interpreter', 'latex');
    ylabel('$r_{_{RC}} - r_{_{MAD}}$', ...
           'interpreter', 'latex', 'FontSize', 15);
    xlabel('Time period (T)', ...
           'interpreter', 'latex');
    axis tight;
end

%% Plot probability plots to check for normality
for i = 1:4
    figure
    normplot(RMAD_performance(i,:));
    title(['Normality plot for $1- \alpha$ = ' num2str(confidence(i))], ...
            'interpreter', 'latex');
    xlabel('$r_{_{RC}} - r_{_{MAD}}$', ...
           'interpreter', 'latex', 'FontSize', 15);
    ylabel('Probability', 'interpreter', 'latex');
end
%% Compute the mean performance increase over the last 2 months for each alpha
mean_performance = mean(RMAD_performance');


