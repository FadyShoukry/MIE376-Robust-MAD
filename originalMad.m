% This script implements a robust version of the Mean Absolute Deviation
% (MAD) financial optimization model.

clear
clc

%% Fetching Required Data

connect = yahoo; %the source
assets = {'GOOGL', 'AAPL', 'ALTR', 'FB', 'YHOO', 'GS', 'TXN', 'IBM', 'SSNLF', 'MSIQX'}; 

% Get the data from Yahoo
for i = 1:length(assets);
    tmp = fetch(connect, assets{i}, 'adj close', 'Aug 3 2013', 'Feb 3 2014', 'd');
    data(:, i) = tmp(:, 2);
end

% Re-order the data from oldest to newest
data = flipud(data);

% Get the size of the data
% N is the number of assets
% T is the number of periods
[T, n] = size(data);
T = T - 1;

%% Calculating the required model parameters

% Calculate period interest rates
r_it = (data(2:end,:)./data(1:end-1,:))-1;

% Calculate mu using geometric mean
mu = geomean(1 + r_it) - 1;

% Required interest rate
R = 0.05;

%% MAD by linprog
disp('MAD')
c = [zeros(n,1); ones(T,1); ones(T,1)];% [x_i; y_t; z_t]
Aeq = [r_it-repmat(mu,T,1) -eye(T) eye(T); 
    ones(1,n) zeros(1,2*T);];% the constraint coefficient of MAD
beq=[zeros(T,1); 1;];
lb =[ones(n, 1)*-inf; zeros(T+T,1)];         % the lower bound of the variables

A=-[mu zeros(1,2*T);];
b=-R;

%% Optimize the model

[x, fval] = linprog(c, A, b, Aeq, beq, lb, []);

