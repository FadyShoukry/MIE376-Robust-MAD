function [x, fval] = RCMAD(r_it, confidence, R)
%RCMAD implements a robust version of the Mean Absolute Deviation
% (MAD) financial optimization model.

    % Get the size of the data
    % N is the number of assets
    % T is the number of periods
    [T, N] = size(r_it);

    %% Calculating the required model parameters

    % Confidence level
    alpha = 1-confidence;

    % Calculate mu using geometric mean
    mu = geomean(1 + r_it) - 1;

    % Calculate the sample standard deviation for each asset
    sigma = std(r_it);

    % Calculate the standard error in mean for a (1-alpha) confidence interval
    s_e = tinv(alpha/2, T-1)*sigma./sqrt(T);

    %% Formulate the model

    % Cost function
    f = [ones(T, 1); zeros(N, 1)];

    % Inequality Constraints
    % The following are repitition of my and s_e needed for
    % the construction of the A Matrix of inequality constraints
    mu_rep = repmat(mu, T, 1);
    s_e_rep = repmat(s_e, T, 1);

    A = [-eye(T), -r_it + mu_rep + s_e_rep;
         -eye(T), -r_it + mu_rep - s_e_rep;
         -eye(T),  r_it - mu_rep - s_e_rep;
         -eye(T),  r_it - mu_rep + s_e_rep;
         zeros(1, T), -mu - s_e;
         zeros(1, T), -mu + s_e];

    b = [zeros(T*4, 1); -R; -R];
    % Euqality constraint
    % sum(xi) = 1
    Aeq = [zeros(1, T), ones(1, N)];
    beq = 1;

    % Lower and Upper bounds
    lb = [zeros(T, 1); -inf*ones(N, 1)];

    %% Optimize the model
    
    %set the options
    options.Display = 'off';

    [x, fval] = linprog(f, A, b, Aeq, beq, lb, [], 0, options);
    x = x(T+1 : T+N);

end

