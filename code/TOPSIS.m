% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% TOPSIS and Parallel Coordinate Plot
% 30/07/25

% Script to apply TOPSIS and create a Parallel Coordinate Plot for four scenarios
% Data: All_data (80x9) with columns [Nt, As, TR, Vc, N_PV, N_ro, Qp, LCOW, GWP]
% Scenarios: 3 with 80% weight on one objective (Qp, LCOW, GWP) and 10% on others,
%            1 with equal weights (0.333 each)

clc; clear;
% Load your Pareto front solutions (adjust filename if needed)
load('All_data_MOO.mat');

% Define variable names for clarity
var_names = {
    'VAWTs No.', ...      % 1. Number of VAWTs
    'VAWTs Size', ...     % 2. VAWT size scale
    'TR', ...             % 3. Transmission ratio
    'Air comp. Size', ... % 4. Air compressor size scale
    'PV No.', ...         % 5. Number of pressure vessel
    'ROs No.', ...        % 6. Number of ROs
    'Qp', ...             % 7. Annual water production
    'LCOW', ...           % 8. Levelised cost of water
    'GWP'                 % 9. Global warming potential
};

% --- TOPSIS Analysis ---
% Extract objectives: Qp (7), LCOW (8), GWP (9)
objectives = All_data(:, [7, 8, 9]); % 80x3 matrix
obj_names = {'Qp', 'LCOW', 'GWP'};
n_solutions = size(objectives, 1); % 80
n_objectives = size(objectives, 2); % 3

% Define scenarios (weights for Qp, LCOW, GWP)
scenarios = {
    'Prioritise Qp', [0.8, 0.1, 0.1];
    'Prioritise LCOW', [0.0, 1, 0.0];
    'Prioritise GWP', [0.1, 0.1, 0.8];
    'Equal Weights', [0.333, 0.333, 0.333]
};

% Define optimisation direction: 1 for maximise, -1 for minimise
direction = [1, -1, -1]; % Maximise Qp, Minimise LCOW, Minimise GWP


% Step 1: Normalise the decision matrix
norm_matrix = zeros(n_solutions, n_objectives);
for j = 1:n_objectives
    col = objectives(:, j);
    norm_matrix(:, j) = col / norm(col); % Vector normalisation
end

% Initialise results storage
results = zeros(size(scenarios, 1), 5); % [Index, Qp, LCOW, GWP, Closeness]

% Process each scenario
for s = 1:size(scenarios, 1)
    weights = scenarios{s, 2};
    
    % Step 2: Apply weights to normalised matrix
    weighted_matrix = norm_matrix .* weights;
    
    % Step 3: Determine ideal and negative-ideal solutions
    ideal = zeros(1, n_objectives);
    neg_ideal = zeros(1, n_objectives);
    for j = 1:n_objectives
        if direction(j) == 1 % Maximise
            ideal(j) = max(weighted_matrix(:, j));
            neg_ideal(j) = min(weighted_matrix(:, j));
        else % Minimise
            ideal(j) = min(weighted_matrix(:, j));
            neg_ideal(j) = max(weighted_matrix(:, j));
        end
    end
    
    % Step 4: Calculate Euclidean distances
    dist_ideal = zeros(n_solutions, 1);
    dist_neg_ideal = zeros(n_solutions, 1);
    for i = 1:n_solutions
        dist_ideal(i) = sqrt(sum((weighted_matrix(i, :) - ideal).^2));
        dist_neg_ideal(i) = sqrt(sum((weighted_matrix(i, :) - neg_ideal).^2));
    end
    
    % Step 5: Calculate closeness coefficient
    closeness = dist_neg_ideal ./ (dist_ideal + dist_neg_ideal);
    
    % Step 6: Find top solution
    [~, top_idx] = max(closeness);
    results(s, :) = [top_idx, objectives(top_idx, :), closeness(top_idx)];
end

% Display results in a table
fprintf('TOPSIS Results for Four Scenarios\n');
fprintf('---------------------------------\n');
fprintf('Scenario\tIndex\tQp\tLCOW\tGWP\tCloseness\n');
for s = 1:size(scenarios, 1)
    fprintf('%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        scenarios{s, 1}, results(s, 1), results(s, 2:4), results(s, 5));
end

% Extract relevant data for the top solutions
topSolutions = All_data(results(:, 1), :);
% Display the top solutions for each scenario
fprintf('Top Solutions for Each Scenario:\n');
disp(array2table(topSolutions, 'VariableNames', var_names));

%% --- Parallel Coordinate Plot ---
% Create GroupData to highlight the four TOPSIS solutions
group_data = categorical(repmat({'Other'}, n_solutions, 1));
for s = 1:size(scenarios, 1)
    group_data(results(s, 1)) = scenarios{s, 1};
end

% Create parallel coordinate plot with raw data
figure('Name', 'Parallel Coordinate Plot for TOPSIS Solutions', 'Position', [100, 100, 1200, 600]);
p = parallelplot(All_data, 'GroupData', group_data);

% Customize plot
p.CoordinateTickLabels = var_names; % Set axis labels
p.LineStyle = '-'; % Solid line for all

% Set colors, line widths, and transparency for groups
% Note: categorical sorts groups alphabetically: ['Other','Prioritise LCOW','Prioritise Qp','Prioritise GWP','Equal Weights']
p.Color = [0.7 0.7 0.7;         % Gray for 'Other' 
           0.0 0.0 0.0;         % Black for 'Equal Weights'
           0.494, 0.184, 0.556; % Purple for 'Prioritise GWP'
           0, 0.447, 0.741;     % Orange for 'Prioritise Qp'
           0.85, 0.325, 0.098]; % Blue for 'Prioritise LCOW' 
p.LineWidth = [1, 2.5, 2.5, 2.5, 2.5]; % Line widths for ['Other','Prioritise LCOW','Prioritise Qp','Prioritise GWP','Equal Weights']
p.LineAlpha = [0.5, 1, 1, 1, 1]; % Transparency for ['Other','Prioritise LCOW','Prioritise Qp','Prioritise GWP','Equal Weights']

% To remove legend,
%p.LegendVisible = 'off';