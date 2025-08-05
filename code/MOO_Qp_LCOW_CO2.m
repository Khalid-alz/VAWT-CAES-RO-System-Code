% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% multi-objective optimisation
% 30/07/25

% This code to apply multi-objective optimisation using the three surrogate models
% and plotting all solutions alongside the final Pareto front solutions.

clear;clc;tic;

% Initialize global variables to store all solutions and objective values
global allSolutions allObjectives;
% Create empty arrays to hold solutions and objective values
allSolutions = [];
allObjectives = [];

% Create optimisation variables
Nt = optimvar("Nt","Type","integer","LowerBound",10,"UpperBound",80);
As = optimvar("As","LowerBound",2,"UpperBound",4.5);
TR = optimvar("TR","LowerBound",3,"UpperBound",10);
Vc = optimvar("Vc","LowerBound",1,"UpperBound",3);
N_PV = optimvar("N_PV","Type","integer","LowerBound",1,"UpperBound",2);
N_ro = optimvar("N_ro","Type","integer","LowerBound",2,"UpperBound",6);

% Set initial starting point for the solver
initialPoint.Nt = repmat(35,size(Nt));
initialPoint.As = repmat(3.5,size(As));
initialPoint.TR = repmat(6.5,size(TR));
initialPoint.Vc = repmat(1.7,size(Vc));
initialPoint.N_PV = repmat(2,size(N_PV));
initialPoint.N_ro = repmat(6,size(N_ro));

% Create problem
problem = optimproblem;

% Define problem objective
problem.Objective = fcn2optimexpr(@Qp_LCOW_CO2_objectives,Nt,As,TR,Vc,N_PV,...
    N_ro);

% Set nondefault solver options
options = optimoptions("gamultiobj",...
    'PopulationSize', 100, ... % Increase population size
    'MaxGenerations', 100, ... % Increase maximum generations
    'ParetoFraction', 0.35, ... % Increase Pareto fraction
    'CrossoverFraction', 0.8, ... % Increase crossover fraction
    'MutationFcn', @mutationadaptfeasible, ... % Use adaptive mutation
    'FunctionTolerance', 1e-4);  % Increase function tolerance

% Display problem information
show(problem);

% Solve problem
[solution,objectiveValue,reasonSolverStopped] = solve(problem,initialPoint,...
    "Solver","gamultiobj","Options",options);

% Display results
solution
reasonSolverStopped
objectiveValue

% Clear variables
clearvars Nt As TR Vc N_PV N_ro initialPoint options reasonSolverStopped...
    objectiveValue

% Display all solutions
All_data=[solution.Nt',solution.As',solution.TR',solution.Vc',solution.N_PV',solution.N_ro', -solution.Objective(1,:)', solution.Objective(2,:)', solution.Objective(3,:)'];

figure;
% Plot all solutions with filled gray circles
scatter3(-allObjectives(:,1), allObjectives(:,2), allObjectives(:,3), 15, 'o', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5], ... % gray fill
    'MarkerEdgeColor', [0, 0, 0], ...;      % Black border
    'MarkerFaceAlpha', 0.2, ...             % Transparency for fill
    'MarkerEdgeAlpha', 0.2)                 % Transparency for border
xlabel('Annual water production (Q_p, m^3/year)');
ylabel('Levelised cost of water (LCOW, US$/m³)');
zlabel('Global warming potential (GWP, kg CO₂eq/m³)');
grid on;
xlim([0 55000])
ylim([1 10])
zlim([0 6])
hAxes.XDir = "reverse";
hAxes.Box = "on";
hold on;
% Highlight the final Pareto front solutions with filled red circles
scatter3(-solution.Objective(1,:), solution.Objective(2,:), solution.Objective(3,:), 100, 'o', ...
    'MarkerFaceColor', [1, 0, 0], ... % Red fill
    'MarkerEdgeColor', [0, 0, 0]);    % Black border
legend('All Solutions','Final Pareto Front');
hold off;

% Create a new figure for the 2D scatter plots
figure;
% Plot 2D scatter Qp vs LCOW
subplot(2, 1, 1); % Create a subplot in a 2-row, 1-column layout
scatter(-solution.Objective(1,:), solution.Objective(2,:), 'o', ...
    'MarkerFaceColor', [1, 0, 0], ... % Red fill
    'MarkerEdgeColor', [0, 0, 0]);    % Black border
xlabel('Annual water production (Q_p, m^3/year)');
ylabel('Levelised cost of water (LCOW, US$/m³)');
grid on;
% Plot 2D scatter Qp vs GWP
subplot(2, 1, 2); % Create a subplot in a 2-row, 1-column layout
scatter(-solution.Objective(1,:), solution.Objective(3,:), 'o', ...
    'MarkerFaceColor', [1, 0, 0], ... % Red fill
    'MarkerEdgeColor', [0, 0, 0]);    % Black border
xlabel('Annual water production (Q_p, m^3/year)');
ylabel('Global warming potential (GWP, kg CO₂eq/m³)');
grid on;

toc
%% Helper Functions

function  outputs = Qp_LCOW_CO2_objectives(Nt, As, TR, Vc, N_PV, N_ro)

% Load models only if they have not been loaded before
persistent Qp_MODEL;
if isempty(Qp_MODEL)
    load O_GPR_Qp_10F.mat; % Load the Qp model from file
    Qp_MODEL = O_GPR_Qp_10F;
end

persistent LCOW_MODEL;
if isempty(LCOW_MODEL)
    load ONN_LCOW_10F.mat; % Load the LCOW model from file
    LCOW_MODEL = ONN_LCOW_10F;
end

persistent CO2_MODEL;
if isempty(CO2_MODEL)
    load ONN_CO2_10F.mat; % Load the CO2 model from file
    CO2_MODEL = ONN_CO2_10F_V2;
end

% Make predictions using the loaded models
P_Qp = Qp_MODEL.predictFcn([Nt, As, TR, Vc, N_PV, N_ro]);
P_LCOW = LCOW_MODEL.predictFcn([Nt, As, TR, Vc, N_PV, N_ro]);
P_CO2 = CO2_MODEL.predictFcn([Nt, As, TR, Vc, N_PV, N_ro]);

% Store the solution and objective values
global allSolutions allObjectives;
allSolutions = [allSolutions; Nt, As, TR, Vc, N_PV, N_ro]; % Store the rounded solution
allObjectives = [allObjectives; [-P_Qp, P_LCOW, P_CO2]];   % Ensure this is a row vector

% Combine output variables into a single array
outputs = [-P_Qp P_LCOW P_CO2];
end