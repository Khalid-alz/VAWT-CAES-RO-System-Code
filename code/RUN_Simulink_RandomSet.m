% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% Latin Hypercube Samples for the test set for Simulink model
% 30/07/25

% This code runs the Simulink model of the VAWT-CASE-RO system in parallel computing
% for test set (random using Latin Hypercube Samples) cases from the generated data or from a file.

%% start the pararall pool
clc; clear;tic;
%parpool(8) % Change the number of required cores

%% Define Ranges for Design Variables
Nt_range = [10, 80];
As_range = [2, 4.5];
TR_range = [3, 10];
Vc_range = [1, 3];
N_ro_range = [2, 6];
N_PV_range = [1, 2];

%% Generate Latin Hypercube Samples

num_random_points = 10000; % Adjust as needed

% Define the ranges for the parameters
ranges = [Nt_range; As_range; TR_range; Vc_range; N_ro_range; N_PV_range];

% Generate Latin Hypercube Samples for the defined ranges
lhs_samples = lhsdesign(num_random_points, size(ranges, 1));

% Scale the samples to the defined ranges for each parameter
Nt_values = round(ranges(1, 1) + lhs_samples(:, 1) * (ranges(1, 2) - ranges(1, 1)));
As_values = ranges(2, 1) + lhs_samples(:, 2) * (ranges(2, 2) - ranges(2, 1));
TR_values = ranges(3, 1) + lhs_samples(:, 3) * (ranges(3, 2) - ranges(3, 1));
Vc_values = ranges(4, 1) + lhs_samples(:, 4) * (ranges(4, 2) - ranges(4, 1));
N_ro_values = round(ranges(5, 1) + lhs_samples(:, 5) * (ranges(5, 2) - ranges(5, 1)));
N_PV_values = round(ranges(6, 1) + lhs_samples(:, 6) * (ranges(6, 2) - ranges(6, 1)));

% Combine the generated values into a single matrix of random inputs
random_inputs = [Nt_values, As_values, TR_values, Vc_values, N_PV_values, N_ro_values];

%% Air compressor parameters
D = 0.059;         % (m) bore
L = 0.05;          % (m) stroke
Pa = 101325;       % (Pa) atmosperic pressure
P1 = Pa;           % (Pa) inlet pressure
P2 = 800000;       % (Pa) outlet pressure
T1 = 298;          % (K) inlet temp
n = 1.25;          % polytrophic exponent
C = 0.07;          % clearance volume
Vt = 2000;         % (litre) valume for the air compressed tank
CutIn = 3;         % Cut in speed for VAWT

% pre-calculations
Vs = 1*(pi/4)*D^2*L;                         % Swept volume (Displacement Volume) (m3)
e_vol = 1 + C - C*(P2/P1)^(1/n);             % volumetric efficiency
T2 = T1*(P2/P1)^((n-1)/n);                   % Delivery temperature
Vf1= ((P1+Pa)*T2/(T1*(P2+Pa)))*1000;         % FAD air flow rate [L/s]
IP1 = (n/(n-1))*P1*( (P2/P1)^((n-1)/n) - 1); % Indicated power of the compressor, Since the compressor is single acting,therefore number of working strokes per minute = Nw = N

%% RO parameters
Cf= 36000;  % Feed Concentration [mg/L]

% Load all trained models for Qp, Cp and Pc
load('RO_Qp.mat');
load('RO_Cp.mat');
load('RO_Pc.mat');
%% Load Simulink Model and Run Simulations
model_system = 'Sim_main';

%simIn = repmat(Simulink.SimulationInput(model_system), num_random_points, 1);

% Assign Global Variables Only Once (to reduce redundancy)
baseSimIn = Simulink.SimulationInput(model_system);
baseSimIn = baseSimIn.setVariable('RO_Qp', RO_Qp);
baseSimIn = baseSimIn.setVariable('RO_Cp', RO_Cp);
baseSimIn = baseSimIn.setVariable('RO_Pc', RO_Pc);
baseSimIn = baseSimIn.setVariable('Vt', Vt);
baseSimIn = baseSimIn.setVariable('CutIn', CutIn);
baseSimIn = baseSimIn.setVariable('P2', P2);
baseSimIn = baseSimIn.setVariable('IP1', IP1);
baseSimIn = baseSimIn.setVariable('Cf', Cf);
baseSimIn = baseSimIn.setVariable('Vf1', Vf1);
baseSimIn = baseSimIn.setModelParameter('StopTime', '8760');

% Preallocate Simulink Input Array
simIn = repmat(baseSimIn, num_random_points, 1);

for i = 1:num_random_points

    disp(['Running case ', num2str(i), ' of ', num2str(num_random_points)]);
    
    % Assign design variables
    Nt = Nt_values(i);
    As = As_values(i);
    TR = TR_values(i);
    Vc = Vc_values(i);
    N_ro = N_ro_values(i);
    N_PV = N_PV_values(i);

    % Scaling calculations
    R = 0.515 * As;  % Turbine radius
    H = 1.456 * As;  % Turbine height
    A = R * 2 * H;   % Swept area
    J1 = 1.5 * As;   % Turbine moment of inertia

    % Air compressor displacement volume scaling
    Vs = Vs * Vc;                       % Apply scaling factor
    Va = e_vol * Vs;                    % Induced Volume
    Tor_rq = IP1 * Va * 10 / (17 * pi); % Torque required
    Limit_Tor = Tor_rq * 0.5;           % Torque limit

    % Configure SimulationInput object
    simIn(i) = simIn(i).setVariable('Nt', Nt);
    simIn(i) = simIn(i).setVariable('As', As);
    simIn(i) = simIn(i).setVariable('TR', TR);
    simIn(i) = simIn(i).setVariable('Vc', Vc);
    simIn(i) = simIn(i).setVariable('N_PV', N_PV);
    simIn(i) = simIn(i).setVariable('N_ro', N_ro);
    simIn(i) = simIn(i).setVariable('R', R);
    simIn(i) = simIn(i).setVariable('A', A);
    simIn(i) = simIn(i).setVariable('J1', J1);
    simIn(i) = simIn(i).setVariable('Va', Va);
    simIn(i) = simIn(i).setVariable('Tor_rq', Tor_rq);
    simIn(i) = simIn(i).setVariable('Limit_Tor', Limit_Tor);
end

%% Run Simulations in Parallel
disp('Starting parallel simulations...');
simOut = parsim(simIn, 'ShowProgress', 'on');

%% Extract & Save Results
total_Qp = zeros(num_random_points, 1);
Qc= zeros(num_random_points, 1);
total_Qf= zeros(num_random_points, 1);

for i = 1:num_random_points
    % Extract simulation results
    total_Qp(i) = simOut(i).logsout.get('total_Qp').Values.Data(end);
    Qc(i) = max(simOut(i).logsout.get('Qc').Values.Data);
    total_Qf(i) = simOut(i).logsout.get('total_Qf').Values.Data(end);
end

%% Save Results to File

Nt = Nt_values;
As = As_values;
TR = TR_values;
Vc = Vc_values;
N_ro = N_ro_values;
N_PV = N_PV_values;

% Convert results to table and save
results_table = table(Nt, As, TR, Vc, N_PV, N_ro, total_Qp, Qc, total_Qf);
writetable(results_table, 'simulink_LHS_random_results.csv');

time = toc;
disp(['Simulink simulations completed. Results saved to simulink_random_resultsV3.csv. Time = ', num2str(time), ' sec']);