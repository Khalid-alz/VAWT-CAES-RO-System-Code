% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% Run all cases from the Simulink model
% 30/07/25

% This code runs the Simulink model of the VAWT-CASE-RO system in parallel computing
% for all cases from the generated grid or from a file.

%% start the pararall pool
clc; clear;tic;
parpool(8) % Change the number of required cores

%% Define Ranges for Design Variables
[Nt_grid, As_grid, TR_grid, Vc_grid, N_ro_grid, N_PV_grid] = ndgrid(...
    10:5:80, ...    % Nt: Number of vertical axis wind turbines (VAWTs)
    2:0.25:4.5, ... % As: VAWTs Scale factor
    3:10, ...       % TR: Transmission ratio
    1:0.25:3, ...   % Vc: Compressor volume scale
    2:6, ...        % N_ro: Number of RO units
    1:2);           % N_PV: Number of pressure vessels

% Reshape grids into column vectors
Nt_values = Nt_grid(:);
As_values = As_grid(:);
TR_values = TR_grid(:);
Vc_values = Vc_grid(:);
N_ro_values = N_ro_grid(:);
N_PV_values = N_PV_grid(:);

% Total cases
num_cases = numel(Nt_values);
disp(['Total Cases: ', num2str(num_cases)]);

%% Define Start and End Indices (Example: Run first 50000 cases)
start_index = 1;
end_index = 3; % Or any other number less than num_cases

%calculate the number of cases that will be run.
run_cases = end_index - start_index +1;
disp(['Running cases ', num2str(start_index), ' to ', num2str(end_index)]);

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

%% Set Up Result File
results_filename = 'Results.csv';

% Create file if it does not exist
if ~isfile(results_filename)
    headers = {'Nt', 'As', 'TR', 'Vc', 'N_PV', 'N_ro', 'total_Qp', 'total_RR', 'total_Qf', 'total_Qc', 'Qf', 'Qp', 'Qc'};
    writetable(cell2table(cell(0, numel(headers)), 'VariableNames', headers), results_filename);
end

%% Load Simulink Model
model_system = 'Sim_main';

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
simIn = repmat(baseSimIn, run_cases, 1);

for i = start_index:end_index

    disp(['Running case ', num2str(i), ' of ', num2str(end_index)]);

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
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('Nt', Nt);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('As', As);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('TR', TR);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('Vc', Vc);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('N_PV', N_PV);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('N_ro', N_ro);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('R', R);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('A', A);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('J1', J1);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('Va', Va);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('Tor_rq', Tor_rq);
    simIn(i - start_index + 1) = simIn(i - start_index + 1).setVariable('Limit_Tor', Limit_Tor);
end

%% Run Simulations in Parallel
disp('Starting parallel simulations...');
simOut = parsim(simIn, 'ShowProgress', 'on');

%% Extract & Save Results
Results = zeros(run_cases, 13);

for i = start_index:end_index

    % Extract simulation results
    Qf = max(simOut(i - start_index + 1).logsout.get('Qf').Values.Data);
    Qc = max(simOut(i - start_index + 1).logsout.get('Qc').Values.Data);
    Qp = max(simOut(i - start_index + 1).logsout.get('Qp').Values.Data);
    total_Qf = simOut(i - start_index + 1).logsout.get('total_Qf').Values.Data(end);
    total_Qc = simOut(i - start_index + 1).logsout.get('total_Qc').Values.Data(end);
    total_Qp = simOut(i - start_index + 1).logsout.get('total_Qp').Values.Data(end);
    total_RR = max(simOut(i - start_index + 1).logsout.get('total_RR').Values.Data);

    % Store results
    Results(i - start_index + 1, :) = [Nt_values(i), As_values(i), TR_values(i), Vc_values(i), ...
        N_PV_values(i), N_ro_values(i), total_Qp, total_RR, ...
        total_Qf, total_Qc, Qf, Qp, Qc];

end

% Convert results to table and save
ResultsTable = array2table(Results, ...
    'VariableNames', {'Nt', 'As', 'TR', 'Vc', 'N_PV', 'N_ro', 'total_Qp', 'total_RR', ...
    'total_Qf', 'total_Qc', 'Qf', 'Qp', 'Qc'});
writetable(ResultsTable, results_filename, 'WriteMode', 'Append');

time = toc;
disp(['Simulation completed. Results saved. Time = ', num2str(time), ' sec']);
