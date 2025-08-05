% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% single case to run Simulink
% 30/07/25

% This code to run the simlunik for a single case where RUN_Simulink_POOL.m for all cases
clc; clear;
%% inputs for the simulink model

%  input variables for the Equal Weights (33.3%) scenario
Nt = 58;     % Number of VAWTs
As = 2.92;   % WT Scale factor
TR = 9.36;   % Transmission ratio
Vc = 1.48;   % Compressor volume scale
N_PV = 2;    % Number of pressure vessels
N_ro = 6;    % Number of RO units

%%  VAWT calculations after scaling
R = 0.515 * As;  % Turbine radius
H = 1.456 * As;  % Turbine height
A = R * 2 * H;   % Swept area
J1 = 1.5 * As;   % Turbine moment of inertia
CutIn = 3;       % Cut in speed for VAWT
%% Compressor design and tank

% input parameters
D = 0.059;         % 0.063(m) bore
L = 0.05;          % 0.06 (m) stroke
Pa = 101325;       % (Pa) atmosperic pressure
P1 = Pa;           % (Pa) inlet pressure
P2 = 800000;       % (Pa) outlet pressure
T1 = 298;          % (K) inlet temp
n = 1.25;          % polytrophic exponent
C = 0.07;          % clearance volume
Vt = 2000;         % (litre) valume for the air compressed tank

% pre-calculations
Vs = 1*(pi/4)*D^2*L;                         % Swept volume (Displacement Volume) (m3)
e_vol = 1 + C - C*(P2/P1)^(1/n);             % volumetric efficiency 
T2 = T1*(P2/P1)^((n-1)/n);                   % Delivery temperature
Vf1= ((P1+Pa)*T2/(T1*(P2+Pa)))*1000;         % FAD air flow rate [L/s]
IP1 = (n/(n-1))*P1*( (P2/P1)^((n-1)/n) - 1); % Indicated power of the compressor, Since the compressor is single acting,therefore number of working strokes per minute = Nw = N

% Air compressor displacement volume scaling
Vs = Vs * Vc;                       % Apply scaling factor
Va = e_vol * Vs;                    % Induced Volume (m3)
Tor_rq = IP1 * Va * 10 / (17 * pi); % Torque required
Limit_Tor = Tor_rq * 0.5;           % Torque limit

%% RO calculation
Cf= 36000;  % Feed Concentration [mg/L]

% Load all trained models
load('RO_Qp.mat');
load('RO_Cp.mat');
load('RO_Pc.mat');

%% Run the Simulink model

Output = sim("Sim_main.slx", 'StopTime', '8760');