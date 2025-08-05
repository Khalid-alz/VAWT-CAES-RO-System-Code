% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% The Levelised Cost of Water (LCOW) 
% 30/07/25

% This code calculates the Levelised Cost of Water (LCOW) for the water desalination system for VAWT-CAES-RO system.
% The code reads input data from a Excel file ("XXXXX.csv"), calculates the LCOW for each case in the file,
% and saves the results (including the calculated LCOW) to a new CSV file ("XXXXX_with_LCOW.csv").

% Input data includes:
%   - Nt: Number of vertical axis wind turbines (VAWTs)
%   - As: VAWTs Scale factor
%   - TR: Transmission ratio
%   - Vc: Compressor volume scale
%   - N_PV: Number of pressure vessels
%   - N_ro: Number of RO units
%   - total_Qp: Total water production per year
%   - Qc: brine water per hour
%   - D: Air compressor diameter
%   - L: Air compressor length

% The code calculates the following costs:
%   - WT_C: VAWTs cost
%   - AC_C: Air compressor cost
%   - Tank_C: Tank cost
%   - HP_pump_C: High-pressure pump cost
%   - LP_pump_C: Low-pressure pump cost
%   - RO_C: RO unit cost
%   - pre_post_C: Pre/post-treatment cost
%   - PX_C: Pressure exchanger cost
%   - II_C: Installation & infrastructure cost
%   - CAPEX: Capital cost of system
%   - ALC: Annual labor cost
%   - ACC: Annual chemical cost
%   - ARC: Annual replacement cost
%   - OPEX: Annual operation and maintenance cost
%   - LCOW: Levelised cost of water

% The LCOW is calculated using the following formula:
%   LCOW = (CAPEX * CFR + OPEX) / total_Qp
%   where CFR is the capital recovery factor.


%% Define file name

filename = 'ALL_VAWTs_CAES_RO_Cases.xlsx';

% Load the data from the Excel file for all cases
data = readtable(filename);

% Number of cases
numCases = size(data, 1);

% Initialize an array to store LOCW results
LCOW = zeros(numCases, 1);

% Loop through each case
for i = 1:numCases

    % Extract input variables for the current case
    Nt = data.Nt(i);        % Number of VAWTs
    As = data.As(i);        % WT Scale factor
    TR = data.TR(i);        % Transmission ratio
    Vc = data.Vc(i);        % Compressor volume scale
    N_PV = data.N_PV(i);    % Number of pressure vessels
    N_ro = data.N_ro(i);    % Number of RO units

    % Extract output variables for the current case
    total_Qp = data.total_Qp(i);  % Total water production per year
    Qc = data.Qc(i);              % brine water per hour

    %% wind turbine cost

    % The first column is the swept area in m^2, the second column is the power in kW, and the third column is the cost in USD.
    data_WT = [
        1.92, 0.3, 1800.00;
        3.2, 0.6, 2450.00;
        5.6, 1, 4070.00;
        7.04, 2, 4740.00;
        10.08, 3, 6600.00;
        21.6, 5, 12250.00;
        33, 10, 16770.00];

    SweptArea = data_WT(:, 1);
    power = data_WT(:, 2);
    Cost = data_WT(:, 3);

    p_WT_C = polyfit(SweptArea, Cost, 2); % Fit a Linear Model
    
    % Scaling calculations
    R = 0.515 * As;  % Turbine radius
    H = 1.456 * As;  % Turbine height
    A = R * 2 * H;   % Swept area

    WT_C = polyval(p_WT_C, A);

    WT_C = WT_C * Nt; % VAWTs cost

    %% Air compressor cost
    D = 0.059;         % (m) bore
    L = 0.05;          % (m) stroke

    Vs = 1 * (pi / 4) * D^2 * L;  % Base swept volume
    Vs = Vs * Vc;                 % Apply scaling factor
    Va_cm3 = Vs * 1000000;        % Induced Volume [cm3]
    
    AC_C = 223.40 * (Va_cm3 / 226.1)^0.7; % scaling factor method C2​=C1​(V2​/V1​​)^n
    AC_C = AC_C * Nt;                     % air compressor cost

    %% other system costs
    Tank_C = 3725; HP_pump_C = 11371; LP_pump_C = 735;
    RO_C = 842 * N_ro*N_PV;
    pre_post_C = (12/100) * (HP_pump_C + LP_pump_C + RO_C); % 12% of RO cost for pre/post-treatment units
    PX_C = 3134.7 * Qc^0.58;  % cost of the pressure exchanger

    %% Installation & infrastructure cost (30% of components costs) 
    II_C = (WT_C + AC_C + Tank_C + HP_pump_C + LP_pump_C + RO_C + pre_post_C + PX_C) * 0.3;

    %% The capital cost of system
    CAPEX = WT_C + AC_C + Tank_C + HP_pump_C + LP_pump_C + RO_C + pre_post_C + PX_C + II_C;

    %% annual operation and maintenance (O&M) costs

    % annual labour cost
    ALC = total_Qp * 0.05; % Annual Water Production [m^3]*operating labour cost (US$/m3)
    % annual  chemical cost
    ACC = total_Qp * 0.033; % Annual Water Production [m^3]*cost of chemicals (US$/m3)
    % annual replacement cost
    ARC = RO_C * (1/4) + HP_pump_C * (1/20) + LP_pump_C * (1/20) + AC_C * (1/3) ...
        + PX_C * (1/15) + Tank_C * (1/20); %+ WT_C * (1/20);

    OPEX = ALC + ACC + ARC;

    %% levelised cost of water

    r= 8/100; % discount rate(r)
    n= 25;    % system lifetime (n)
    CFR= ((1+r)^n)*r/(((1+r)^n)-1); % Capital recovery factor CRF [%/Year]

    LCOW(i) = (CAPEX * CFR + OPEX) / total_Qp;

end

% Add LCOW results to the table
data.LCOW = LCOW;

% Save the updated table to a new CSV file
writetable(data, 'ALL_VAWTs_CAES_RO_Cases_LCOW.xlsx');