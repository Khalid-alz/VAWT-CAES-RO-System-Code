% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% comparing GWP model against real/test data
% 30/07/25

% This code compares the surrogate model of GWP (CO2) against real/test data
% by calculating RMSE and plotting the results to compare the surrogate model against the real data. 
% Change the name of the test data, and be sure the name variables are the same as the ones in the file.

%% Load real/test dat from Simulink results

% Load the test data from the Excel file
results_table = readtable('simulink_random_LCOW_CO2_filtered.xlsx');
random_Nt = results_table.Nt;
random_As = results_table.As;
random_TR = results_table.TR;
random_Vc = results_table.Vc;
random_N_PV = results_table.N_PV;
random_N_ro = results_table.N_ro;
Qp = results_table.total_Qp;
LCOW = results_table.LCOW;
CO2 = results_table.KgCO2_m3;

random_inputs = [random_Nt, random_As, random_TR, random_Vc, random_N_PV, random_N_ro];
num_random_points = length(CO2);
CO2 = CO2(1:num_random_points);

%% Load surrogate model of CO2
load('ONN_CO2_10F.mat'); % Load your ML model CO2 

%% Evaluate Models
ml1_CO2 = zeros(num_random_points, 1);

parfor i = 1:num_random_points
    ml1_CO2(i) = ONN_CO2_10F_V2.predictFcn(random_inputs(i, :));
end

%% Calculate Errors and Metrics
error_ml1 = CO2 - ml1_CO2;
rmse_ml1 = sqrt(mean(error_ml1.^2, 'omitnan'));
disp(['RMSE surrogate model of GWP: ', num2str(rmse_ml1)]);

%% Plotting
% Predicted vs. Actual
figure;
scatter(CO2, ml1_CO2, 'DisplayName', 'surrogate model of GWP');
hold on;
plot([min(CO2), max(CO2)], [min(CO2), max(CO2)]); % Diagonal line
xlabel('Actual GWP');
ylabel('Predicted GWP');
title('Predicted vs. Actual GWP');
legend('show');
hold off;

% Error Distribution
figure;
histogram(error_ml1, 'Normalization', 'probability', 'DisplayName', 'ML1 5F');
xlabel('Error (Actual - Predicted)');
ylabel('Probability');
title('Error Distribution');
legend('show');
hold off;

% Error Box Plot
figure;
boxplot(error_ml1, {'surrogate model of GWP'});
ylabel('Error (Actual - Predicted)');
title('Error Comparison (Box Plot)');

% the real values plot for the Predictions vs. Actual 

% Sorted Predictions vs. Actual
[sorted_simulink_CO2, sort_indices] = sort(CO2);
sorted_ml1 = ml1_CO2(sort_indices);
figure;
plot(sorted_simulink_CO2, 'w-', 'DisplayName', 'Actual');
hold on;
plot(sorted_ml1, 'DisplayName', 'surrogate model of GWP');
xlabel('Sorted Cases');
ylabel('GWP');
title('Sorted Predictions vs. Actual');
legend('show');
grid on
hold off;