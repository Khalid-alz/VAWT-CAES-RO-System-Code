% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% comparing LCOW model against real/test data
% 30/07/25

% This code compares the surrogate model of LCOW against real/test data
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

random_inputs = [random_Nt, random_As, random_TR, random_Vc, random_N_PV, random_N_ro];
num_random_points = length(Qp);
Qp = Qp(1:num_random_points);
LCOW = LCOW(1:num_random_points);

%% Load surrogate model of LCOW
load('ONN_LCOW_10F.mat'); % Load your ML model LCOW 

%% Evaluate Models
ml1_LCOW = zeros(num_random_points, 1);

parfor i = 1:num_random_points
    ml1_LCOW(i) = ONN_LCOW_10F.predictFcn(random_inputs(i, :));
end

%% Calculate Errors and Metrics
error_ml1 = LCOW - ml1_LCOW;
rmse_ml1 = sqrt(mean(error_ml1.^2, 'omitnan'));
disp(['RMSE surrogate model of LCOW: ', num2str(rmse_ml1)]);

%% Plotting
% Predicted vs. Actual
figure;
scatter(LCOW, ml1_LCOW, 'DisplayName', 'surrogate model of LCOW');
hold on;
plot([min(LCOW), max(LCOW)], [min(LCOW), max(LCOW)], 'DisplayName', 'y=x'); % Diagonal line
xlabel('Actual LCOW');
ylabel('Predicted LCOW');
title('Predicted vs. Actual LCOW');
legend('show');
hold off;

% Error Distribution
figure;
histogram(error_ml1, 'Normalization', 'probability', 'DisplayName', 'surrogate model of LCOW');
xlabel('Error (Actual - Predicted)');
ylabel('Probability');
title('Error Distribution');
legend('show');

% Error Box Plot
figure;
boxplot(error_ml1, {'surrogate model of LCOW'});
ylabel('Error (Actual - Predicted)');
title('Error Comparison (Box Plot)');

% the real values plot for the Predictions vs. Actual 

% Sorted Predictions vs. Actual
[sorted_simulink_LCOW, sort_indices] = sort(LCOW);
sorted_ml1 = ml1_LCOW(sort_indices);
figure;
plot(sorted_simulink_LCOW, 'w-', 'DisplayName', 'Simulink');
hold on;
plot(sorted_ml1, 'DisplayName', 'surrogate model of LCOW');
xlabel('Sorted Cases');
ylabel('LCOW');
title('Sorted Predictions vs. Actual');
legend('show');
grid on
hold off;