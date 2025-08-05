% Khalid Alzahrani
% k.m.alzahrani@tu.edu.sa
% VAWT-CASE-RO system
% comparing Qp model against real/test data
% 30/07/25

% This code compares the surrogate model of Qp against real/test data
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

%% Load surrogate model of Qp
load('O_GPR_Qp_10F.mat'); % Load your ML model

%% Evaluate Models
ml1_Qp = zeros(num_random_points, 1);

parfor i = 1:num_random_points
    ml1_Qp(i) = O_GPR_Qp_10F.predictFcn(random_inputs(i, :)); % be sure the predictor name
end

%% Calculate Errors and Metrics
error_ml1 = Qp - ml1_Qp;
rmse_ml1 = sqrt(mean(error_ml1.^2, 'omitnan'));
disp(['RMSE surrogate model of Qp: ', num2str(rmse_ml1)]);

%% Plotting
% Predicted vs. Actual
figure;
scatter(Qp, ml1_Qp, 'DisplayName', 'surrogate model of Qp');
hold on;
plot([min(Qp), max(Qp)], [min(Qp), max(Qp)], 'DisplayName', 'y=x'); % Diagonal line
xlabel('Actual total\_Qp');
ylabel('Predicted total\_Qp');
title('Predicted vs. Actual total\_Qp');
legend('show');
hold off;

% Error Distribution
figure;
histogram(error_ml1, 'Normalization', 'probability', 'DisplayName', 'surrogate model of Qp');
xlabel('Error (Actual - Predicted)');
ylabel('Probability');
title('Error Distribution');
legend('show');

% Error Box Plot
figure;
boxplot(error_ml1, {'surrogate model of Qp'});
ylabel('Error (Actual - Predicted)');
title('Error Comparison (Box Plot)');

% the real values plot for the Predictions vs. Actual 

% Sorted Predictions vs. Actual
[Qp, sort_indices] = sort(Qp);
sorted_ml1 = ml1_Qp(sort_indices);
figure;
plot(Qp, 'w-', 'DisplayName', 'Actual');
hold on;
plot(sorted_ml1, 'DisplayName', 'surrogate model of Qp');
xlabel('Sorted Cases');
ylabel('total\_Qp');
title('Sorted Predictions vs. Actual');
legend('show');
grid on
hold off;