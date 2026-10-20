% --- MAIN SCRIPT ---
clear; close all;

% 1. Define the target percentages based on your clinical data
% Bins: 0-4, 5-9, 10-20, >20
targets = [16.657, 50.903, 15.259, 17.181]';

% 2. Setup Genetic Algorithm bounds and options
nvars = 2; % Optimizing two parameters: [mu, sigma]

% Define Lower Bounds (lb) and Upper Bounds (ub)
% Because of the formula X = mu * exp(randn * sigma), mu is now the median on the linear scale.
% Therefore, mu MUST be > 0.
lb = [1e-4, 1e-4]; 
ub = [100, 5]; 

% Set GA options for display and population size
options = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 54, ...
    'UseParallel', true, 'PlotFcn', {@gaplotbestf, @gaplotdistance, @gaplotscores});

% 3. Run the optimizer (Genetic Algorithm)
fprintf('Starting Genetic Algorithm optimization for Custom Log-Normal...\n');
% Reset global RNG so the GA's internal random operations are reproducible
rng('default'); 
best_params = ga(@(params) sim_objfun(params, targets), nvars, [], [], [], [], lb, ub, [], options);

% Extract and display the best parameters
mu_fit = best_params(1);
sigma_fit = best_params(2);
fprintf('\nOptimized Custom Log-Normal parameters: mu = %.4f, sigma = %.4f\n', mu_fit, sigma_fit);

% 4. Verify the fit by generating a final set of virtual patients
% Use the custom formula and round the output to mimic discrete count data
virtual_patients = round(mu_fit .* exp(randn(10000, 1) .* sigma_fit));
final_pct = [
    sum(virtual_patients == 0);
    sum(virtual_patients > 0 & virtual_patients <= 9);
    sum(virtual_patients >= 10 & virtual_patients <= 20);
    sum(virtual_patients > 20)
] / 10000 * 100;

fprintf('\nTarget Percentages vs. Generated Percentages:\n');
disp(table(targets, final_pct, 'VariableNames', {'Target_Pct', 'Generated_Pct'}));


% 6. Visualization: Binned Clinical Targets vs. Simulated Data
figure;

% Combine the target percentages and generated percentages into a 4x2 matrix
comparison_data = [targets, final_pct]; 

% Define the bin labels for the x-axis
bin_labels = {'0-4', '5-9', '10-20', '> 20'};
x_categorical = categorical(bin_labels, bin_labels);

% Create a grouped bar chart
b = bar(x_categorical, comparison_data, 'grouped');

% Customize colors
b(1).FaceColor = [0.2 0.4 0.7]; % Clinical Data - Dark Blue
b(2).FaceColor = [0.8 0.4 0.2]; % Simulated Data - Dark Orange

% Add the exact percentage values on top of the bars for quick reading
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1, ytips1 + 1, string(round(b(1).YData, 1)) + "%", ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
text(xtips2, ytips2 + 1, string(round(b(2).YData, 1)) + "%", ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

% Formatting the plot
title('Clinical Targets vs. Simulated Virtual Patients (Custom Log-Normal)');
xlabel('Number of Bone Metastases (Bins)');
ylabel('Percentage of Patients (%)');
legend('Clinical Data (Target)', 'Simulated Cohort (Generated)', 'Location', 'northeast');
ylim([0, max(max(comparison_data)) + 10]); 
grid on; pbaspect([4 2 1]);

% --- OBJECTIVE FUNCTION ---
function sse = sim_objfun(params, targets)
    mu = params(1);
    sigma = params(2);
    
    % CRITICAL: Lock the random seed so the optimizer sees a stable error surface
    rng(0); 
    
    % Generate samples using the custom formula X = mu * exp(Z * sigma)
    num_samples = 100000;
    
    % Generate standard normal values, multiply by sigma, exponentiate, scale by mu, and round
    samples = round(mu .* exp(randn(num_samples, 1) .* sigma));
    
    % Count frequencies for each bin
    bin1 = sum(samples == 0);
    bin2 = sum(samples > 0 & samples <= 9);
    bin3 = sum(samples >= 10 & samples <= 20);
    bin4 = sum(samples > 20);
    
    % Calculate percentages
    pct = [bin1; bin2; bin3; bin4] / num_samples * 100;
    
    % Calculate Chi-Square style Sum of Squared Errors
    sse = sum((pct - targets).^2 ./ targets);
end