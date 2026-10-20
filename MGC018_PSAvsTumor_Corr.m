clear; clc; close all;

% --- 1. Define Data ---
% Structure: .t_psa (time), .v_psa (values), .t_les (time), .v_les (values)

% Patient 1
p(1).t_psa = [0; 4.26; 8.03];
p(1).v_psa = [0; 61.79; 157.72];
p(1).t_les = [0; 9.09];
p(1).v_les = [0; 30.49];

% Patient 2
p(2).t_psa = [0; 3.05; 5.91; 7.39; 9.00; 12.12; 14.96; 16.06; 18.94];
p(2).v_psa = [0; -21.13; -59.35; -67.48; -71.54; -75.61; -67.48; -58.54; -62.60];
p(2).t_les = [0; 10; 18.09; 27.12];
p(2).v_les = [0; -13.19; -13.19; 4.95];

% Patient 3
p(3).t_psa = [0; 4.00; 9.07; 14.06];
p(3).v_psa = [0; -84.80; -98.26; -98.26];
p(3).t_les = [0; 9.00; 15.38];
p(3).v_les = [0; -100; -100];

% --- 2. Interpolation (Inverse Strategy) ---
matched_psa = [];
matched_lesion = [];
patient_labels = []; 

for i = 1:length(p)
    % We now iterate through PSA times (dense) and interpolate Lesion (sparse)
    target_times = p(i).t_psa;
    
    % Interpolate Lesion values at PSA times
    % 'linear' is standard. 'extrap' handles cases where PSA time is slightly
    % outside Lesion bounds (though your data is mostly well-bounded).
    est_lesion = interp1(p(i).t_les, p(i).v_les, target_times, 'linear', 'extrap');
    
    % Store the raw PSA values and the ESTIMATED Lesion values
    matched_psa = [matched_psa; p(i).v_psa];
    matched_lesion = [matched_lesion; est_lesion];
    
    % Create label vector for plotting
    patient_labels = [patient_labels; repmat(i, length(target_times), 1)];
end

% --- 3. Filter to 14 Points (Consolidate Baselines) ---
% Find indices where both values are exactly 0 (The Baselines)
baseline_indices = find(matched_psa == 0 & matched_lesion == 0);

% Remove all baselines first
matched_psa(baseline_indices) = [];
matched_lesion(baseline_indices) = [];
patient_labels(baseline_indices) = [];

% Add back exactly ONE baseline point (0,0)
matched_psa = [0; matched_psa];
matched_lesion = [0; matched_lesion];
patient_labels = [0; patient_labels]; % Label 0 for the shared origin

% --- 4. Calculate Correlation ---
[r_pearson, p_pearson] = corr(matched_psa, matched_lesion, 'Type', 'Pearson');
[r_spearman, p_spearman] = corr(matched_psa, matched_lesion, 'Type', 'Spearman');

% --- 5. Output & Visualization ---
fprintf('\n----------------------------------------\n');
fprintf('Analysis Results (Matched Pairs: n=%d)\n', length(matched_psa));
fprintf('----------------------------------------\n');
fprintf('Pearson Correlation (r):  %.4f  (p-value: %.4f)\n', r_pearson, p_pearson);
fprintf('Spearman Correlation (rho): %.4f  (p-value: %.4f)\n', r_spearman, p_spearman);

figure('Color','w');
hold on; box on;

% Define Colors (0=Black for Baseline, 1=Red, 2=Green, 3=Blue)
cmap = [0 0 0; 1 0 0; 0 0.5 0; 0 0 1]; 

% Scatter Plot
for i = 0:3
    idx = (patient_labels == i);
    if sum(idx) > 0
        scatter(matched_psa(idx), matched_lesion(idx), 80, cmap(i+1,:), 'filled', ...
            'DisplayName', sprintf('Patient %d', i));
    end
end

% Regression Line
p_fit = polyfit(matched_psa, matched_lesion, 1);
x_fit = linspace(min(matched_psa), max(matched_psa), 100);
y_fit = polyval(p_fit, x_fit);
plot(x_fit, y_fit, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Regression');

title({'Correlation: PSA vs Interpolated Lesion Size', ...
       '(Interpolating Lesion at PSA timepoints)'});
xlabel('PSA Change (%)','FontWeight','bold');
ylabel('Estimated Lesion Change (%)','FontWeight','bold');
legend({'Shared Baseline', 'Patient 1', 'Patient 2', 'Patient 3', 'Fit Line'}, 'Location', 'Best');
xlim([-105, 200]);
ylim([-105, 70]);
grid off;
box on;
pbaspect([3 2 1]);
fontsize(15,'points');
fontname('Consolas');


% Define a fixed color map for consistency
% Patient 1 = Red, Patient 2 = Green, Patient 3 = Blue
patient_colors = [0.85 0.32 0.10;  % Red/Orange
                  0.47 0.67 0.19;  % Green
                  0.00 0.45 0.74]; % Blue

% ---------------------------------------------------------
% VISUALIZATION 1: Side-by-Side Comparison
% ---------------------------------------------------------
figure('Name', 'Side-by-Side Trajectories', 'Color', 'w');

% -- Subplot 1: PSA Trajectories --
subplot(1, 2, 1); hold on; box on; grid off;
for i = 1:3
    plot(p(i).t_psa, p(i).v_psa, 'o-', ...
        'Color', patient_colors(i,:), ...
        'LineWidth', 2, ...
        'MarkerFaceColor', patient_colors(i,:), ...
        'DisplayName', sprintf('Pt %d', i));
end
title('PSA Trajectories');
xlabel('Time (Weeks)','FontWeight','bold');
ylabel('% Change from Baseline','FontWeight','bold');
fontsize(15,"points");
legend('show', 'Location', 'northeast');
pbaspect([1 1 1]);
ylim([-110 200]); % Set limits to make comparison easier

% -- Subplot 2: Lesion Trajectories --
subplot(1, 2, 2); hold on; box on;
for i = 1:3
    plot(p(i).t_les, p(i).v_les, 's-', ... % Square marker for lesion
        'Color', patient_colors(i,:), ...
        'LineWidth', 2, ...
        'MarkerFaceColor', patient_colors(i,:), ...
        'DisplayName', sprintf('Pt %d', i));
end
title('Lesion Size Trajectories');
xlabel('Time (Weeks)','FontWeight','bold');
ylabel('% Change from Baseline','FontWeight','bold');
fontsize(15,"points");
legend('show', 'Location', 'northeast');
pbaspect([1 1 1]);
ylim([-110 200]); % Match Y-limits with PSA plot for fair comparison
fontname('Consolas');

% ---------------------------------------------------------
% VISUALIZATION 2: Combined Plot (All 6 Curves)
% ---------------------------------------------------------
figure('Name', 'Combined Trajectories', 'Color', 'w');
hold on; 
% Dummy plots for Legend creation (to explain Line Styles)
h1 = plot(nan, nan, 'k-', 'LineWidth', 2, 'DisplayName', 'PSA (Solid)');
h2 = plot(nan, nan, 'k--', 'LineWidth', 2, 'DisplayName', 'Lesion (Dashed)');

for i = 1:3
    % Plot PSA (Solid Line)
    plot(p(i).t_psa, p(i).v_psa, 'o-', ...
        'Color', patient_colors(i,:), ...
        'LineWidth', 2, ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', patient_colors(i,:), ...
        'HandleVisibility', 'off'); % Hide from legend to prevent clutter
    
    % Plot Lesion (Dashed Line)
    plot(p(i).t_les, p(i).v_les, 's--', ...
        'Color', patient_colors(i,:), ...
        'LineWidth', 2, ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', 'w', ... % Hollow marker for contrast
        'HandleVisibility', 'off'); % Hide from legend
end


title('Combined Dynamics: PSA vs Lesion Size');
xlabel('Time (Weeks)','FontWeight','bold');
ylabel('% Change from Baseline','FontWeight','bold');
fontsize(15,'points');
legend([h1, h2], 'Location', 'NorthWest', 'FontSize',12); % Shows Line Style key
text(17, 150, {'Color Key:', 'Red: Patient 1', 'Green: Patient 2', 'Blue: Patient 3'},'FontSize',12);
box on; 
ylim([-120 220]);
pbaspect([3 2 1]);
fontname('Consolas');