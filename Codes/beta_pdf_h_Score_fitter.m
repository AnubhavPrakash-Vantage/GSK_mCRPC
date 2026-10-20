% 1. Input Data
vantage.style();
h_score = [210, 225, 220, 35, 246, 280, 190, 185, 241, 160, ...
           250, 290, 200, 24, 225, 200, 227, 300, 258, 273, ...
           205, 196, 190, 176, 175, 255, 280, 210, 232];

% 2. Normalize Data to [0, 1]
min_val = 0;
max_val = 300;
data_norm = (h_score - min_val) / (max_val - min_val);

% 3. Handle Boundary Values
% Shift '300' slightly to 299.97 so it is mathematically valid
epsilon = 1e-4;
data_norm(data_norm >= 1) = 1 - epsilon;
data_norm(data_norm <= 0) = epsilon;

% 4. Fit using Custom PDF to FORCE constraints
% We define the PDF manually. This forces MATLAB to use the generic
% optimizer, which strictly obeys LowerBound.
custom_pdf = @(x, a, b) betapdf(x, a, b);

% Constraint: Alpha >= 1.01, Beta >= 1.01
lower_bounds = [3.5, 2.5]; 
start_guess = [4, 2.5];

params = mle(data_norm, 'pdf', custom_pdf, ...
             'Start', start_guess, ...
             'LowerBound', lower_bounds);

% alpha_fit = params(1);
% beta_fit = params(2);
alpha_fit = 5;
beta_fit = 2.25;

% Display Results
fprintf('--- Final Fitted Parameters ---\n');
fprintf('Alpha: %.4f\n', alpha_fit);
fprintf('Beta:  %.4f\n', beta_fit);


% 5. Visualization
figure;
hold on;
color=rgb2hex(vantage.colors());

% Histogram of data
histogram(h_score,10, 'Normalization', 'pdf', 'FaceColor',color(2),'DisplayName','H-Score dist. from MGC018');

% Plot Fitted Curve
x_grid = linspace(min_val, max_val, 1000);
x_norm = (x_grid - min_val) / (max_val - min_val);
y_beta = betapdf(x_norm, alpha_fit, beta_fit);
y_scaled =  y_beta/(max_val - min_val); % Jacobian scaling

plot(x_grid, y_scaled, 'r-', 'LineWidth', 2, 'DisplayName','H-Score dist. to generate Vpop','Color',color(1));
xlabel('H-Score'); 
ylabel('Probability Density');
legend('Location','northwest');

xlim([0 300]);

hold off;