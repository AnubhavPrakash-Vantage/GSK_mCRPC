
sz = [1000000, 1]; % Generate random numbers
r = 140*exp(randn(sz)*0.475);
r=r(r<=300);
% Calculate percentages in specific bins
num_samples = length(r);
count_0_100 = sum(r >= 0 & r < 100);
count_100_200 = sum(r >= 100 & r < 200);
count_200_300 = sum(r >= 200 & r <= 300); 

pct_0_100 = (count_0_100 / num_samples) * 100;
pct_100_200 = (count_100_200 / num_samples) * 100;
pct_200_300 = (count_200_300 / num_samples) * 100;


% Plot
figure;
gem_hex_codes = rgb2hex(orderedcolors("gem"));
histogram(r, 100, 'Normalization', 'pdf','EdgeColor',gem_hex_codes(1),'FaceColor',gem_hex_codes(1),'FaceAlpha',1);
xlabel('H-Score', 'FontWeight','bold');
ylabel('Probability Density', 'FontWeight','bold');
box on; set(gca,'LineWidth',1.25);
hold on;
xline(100, 'k--', 'LineWidth', 2, 'Label', '100');
xline(200, 'k--', 'LineWidth', 2, 'Label', '200');
fontsize(15,'points');

str_line1 = sprintf('%% in [0, 100)   = %8.4f%%', pct_0_100);
str_line2 = sprintf('%% in [100, 200) = %8.4f%%', pct_100_200);
str_line3 = sprintf('%% in [200, 300] = %8.4f%%', pct_200_300);
str = {str_line1; str_line2; str_line3};
text(110, 0.006, str,'FontSize', 15, 'Color', gem_hex_codes(7));
title(sprintf('H-Score distribution for Colon Cancer (n = %.f)', num_samples),'FontSize', 16, 'FontWeight','bold');
fontname('Consolas');
pbaspect([3 2 1]);
hold off;


% Display the results
fprintf('\n--- Verification Bin Statistics ---\n');
fprintf('Total Samples: %d\n', num_samples);
fprintf('--------------------------------------\n');
fprintf('Target [0, 100)  :    25.0000%%\n');
fprintf('Obtained [0, 100):   %8.4f%%\n', pct_0_100);
fprintf('--------------------------------------\n');
fprintf('Target [100, 200)  :  57.0000%%\n');
fprintf('Obtained [100, 200): %8.4f%%\n', pct_100_200);
fprintf('--------------------------------------\n');
fprintf('Target [200, 300]  :   18.0000%%\n');
fprintf('Obtained [200, 300]: %8.4f%%\n', pct_200_300);
fprintf('--------------------------------------\n');
fprintf(' Median H-Score = %.2f\n', median(r));
fprintf(' Mean H-Score = %.2f\n', mean(r));
fprintf('--------------------------------------\n');

%%

figure;
swarmchart(ones(1,num_samples),round(r),'blue','filled','o','LineWidth',2,'XJitterWidth',0.75,'MarkerEdgeColor','none','MarkerFaceColor',gem_hex_codes(1),'MarkerFaceAlpha',0.5);
xlabel(' ');
ylabel('H-Score', 'FontWeight','bold');
box on; set(gca,'LineWidth',1.25);
ylim([-10,310]); xlim([0, 2]);
xticks(1); xticklabels('Colon Cancer');
fontsize(15,'points');
title(sprintf('H-Score distribution for Colon Cancer (n = %.f)', num_samples),'FontSize', 16, 'FontWeight','bold');
fontname('Consolas');
pbaspect([1 1 1]);

figure;

% Create the boxchart using the same sampled data
b = boxchart(round(r));
b.BoxFaceColor = gem_hex_codes(1);   
b.BoxEdgeColor = gem_hex_codes(1);
b.LineWidth = 2;
b.BoxWidth =0.1;
b.WhiskerLineColor = gem_hex_codes(1);
b.MarkerStyle = '.'; 
b.MarkerSize = 15;
b.MarkerColor = gem_hex_codes(1);    

xlabel(' ');
ylabel('H-Score', 'FontWeight','bold');
box on; 
set(gca,'LineWidth',1.25);
ylim([-10,310]);
xticklabels('Colon Cancer');
fontsize(15,'points');
title(sprintf('H-Score distribution for Colon Cancer (n = %.f)', num_samples), ...
      'FontSize', 16, 'FontWeight','bold');
fontname('Consolas');
pbaspect([1 1 1]);




function F = beta_error(p)
    alpha = p(1);
    beta = p(2);

    % Constraint the solver to keep alpha and beta in their allowed range
    if alpha <= 0 || beta <= 0
        F = 1e10; % Return a very large error
        return;
    end

    % Define our targets
    target_1_pct = 0.76; % P(X < 1/3)
    target_2_pct = 0.99; % P(X < 2/3)

    % Calculate the error for each target 
    err1 = betacdf(1/3, alpha, beta) - target_1_pct;
    err2 = betacdf(2/3, alpha, beta) - target_2_pct;

    % Return the total squared error
    F = err1^2 + err2^2;
end