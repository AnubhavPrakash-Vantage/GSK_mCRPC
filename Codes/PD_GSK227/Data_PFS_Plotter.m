load("Lesion_Waterfall_common.mat");
raw=readtable('G:Shared drives\GSK_Oncology\B7_H3_Pilot_Project\06 Slides\Shared_by_client_(No Downloading)\CRPC_HNSCC_11thMay\05312026\nm-gsk5764227-223054-slee-v3-hnscc-mcrpc.csv');
t_variables = {'SUBJIDN', 'PFSD', 'PFSCNSRN', 'BCRINV', 'RSPPSA50', 'DGRPN', 'RSPPCWG3'};
summary_data = raw(ismember(raw.SUBJIDN,unique(Lesion_Waterfall.SUBJID)),t_variables);


[f, x, flo, fup] = ecdf(summary_data.PFSD, ...
                        'Censoring', summary_data.PFSCNSRN, ...
                        'Function', 'survivor');

% At Day 0 (where survival f == 1), bounds are logically 1
flo(isnan(flo) & f == 1) = 1;
fup(isnan(fup) & f == 1) = 1;

% If survival drops to 0 at the end, bounds are logically 0
flo(isnan(flo) & f == 0) = 0;
fup(isnan(fup) & f == 0) = 0;

% For trailing censored edge cases, carry the last known valid bound forward/backward
flo = fillmissing(flo, 'previous');
flo = fillmissing(flo, 'next');
fup = fillmissing(fup, 'previous');
fup = fillmissing(fup, 'next');
flo = max(0, min(1, flo));
fup = max(0, min(1, fup));

gem_palette = orderedcolors('gem');
line_color  = gem_palette(7, :); 
shade_color = gem_palette(7, :); 

[x_steps, y_upper] = stairs(x, fup);
[~,       y_lower] = stairs(x, flo); 

X_poly = [x_steps; flipud(x_steps)];
Y_poly = 100.*[y_upper; flipud(y_lower)];

figure(); hold on;
hCI = fill(X_poly, Y_poly, shade_color, ...
           'FaceAlpha', 0.25, ...  % Balanced transparency
           'EdgeColor', 'none');   % Removes solid perimeter lines

hPFS = stairs(x, 100*f, '-', 'Color', line_color, 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('PFS Probability (%)', 'FontSize', 13, 'FontWeight', 'bold');
title('Progression-Free Survival Plot', 'FontSize', 15, 'FontWeight', 'bold');

ylim([0 105]);
xlim([0 max(summary_data.PFSD) + 10]);

legend([hPFS, hCI], {'PFS Estimate', '95% Confidence Interval'}, ...
       'Location', 'Best', 'FontSize', 11);

hold off;