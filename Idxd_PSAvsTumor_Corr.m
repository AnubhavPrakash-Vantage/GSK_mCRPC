clc; clear;
%% Setup and Data Loading
filename = 'Idxd data.xlsx';
try
    sheetNames = sheetnames(filename);
catch
    error('Could not read file. Make sure "Idxd data.xlsx" is in your current folder.');
end

% Storage for results
meanCurves = cell(1, 2);
timeGrids = cell(1, 2);
sheetLabels = cell(1, 2);

%% Processing Loop (for both sheets)
for s = 1:2
    % 1. Read Data
    currSheet = sheetNames(s);
    fprintf('\n--- Processing Sheet: %s ---\n', currSheet);
    
    opts = detectImportOptions(filename, 'Sheet', currSheet);
    opts.VariableNamingRule = 'preserve';
    % Force first column to be string to handle 'P1' and empties correctly
    opts = setvartype(opts, 1, 'string'); 
    T = readtable(filename, opts);
    
    % 2. Extract and Clean Data
    % Assumes Col 1 = ID, Col 2 = Time, Col 3 = Value
    rawID = T{:, 1};
    timeData = T{:, 2};
    valData = T{:, 3};
    
    % Fill missing IDs (simulate "Fill Down" for merged cells)
    % If a row has a valid time/value but missing ID, assign previous ID
    for r = 2:length(rawID)
        if ismissing(rawID(r)) || rawID(r) == ""
            rawID(r) = rawID(r-1);
        end
    end
    
    % Remove rows that are still invalid (e.g., NaN time or empty ID)
    validRows = ~ismissing(rawID) & ~isnan(timeData) & ~isnan(valData);
    pID = rawID(validRows);
    time = timeData(validRows);
    val = valData(validRows);
    
    % 3. Identify Unique Patients and Their Max Times
    [uPatients, ~, idx] = unique(pID);
    numPatients = length(uPatients);
    
    maxTimePerPatient = zeros(numPatients, 1);
    for i = 1:numPatients
        patientTime = time(idx == i);
        if ~isempty(patientTime)
            maxTimePerPatient(i) = max(patientTime);
        end
    end
    
    % 4. Define Time Grid
    % "Find patient with longest time series... divide into 20 equidistant points"
    globalMaxTime = max(maxTimePerPatient);
    tQuery = linspace(0, globalMaxTime, 21); % 21 points = 20 intervals
    timeGrids{s} = tQuery;
    
    % 5. Linear Interpolation
    % Rows = Patients, Cols = TimePoints
    interpData = nan(numPatients, length(tQuery));
    
    for i = 1:numPatients
        mask = (idx == i);
        p_t = time(mask);
        p_v = val(mask);
        
        % Sort and Unique
        [p_t, sortIdx] = sort(p_t);
        p_v = p_v(sortIdx);
        [p_t, uniqIdx] = unique(p_t);
        p_v = p_v(uniqIdx);
        
        if length(p_t) < 2
            continue; 
        end
        
        % Interpolate only up to this patient's max time
        p_max = max(p_t);
        validQueryMask = tQuery <= p_max;
        
        if any(validQueryMask)
            interpData(i, validQueryMask) = interp1(p_t, p_v, tQuery(validQueryMask), 'linear');
        end
    end
    
    % --- REQUESTED ADDITION: Print Time Point Arrays ---
    fprintf('Printing data arrays for each time point before averaging:\n');
    for t = 1:length(tQuery)
        % Get column for this time point
        colData = interpData(:, t);
        % Keep only valid numbers (exclude NaNs from patients who dropped out)
        validData = colData(~isnan(colData));
        
        fprintf('  Time Point %2d (t=%.2f): [%s]\n', ...
            t, tQuery(t), strjoin(string(validData'), ', '));
    end
    % ---------------------------------------------------
    
    % 6. Calculate Mean Profile
    meanProfile = mean(interpData, 1, 'omitnan');
    meanCurves{s} = meanProfile;
    sheetLabels{s} = char(currSheet);
end

%% Correlation and Plotting

% Calculate Correlation (shape similarity)
correlationVal = corr(meanCurves{1}', meanCurves{2}', 'Rows', 'complete');
fprintf('\nCorrelation between the two mean curves: %.4f\n', correlationVal);

% Create Plot
fig = figure('Color', 'w');
ax = gca;
hold(ax, 'on');

% Set Colors
try
    colors = orderedcolors('gem');
catch
    warning('orderedcolors not found. Using standard lines.');
    colors = lines(2);
end

% Plot Lines
p1 = plot(timeGrids{1}, meanCurves{1}, '-o', 'Color', colors(1,:), ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'MarkerFaceColor', colors(1,:));
p2 = plot(timeGrids{2}, meanCurves{2}, '-^', 'Color', colors(2,:), ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'MarkerFaceColor', colors(2,:));

% Styling
set(ax, 'FontName', 'Consolas');
set(ax, 'FontSize', 15);
set(ax, 'Box', 'on');
set(ax, 'LineWidth', 1.25); % Plot box line width
grid(ax, 'off');

% Labels
xlabel('Time (Weeks)', 'FontWeight', 'bold');
ylabel('% Change', 'FontWeight', 'bold');
title(sprintf('Mean Profiles (Correlation: %.4f)', correlationVal), 'FontWeight', 'bold');

% Legend
lgd = legend([p1, p2], sheetLabels{1}, sheetLabels{2});
set(lgd, 'FontSize', 12);
set(lgd, 'FontName', 'Consolas');
set(lgd, 'Box', 'off');

hold(ax, 'off');

%% --- APPEND THIS TO THE END OF THE SCRIPT ---

fprintf('\n--- Calculating Correlation on Shared Time Domain ---\n');

% 1. Check which curve is shorter in time (Line A)
maxTime1 = max(timeGrids{1});
maxTime2 = max(timeGrids{2});

if maxTime1 < maxTime2
    % Sheet 1 is the shorter curve (Line A)
    timeA = timeGrids{1};
    valA  = meanCurves{1};
    
    timeB = timeGrids{2};
    valB  = meanCurves{2};
    
    nameA = sheetLabels{1};
    nameB = sheetLabels{2};
else
    % Sheet 2 is the shorter curve (Line A)
    timeA = timeGrids{2};
    valA  = meanCurves{2};
    
    timeB = timeGrids{1};
    valB  = meanCurves{1};
    
    nameA = sheetLabels{2};
    nameB = sheetLabels{1};
end

% 2. Interpolate Line B at Line A's time points
% We use 'linear' interpolation. 'extrap' is generally not needed 
% because Line B covers a larger range (0 to MaxB) than Line A (0 to MaxA).
valB_interpolated = interp1(timeB, valB, timeA, 'linear');

% 3. Calculate Correlation
% Transpose (') to ensure we pass column vectors to the corr function
sharedTimeCorr = corr(valA', valB_interpolated', 'Rows', 'complete');

% 4. Print Results
fprintf('Shorter Curve (A): %s (Duration: %.2f weeks)\n', nameA, max(timeA));
fprintf('Longer Curve  (B): %s (Duration: %.2f weeks)\n', nameB, max(timeB));
fprintf('Correlation (A vs Interpolated B): %.4f\n', sharedTimeCorr);

% 5. (Optional) Update the Plot Title with this new correlation
title(ax, sprintf('Mean Profiles (Shared Time Corr: %.4f)', sharedTimeCorr), 'FontWeight', 'bold');

%% --- APPEND THIS FOR SCATTER PLOT AND LINEAR FIT (Y = mx) ---

fprintf('\n--- Generating Scatter Plot with Linear Fit (y=mx) ---\n');

% 1. Assign X and Y based on Sheet Names
% We need "Lesion" on X and "Sheet1" (or the other dataset) on Y.
% We check the names established in the previous block (nameA, nameB).

if contains(nameA, 'Lesion', 'IgnoreCase', true)
    X_data = valA;
    Y_data = valB_interpolated;
    x_label_str = nameA;
    y_label_str = nameB;
elseif contains(nameB, 'Lesion', 'IgnoreCase', true)
    X_data = valB_interpolated;
    Y_data = valA;
    x_label_str = nameB;
    y_label_str = nameA;
else
    % Fallback if "Lesion" isn't explicitly found in names
    warning('Could not auto-detect "Lesion" sheet. Assigning X to first curve.');
    X_data = valA;
    Y_data = valB_interpolated;
    x_label_str = nameA;
    y_label_str = nameB;
end

% Ensure column vectors
X_data = X_data(:);
Y_data = Y_data(:);

% 2. Calculate Linear Fit forced through origin (y = mx)
% Formula: m = (x' * y) / (x' * x)
slope_m = (X_data' * Y_data) / (X_data' * X_data);
fprintf('Linear Fit Equation: y = %.4f * x\n', slope_m);

% 3. Create Scatter Plot
fig2 = figure('Color', 'w');
ax2 = gca;
hold(ax2, 'on');

% Plot Scatter Points
scatter(X_data, Y_data, 50, 'Filled', ...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'k');

% Plot Fitted Line
% Generate x range for the line (from min X to max X)
x_line = linspace(min([0; X_data]), max([0; X_data]), 100);
y_line = slope_m * x_line;

plot(x_line, y_line, 'r--', 'LineWidth', 2);

% 4. Styling
set(ax2, 'FontName', 'Consolas');
set(ax2, 'FontSize', 15);
set(ax2, 'Box', 'on');
set(ax2, 'LineWidth', 1.25);
grid(ax2, 'off');

% Labels
xlabel(['% Change (' x_label_str ')'], 'FontWeight', 'bold');
ylabel(['% Change (' y_label_str ')'], 'FontWeight', 'bold');
title(sprintf('Linear Fit: y = %.4fx', slope_m), 'FontWeight', 'bold');

% Legend
lgd2 = legend('Data Points', 'Fit (y=mx)', 'Location', 'southeast');
set(lgd2, 'FontSize', 12);
set(lgd2, 'FontName', 'Consolas');
set(lgd2, 'Box', 'off');

hold(ax2, 'off');