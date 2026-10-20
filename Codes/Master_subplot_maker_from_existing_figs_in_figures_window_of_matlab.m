%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file can take multiple existing figures from the figure windows and 
% plot them combined in a single master plot as subplots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masterFig = figure('Name', 'PSA Spiders', 'NumberTitle', 'off');
tlo = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Change as per the number of figures

% Define figure indices to be collecated in the master figure
Indices_to_use = 4:4:36;
Doses = [1 2 4 6 8 10 12 16 20]; % This is just for adding titles to each subplot
for i = 1:length(Indices_to_use)
    figNum = Indices_to_use(i);
    
    % Check if the figure actually exists to avoid errors
    if ishandle(figNum)
        srcAx = findobj(figure(figNum), 'type', 'axes');
        
        if ~isempty(srcAx)
            % Create destination tile
            destTile = nexttile(tlo);
            
            % copy the children of the source axes into the destination tile
            copyobj(srcAx.Children, destTile);
            
            % Carry over labels
            title(destTile, sprintf('Dose %d mpk', Doses((figNum)/4)));
            xlabel(destTile, srcAx.XLabel.String);
            ylabel(destTile, srcAx.YLabel.String);
            % Match the axis limits and scales (log/linear)
            destTile.XScale = srcAx.XScale;
            destTile.YScale = srcAx.YScale;
            destTile.XLim = srcAx.XLim;

            % ylim(destTile, [-105, 200]);
            destTile.YLim = srcAx.YLim;
        end
    else
        warning('Figure %d does not exist.', figNum);
    end
end
