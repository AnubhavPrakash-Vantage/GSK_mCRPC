videoName = 'TumorSpider.mp4';
v = VideoWriter(videoName, 'MPEG-4');
v.FrameRate = 3;
open(v);

FigIndices = 3:4:120;

for i = 1:length(FigIndices)
    figNum = FigIndices(i);
    
    if ishandle(figNum)
        figure(figNum);
        
        ax = findobj(gcf, 'type', 'axes');
        if ~isempty(ax)
            title(ax(1), sprintf('Tumor Trajectories | Dose = %d mpk', i), 'FontSize', 18); 
        end
        drawnow; 
        
        frame = getframe(gcf);
        writeVideo(v, frame);
    else
        warning('Figure %d not found. Skipping...', figNum);
    end
end

close(v);