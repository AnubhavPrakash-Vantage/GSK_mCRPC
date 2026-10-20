clear; close all; clc; warning off;
vantage.style();

%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This path may be changed to ensure raw data is properly loaded
FilePath = "G:\Shared drives\GSK_Oncology\B7_H3_mCRPC\07 Team Working Folders\Suba working folder\Idxd data.xlsx";
LesionSheetname = 'Lesion Waterfall';
PSASheetname = 'PSA Waterfall';
Lesion_filename_to_save = 'LesionWaterfall';
PSA_filename_to_save = 'PSAWaterfall';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Extracting the PD Data and saving........\n");
[LesionWaterfall, PSAWaterfall] = PD_data_viz(FilePath, LesionSheetname, PSASheetname , true); % Analysies and processes the raw data. 
save([Lesion_filename_to_save,'.mat'],"LesionWaterfall");
fprintf(['Saved processed data as ' Lesion_filename_to_save '.mat!\n']);
save([PSA_filename_to_save,'.mat'],"PSAWaterfall");
fprintf(['Saved processed data as ' PSA_filename_to_save '.mat!\n']);
% writetable(Waterfall,[Lesion_filename_to_save,'.xlsx']);
% fprintf(['Saved processed data as ' Lesion_filename_to_save '.xlsx!\n']);
% 




























%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LesionWaterfall, PSAWaterfall] = PD_data_viz(filepath, lesionsheet, psasheet, plot_flag)
    %
    % Extracts, processes, and visualizes clinical PD (tumor size) data from a given file.
    %
    % INPUTS:
    %   filepath:       Location of the Raw tumor growth data.
    %   lesionsheet:
    %   psasheet:
    %   plot_flag:      Boolean value deciding whether to generate plots or just process the data without plotting anything
    %
    % OUTPUTS:
    %   Waterfall:      A table containing the processed best response data.
    %

    LesionWaterfall = readtable(filepath, Sheet=lesionsheet);
    LesionWaterfall.Properties.VariableNames = {'SUBJID','BOR','Dose_mpk'};
    
    PSAWaterfall = readtable(filepath, Sheet=psasheet);
    PSAWaterfall.Properties.VariableNames = {'SUBJID','BOR','Dose_mpk'};

    if plot_flag 
        LesionWaterfall = sortrows(LesionWaterfall,'BOR','descend');
        PSAWaterfall = sortrows(PSAWaterfall,'BOR','descend');
        

        dose_colors = array2table(rgb2hex(vantage.colors())');
        dose_colors.Properties.VariableNames = {'12','8','d3','d2','16','d1','6.4'};

        % Plot Lesion Waterfall     
        figure();
        hold on;
        box on; set(gca,'LineWidth',1.25);
        ylabel('BOR (%)','FontWeight','bold');
        xlabel('Patient No.','FontWeight','bold');
        fontsize(15,'points');
        doses_in_legend =[];
        for patient_id = 1:size(LesionWaterfall,1)
            current_dose = LesionWaterfall.Dose_mpk(patient_id);
            if ismember(current_dose, doses_in_legend)
                h_vis = 'off'; 
            else
                h_vis = 'on';
                doses_in_legend = [doses_in_legend; current_dose]; %#ok<AGROW>
            end
            bar(patient_id, LesionWaterfall.BOR(patient_id),'FaceColor', dose_colors.(num2str(current_dose)),'BarWidth',1,'EdgeColor','k','HandleVisibility',h_vis,'DisplayName',[num2str(current_dose) ' mpk']);
            
        end
        legend('Box','off', 'Location','southwest');
        title('I-Dxd Waterfall plot for \itLesion Data','FontSize',16,'FontWeight','bold','Interpreter','tex');
        fontname('Consolas');
        pbaspect([3 2 1]);
        hold off;   
        

        dose_colors.Properties.VariableNames = {'12','16','d3','d2','6.4','d1','8'};
        % Plot PSA Waterfall     
        figure();
        hold on;
        box on; set(gca,'LineWidth',1.25);
        ylabel('BOR (%)','FontWeight','bold');
        xlabel('Patient No.','FontWeight','bold');
        fontsize(15,'points');
        doses_in_legend =[];
        for patient_id = 1:size(PSAWaterfall,1)
            current_dose = PSAWaterfall.Dose_mpk(patient_id);
            if ismember(current_dose, doses_in_legend)
                h_vis = 'off'; 
            else
                h_vis = 'on';
                doses_in_legend = [doses_in_legend; current_dose]; %#ok<AGROW>
            end
            bar(patient_id, PSAWaterfall.BOR(patient_id),'FaceColor', dose_colors.(num2str(current_dose)),'BarWidth',1,'EdgeColor','k','HandleVisibility',h_vis,'DisplayName',[num2str(current_dose) ' mpk']);
            
        end
        legend('Box','off', 'Location','southwest');
        title('I-Dxd Waterfall plot for \itPSA Data','FontSize',16,'FontWeight','bold','Interpreter','tex');
        fontname('Consolas');
        pbaspect([3 2 1]);
        hold off;   
    end
    drawnow
end
 