clear; clc; warning off;% close all;


%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This path may be changed to ensure raw data is properly loaded
FilePath = "G:\Shared drives\GSK_Oncology\B7_H3_mCRPC\02 Data\Sep_2026_DataCut\nm-gsk5764227-223054-tumor-v6-mcrpc.csv";
cols_to_filter = {'SUBJIDN','COHORTC', 'NFRLT','AFRLT', 'DGRPN', ...
                'WTBL','SMDIAMBL','PSABL','DVID', 'DV','PCHBL','NTARGET', 'AMT', ...
                'ADY', 'AFSRLT','DVIDN','LNDV', 'CHBL',  'OCC', 'DOSEA', 'OVLRSPN'};
Lesion_filename_to_save = 'New Analysis-sept\Lesion_Waterfall_common';
PSA_filename_to_save = 'New Analysis-sept\PSA_Waterfall_dec_lesion';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Extracting the PD Data ...\n");
[Tumor_Data, PSA_Data, Lesion_Waterfall, PSA_Waterfall] = process_and_plot_data(FilePath, cols_to_filter,  false); % Analysies and processes the raw data. 

% raw=readtable('G:Shared drives\GSK_Oncology\B7_H3_Pilot_Project\06 Slides\Shared_by_client_(No Downloading)\CRPC_HNSCC_11thMay\05312026\nm-gsk5764227-223054-slee-v3-hnscc-mcrpc.csv');
% t_variables = {'SUBJIDN', 'PFSD', 'PFSCNSRN', 'BCRINV', 'RSPPSA50', 'DGRPN', 'RSPPCWG3'};
% summary_data = raw(ismember(raw.SUBJIDN,unique(PSA_Data.SUBJIDN)),t_variables);
% plot_patient_profiles(PSA_Data, Tumor_Data, summary_data)

Lesion_Waterfall = calculate_TTR(Tumor_Data, Lesion_Waterfall);
Mean_spider_clinical = calculate_median_resp_at_clinical_visits(Tumor_Data, 'LESION');
Median_PSA_spider_clinical = calculate_median_resp_at_clinical_visits(PSA_Data, 'PSA');

save([Lesion_filename_to_save,'.mat'],"Lesion_Waterfall");
% fprintf(['Saved processed data as ' Lesion_filename_to_save '.mat!\n']);
save([PSA_filename_to_save,'.mat'],"PSA_Waterfall");
% fprintf(['Saved processed data as ' PSA_filename_to_save '.mat!\n']);
save(['New Analysis-sept\Mean_spider_clinical','.mat'],"Mean_spider_clinical");
save(['New Analysis-sept\Median_PSA_spider_clinical','.mat'],"Median_PSA_spider_clinical");


%% Extra analysis and plotting functions
% plot_dose_timings(FilePath, cols_to_filter);
% plot_dose_modifications(FilePath, cols_to_filter);
% psa_sld_correlation(FilePath);
% psa_ntarget_correlation(FilePath);
% psa_sld_change_auc_corr(Tumor_Data, PSA_Data);

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tumor_Data, PSA_Data, Lesion_Waterfall, PSA_Waterfall] = process_and_plot_data(filepath, cols_to_filter, plot_flag)
    %
    % Extracts, processes, and visualizes clinical PD (tumor size) data from a given file.
    %
    % INPUTS:
    %   filepath:       Location of the Raw tumor growth data.
    %   plot_flag:      Boolean value deciding whether to generate plots or just process the data without plotting anything
    %
    % OUTPUTS:
    %   Waterfall:      A table containing the processed best response data.
    %

    Raw_Data = readtable(filepath);
    Raw_Data = Raw_Data(:,cols_to_filter);
    % Raw_Data = Raw_Data(Raw_Data.SUBJIDN == 161,:);
    % Raw_Data = Raw_Data(Raw_Data.DGRPN == 8 , :);
    Tumor_Data = Raw_Data(strcmp(Raw_Data.DVID,'SUM_OF_DIAMETERS_-_INVESTIGATOR'),:);
    Tumor_Data.PCHBL(Tumor_Data.PCHBL == -999) = 0;
    PSA_Data = Raw_Data(strcmp(Raw_Data.DVID,'PROSTATE_SPECIFIC_ANTIGEN'),:);
    PSA_Data.PCHBL(PSA_Data.PCHBL == -999) = 0;
    Raw_Data.Dose_mpk = round(Raw_Data.DOSEA ./ Raw_Data.WTBL);
    
    

    fields_to_keep = {'SUBJIDN', 'AFRLT', 'PCHBL','DGRPN', 'OVLRSPN'};
    Tumor_Data = sortrows(Tumor_Data,'SUBJIDN');
    PSA_Data = sortrows(PSA_Data,'SUBJIDN');
    rows_to_delete = (PSA_Data.AFRLT(1:end-1) == 0) & (PSA_Data.AFRLT(2:end) == 0);
    PSA_Data(rows_to_delete, :) = [];
    

    Tumor_Pts = unique(Tumor_Data.SUBJIDN);
    cnt = 0; 
    Lesion_Waterfall = table(); % Initialize Waterfall as a table
    
    for patient_id = 1:size(Tumor_Pts,1)
        Temp_tumor = Tumor_Data(Tumor_Data.SUBJIDN == Tumor_Pts(patient_id,1), fields_to_keep);
        
        if max(Temp_tumor.AFRLT) > 0
            cnt = cnt + 1;
            Temp_tumor = Temp_tumor(2:end,:);
            Temp1 = Temp_tumor(Temp_tumor.PCHBL == min(Temp_tumor.PCHBL), :);
            Lesion_Waterfall(cnt,:) = Temp1(1,:);
        end
    end
    
    PSA_Pts = unique(PSA_Data.SUBJIDN);
    cnt = 0; % re-initialise for PSA
    PSA_Waterfall = table();
    for patient_id = 1:size(PSA_Pts,1)
        Temp_PSA = PSA_Data(PSA_Data.SUBJIDN == PSA_Pts(patient_id,1), fields_to_keep);
        if max(Temp_PSA.AFRLT) > 0
            if sum(table2array(Temp_PSA(1:2,'AFRLT'))) ~= 0 
                % This is to check if the first time point correspond to 
                % before Dose is administered, we need to ingore that 
                % timepoint as it is pre-treatment
                cnt = cnt + 1;
                Temp_PSA = Temp_PSA(2:end,:);
                Temp1 = Temp_PSA(Temp_PSA.PCHBL == min(Temp_PSA.PCHBL),:);
                PSA_Waterfall(cnt,:) = Temp1(1,:);
            else
                Temp_PSA = Temp_PSA(3:end,:); % Ignoring the first pre-treatment time points
                if ~isempty(Temp_PSA)
                    cnt = cnt + 1;
                    Temp1 = Temp_PSA(Temp_PSA.PCHBL == min(Temp_PSA.PCHBL),:);
                    PSA_Waterfall(cnt,:) = Temp1(1,:);
                end
            end
        end
    end

    % Renaming the columns for better understanding
    PSA_Waterfall.Properties.VariableNames = {'SUBJID','Actual_time', 'BOR','Dose_mpk', 'OVLRSPN'};
    Lesion_Waterfall.Properties.VariableNames = {'SUBJID','Actual_time', 'BOR','Dose_mpk', 'OVLRSPN'};
    PSA_Waterfall = sortrows(PSA_Waterfall,'BOR','descend'); 
    Lesion_Waterfall = sortrows(Lesion_Waterfall,'BOR','descend');

    PSA_Waterfall.BOR = min(100, PSA_Waterfall.BOR);
    common_subjids = intersect(PSA_Waterfall.SUBJID, Lesion_Waterfall.SUBJID);
    PSA_Waterfall = PSA_Waterfall(ismember(PSA_Waterfall.SUBJID, common_subjids), :);
    Lesion_Waterfall = Lesion_Waterfall(ismember(Lesion_Waterfall.SUBJID, common_subjids), :);
    
    % Selecting only patients with decreasing tumor for PSA Waterfall
    decreasing_tumor_patients = unique(Lesion_Waterfall.SUBJID(Lesion_Waterfall.BOR<0));
    PSA_Waterfall = PSA_Waterfall(ismember(PSA_Waterfall.SUBJID, decreasing_tumor_patients), :);
    PSA_Data = PSA_Data(ismember(PSA_Data.SUBJIDN, decreasing_tumor_patients), :);
    PSA_Data = PSA_Data(ismember(PSA_Data.SUBJIDN, common_subjids), :);
    
    Tumor_Data = Tumor_Data(ismember(Tumor_Data.SUBJIDN, common_subjids),:);
    
    PSA_Pts = unique(PSA_Data.SUBJIDN);
    Tumor_Pts = unique(Tumor_Data.SUBJIDN);

    if plot_flag 
        dose_colors = array2table(rgb2hex(orderedcolors('gem'))');
        dose_colors.Properties.VariableNames = {'8','d3','6','d2','d4','d1','d5'};

        % Plot Lesion Waterfall - all doses combined 
        figure();
        hold on;
        box on; set(gca,'LineWidth',1.25);
        ylabel('BOR (%)','FontWeight','bold');
        xlabel('Patient No.','FontWeight','bold');
        fontsize(15,'points');
        doses_in_legend =[];
        for patient_id = 1:size(Lesion_Waterfall,1)
            current_dose = Lesion_Waterfall.Dose_mpk(patient_id);
            if ismember(current_dose, doses_in_legend)
                h_vis = 'off'; 
            else
                h_vis = 'on';
                doses_in_legend = [doses_in_legend; current_dose]; %#ok<AGROW>
            end
            bar(patient_id, Lesion_Waterfall.BOR(patient_id),'FaceColor', dose_colors.(num2str(current_dose)),'BarWidth',1,'EdgeColor','k','HandleVisibility',h_vis,'DisplayName',[num2str(current_dose) ' mpk']);
            
        end
        legend('Box','off', 'Location','southwest');
        title('GSK-227 Waterfall plot for \itLesion Data','FontSize',16,'FontWeight','bold','Interpreter','tex');
        fontname('Poppins');
        pbaspect([3 2 1]);
        ylim([-110, 40]); xlim([0,47])
        hold off;   
        
        % Plot Lesion Waterfall - all doses separately
        doses = unique(Lesion_Waterfall.Dose_mpk);
        for dose_idx = 1:length(doses)
            temp_waterfall = Lesion_Waterfall(Lesion_Waterfall.Dose_mpk == doses(dose_idx),:);
            figure();
            hold on;
            box on; set(gca,'LineWidth',1.25);
            ylabel('BOR (%)','FontWeight','bold');
            xlabel('Patient No.','FontWeight','bold');
            fontsize(15,'points');
            color_dose = dose_colors.(num2str(doses(dose_idx)));
            for patient_id = 1:size(temp_waterfall,1)
                bar(patient_id, temp_waterfall.BOR(patient_id),'FaceColor', color_dose,'BarWidth',1,'EdgeColor','k');
            end
            title(sprintf('GSK-227 Waterfall plot for Lesion Data - %d mpk',doses(dose_idx)),'FontSize',16,'FontWeight','bold');
            fontname('Poppins');
            pbaspect([3 2 1]);
            
            ylim([-110, 40]);
            hold off;   
        end    


        % Plotting Lesion Spider plot - all dose combined
        figure;
        hold on;
        for i = 1:length(Tumor_Pts)
        current_idx = (Tumor_Data.SUBJIDN == Tumor_Pts(i));
        dose = Tumor_Data.DGRPN(current_idx);
        color_dose = dose_colors.(num2str(dose(1,1)));
        current_time = Tumor_Data.AFRLT(current_idx);
        subject_dose = Tumor_Data.PCHBL(current_idx);
        plot(current_time, subject_dose, '-', 'Color', color_dose, 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'HandleVisibility','off');
        end
        plot(NaN, NaN, '-', 'Color', dose_colors.('6'), 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'DisplayName','6 mpk');
        plot(NaN, NaN, '-', 'Color', dose_colors.('8'), 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'DisplayName','8 mpk');
        legend('Location','southeast');
        hold off;
        xlabel('Time (days) - AFRLT');
        ylabel('SLD Change (%)');
        title('Tumor Trajectories | mCRPC Cohort - 6 and 8 mpk combined');

        % Plotting Lesion Spider plot - all dose separately
        doses = unique(Tumor_Data.DGRPN);
        for dose_idx = 1:length(doses)
            temp = Tumor_Data(Tumor_Data.DGRPN==doses(dose_idx),:);
            temp_Pts = unique(temp.SUBJIDN);
            figure;
            hold on;
            for i = 1:length(temp_Pts)
            current_idx = (temp.SUBJIDN == temp_Pts(i));
            current_time = temp.AFRLT(current_idx);
            subject_dose = temp.PCHBL(current_idx);
            plot(current_time, subject_dose, '-', 'Color', dose_colors.(num2str(doses(dose_idx))), 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'HandleVisibility','off');
            end
            
            hold off;
            xlabel('Time (days) - AFRLT');
            ylabel('SLD Change (%)');
            title(sprintf('Tumor Trajectories | mCRPC Cohort - %d mpk combined',doses(dose_idx)));
        end


        % Plot PSA Waterfall - all doses combined 
        figure();
        hold on;
        box on; set(gca,'LineWidth',1.25);
        ylabel('BOR (%)','FontWeight','bold');
        xlabel('Patient No.','FontWeight','bold');
        fontsize(15,'points');
        doses_in_legend =[];
        for patient_id = 1:size(PSA_Waterfall,1)
            current_dose = PSA_Waterfall.Dose_mpk(patient_id);
            if ismember(current_dose, doses_in_legend)
                h_vis = 'off'; 
            else
                h_vis = 'on';
                doses_in_legend = [doses_in_legend; current_dose]; %#ok<AGROW>
            end
            bar(patient_id, PSA_Waterfall.BOR(patient_id),'FaceColor', dose_colors.(num2str(current_dose)),'BarWidth',1,'EdgeColor','k','HandleVisibility',h_vis,'DisplayName',[num2str(current_dose) ' mpk']);
            
        end
        legend('Box','off', 'Location','southwest');
        title('GSK-227 Waterfall plot for \itPSA Data','FontSize',16,'FontWeight','bold','Interpreter','tex');
        fontname('Poppins');
        pbaspect([3 2 1]);
        hold off;   
        
        % Plot PSA Waterfall - all doses separately
        doses = unique(PSA_Waterfall.Dose_mpk);
        for dose_idx = 1:length(doses)
            temp_waterfall = PSA_Waterfall(PSA_Waterfall.Dose_mpk == doses(dose_idx),:);
            figure();
            hold on;
            box on; set(gca,'LineWidth',1.25);
            ylabel('BOR (%)','FontWeight','bold');
            xlabel('Patient No.','FontWeight','bold');
            fontsize(15,'points');
            color_dose = dose_colors.(num2str(doses(dose_idx)));
            for patient_id = 1:size(temp_waterfall,1)
                bar(patient_id, min(temp_waterfall.BOR(patient_id),100),'FaceColor', color_dose,'BarWidth',1,'EdgeColor','k');
            end
            title(sprintf('GSK-227 Waterfall plot for PSA Data - %d mpk',doses(dose_idx)),'FontSize',16,'FontWeight','bold');
            fontname('Poppins');
            pbaspect([3 2 1]);
            ylim([-105, 105]);
            hold off;   
        end    
        
        % Plotting PSA Spider plot - all dose combined
        figure;
        hold on;
        for i = 1:length(PSA_Pts)
        current_idx = (PSA_Data.SUBJIDN == PSA_Pts(i));
        dose = PSA_Data.DGRPN(current_idx);
        color_dose = dose_colors.(num2str(dose(1,1)));
        current_time = PSA_Data.AFRLT(current_idx);
        subject_dose = PSA_Data.PCHBL(current_idx);
        plot(current_time, subject_dose, '-', 'Color', color_dose, 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'HandleVisibility','off');
        end
        plot(NaN, NaN, '-', 'Color', dose_colors.('6'), 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'DisplayName','6 mpk');
        plot(NaN, NaN, '-', 'Color', dose_colors.('8'), 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'DisplayName','8 mpk');
        legend('Location','southeast');
        hold off;
        xlabel('Time (days) - AFRLT');
        ylabel('PSA Change (%)');
        title('PSA Trajectories | mCRPC Cohort - 6 and 8 mpk combined');

        % Plotting PSA Spider plot - all dose separately
        doses = unique(PSA_Data.DGRPN);
        for dose_idx = 1:length(doses)
            temp = PSA_Data(PSA_Data.DGRPN==doses(dose_idx),:);
            temp_Pts = unique(temp.SUBJIDN);
            figure;
            hold on;
            for i = 1:length(temp_Pts)
            current_idx = (temp.SUBJIDN == temp_Pts(i));
            current_time = temp.AFRLT(current_idx);
            subject_dose = temp.PCHBL(current_idx);
            plot(current_time, subject_dose, '-', 'Color', dose_colors.(num2str(doses(dose_idx))), 'LineWidth', 2, 'Marker', '.','MarkerSize',15,'HandleVisibility','off');
            end
            
            hold off;
            xlabel('Time (days) - AFRLT');
            ylabel('PSA Change (%)');
            title(sprintf('PSA Trajectories |  mCRPC Cohort - %d mpk combined',doses(dose_idx)));
        end
        drawnow
    end
end
 
function Lesion_Waterfall = calculate_TTR(Tumor_Data, Lesion_Waterfall)
    
    num_subjects = height(Lesion_Waterfall);
    Lesion_Waterfall.TTR = NaN(num_subjects, 1);
    
    for i = 1:num_subjects
        subject_id = Lesion_Waterfall.SUBJID(i);
        
        subj_data = Tumor_Data(Tumor_Data.SUBJIDN == subject_id, :);
                
        t = subj_data.AFRLT; y = subj_data.PCHBL;
        
        idx = find(y <= -30, 1, 'first');
        
        if ~isempty(idx)
            

            % Get the coordinates for the points immediately before and after crossing
            t0 = t(idx-1);
            y0 = y(idx-1);
            t1 = t(idx);
            y1 = y(idx);
            
            % Safeguard against division by zero in case of identical y values
            if y1 == y0
                Lesion_Waterfall.TTR(i) = t1;
            else
                % Linear interpolation 
                Lesion_Waterfall.TTR(i) = t0 + ((-30 - y0) / (y1 - y0)) * (t1 - t0);
            end
            
        end
        % If idx is empty (trajectory never crossed -30), the TTR value remains NaN
    end
end
 
function Mean_response_table = calculate_median_resp_at_clinical_visits(Data, mode)
    
    
    % Data.PCHBL = min(100, Data.PCHBL);
    % creating a Visit_num column for easier filtering at later point to get mean spider 
    num_rows = height(Data);
    Visit_Num = zeros(num_rows, 1);
    current_subj = Data.SUBJIDN(1);
    visit_count = 0; 
    
    for i = 1:num_rows
        if Data.SUBJIDN(i) == current_subj
            Visit_Num(i) = visit_count;
            visit_count = visit_count + 1;
        else
            % New patient found, reset tracker
            current_subj = Data.SUBJIDN(i);
            visit_count = 0;
            Visit_Num(i) = visit_count;
            visit_count = visit_count + 1;
        end
    end
    Data.Visit_Num = Visit_Num;
    
    unique_visits = unique(Data.Visit_Num);
    mean_resp_list = NaN(length(unique_visits), 1);
    min_resp_list = NaN(length(unique_visits), 1);
    max_resp_list = NaN(length(unique_visits), 1);
    for k = 1:length(unique_visits)
        vis = unique_visits(k);
        responses = Data.PCHBL(Data.Visit_Num == vis);
        
        if strcmp(mode, 'PSA')
            mean_resp_list(k) = median(responses, 'omitnan');
        else
            mean_resp_list(k) = mean(responses, 'omitnan');
        end
        min_resp_list(k) = min(responses);
        max_resp_list(k) = max(responses);
    end
    
    Mean_response_table = table(unique_visits, mean_resp_list, min_resp_list, max_resp_list, 'VariableNames', {'Clinical_visits', 'Mean_Resp_Clinical', 'Min_Resp_Clinical', 'Max_Resp_Clinical'});
end
 
function plot_dose_timings(filepath, cols_to_filter)
    Raw_Data = readtable(filepath);
    Raw_Data = Raw_Data(:, cols_to_filter);
    Raw_Data = Raw_Data(Raw_Data.DGRPN == 6 , :);
    Dose_Data = Raw_Data(strcmp(Raw_Data.DVID,'ACTUAL_DOSE'), :);
    
    unique_patients = unique(Dose_Data.SUBJIDN);
    N = length(unique_patients);
    
    [~, y_idx] = ismember(Dose_Data.SUBJIDN, unique_patients);
    
    figure;
    hold on;     
    scatter(Dose_Data.AFRLT, y_idx, 20, 'k', 'filled');
    max_time = max(Dose_Data.AFRLT);
    if ~isempty(max_time) && max_time > 0
        multiples_of_21 = 21:21:max_time; 
        xline(multiples_of_21, '-b', 'Alpha', 0.5); 
        multiples_of_28 = 28:28:max_time; 
        xline(multiples_of_28, '--r', 'Alpha', 0.5);

    end
    fontsize(8,"points");
    xlabel('Actual Time of Dose administration (days)','FontSize',15);
    ylabel('Patient ID (SUBJIDN)','FontSize',15);
    title('Dose Administration Times per Patient - 6 mpk','FontSize',18);
    
    yticks(1:N);
    yticklabels(string(unique_patients)); 
    pbaspect([1 1 1]);
    grid off;
    hold off;
end
 
function plot_dose_modifications(filepath, cols_to_filter)
    Raw_Data = readtable(filepath);
    Raw_Data = Raw_Data(:, cols_to_filter);
    % Raw_Data = Raw_Data(Raw_Data.DGRPN == 8 , :);
    Dose_Data = Raw_Data(strcmp(Raw_Data.DVID,'ACTUAL_DOSE'), :);
    Dose_Data.Dose_mpk  = Dose_Data.DOSEA./Dose_Data.WTBL;
    unique_patients = unique(Dose_Data.SUBJIDN);
    N = length(unique_patients);
    
    figure('WindowState', 'maximized');
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = gobjects(N, 1);
    for pt_idx = 1:N
        pt = unique_patients(pt_idx);
        ax(pt_idx) = nexttile;
        % title(ax(pt_idx), "Patient " + num2str(pt));       
        temp_data = Dose_Data(Dose_Data.SUBJIDN==pt, :);
        plot(temp_data.AFRLT, temp_data.Dose_mpk, Marker=".", MarkerSize=15, LineStyle="-", DisplayName=num2str(pt));
        ax(pt_idx).Title.String = "Patient " + num2str(pt);
        ax(pt_idx).Title.Visible = 'on';
        fontsize(ax(pt_idx), 9, "points");
        fontname(ax(pt_idx), 'Poppins');
        grid(ax(pt_idx), 'off');
    end

    
    xlabel(t,'Actual Time of Dose administration (days)','FontSize',12,'FontName','Poppins');
    ylabel(t,'Exact Dose administered per kg (mg/kg)','FontSize',12,'FontName','Poppins');
    title(t,'Dose Modifications per Patient','FontSize',18,'FontName','Poppins','FontWeight', 'bold');
    linkaxes(ax, 'xy');
    % ylim([5.2,6.4]);
    % lgd = legend('Location','bestoutside','Box','on','NumColumns',2);
    % title(lgd,"Patients' Subject ID");
    % pbaspect([1 1 1]);
end
 
function plot_patient_profiles(PSA_Data, Tumor_Data, summary_data)
    % PLOT_PATIENT_PROFILES Generates longitudinal % change plots for SLD and PSA
    % across 24 patients split into two 3x4 subplot figures.
    
    colors_matrix = orderedcolors('gem');
    hex_colors = rgb2hex(colors_matrix)';
    dose_colors = array2table(hex_colors);
    dose_colors.Properties.VariableNames = {'8','6','d3','d2','d4','d1','d5'};
    
    unique_patients = unique(summary_data.SUBJIDN, 'stable');
    num_patients = length(unique_patients);
    
    for idx = 1:num_patients
        pt_id = unique_patients(idx);
        
        % Determine figure assignment and local subplot grid position
        fig_num = ceil(idx / 12);
        sub_pos = mod(idx - 1, 12) + 1;
        
        if sub_pos == 1
            figure('Name', sprintf('Patient Profiles - Part %d', fig_num), ...
                   'Units', 'normalized');
        else
            figure(findobj('Name', sprintf('Patient Profiles - Part %d', fig_num)));
        end
        
        ax = subplot(3, 4, sub_pos);
        hold(ax, 'on');
        pt_psa = PSA_Data(PSA_Data.SUBJIDN == pt_id, :);
        [~, sort_idx_psa] = sort(pt_psa.AFRLT);
        pt_psa = pt_psa(sort_idx_psa, :);
  
        pt_tumor = Tumor_Data(Tumor_Data.SUBJIDN == pt_id, :);
        [~, sort_idx_tum] = sort(pt_tumor.AFRLT);
        pt_tumor = pt_tumor(sort_idx_tum, :);
        
        pt_summary = summary_data(summary_data.SUBJIDN == pt_id, :);
        

        dose_group = pt_summary.DGRPN(1);
        if dose_group == 8
            color_hex = dose_colors.('8'){1};
        else 
            color_hex = dose_colors.('6'){1};
        end
        
        
        plot(ax, pt_tumor.AFRLT, pt_tumor.PCHBL, '-', 'Color', color_hex, 'LineWidth', 2);
        plot(ax, pt_psa.AFRLT, pt_psa.PCHBL, ':', 'Color', color_hex, 'LineWidth', 2);

        if pt_summary.RSPPCWG3 == 1
            resp_str = 'Responder';
        else
            resp_str = 'Non-Responder';
        end
        
        if pt_summary.PFSCNSRN == 0
            cens_str = 'Dropped';
        else
            cens_str = 'Censored';
        end
        
        title_text = sprintf('Subj. %d | %s | %s', pt_id, resp_str, cens_str);
        grid(ax, 'off');
        box(ax, 'on');
        title(ax, title_text, 'FontName', 'Poppins', 'FontSize', 10, 'FontWeight', 'bold');
        xlabel(ax, 'Time (AFRLT)', 'FontName', 'Poppins', 'FontSize', 9);
        ylabel(ax, '% Change from Baseline', 'FontName', 'Poppins', 'FontSize', 9);
        set(ax, 'FontName', 'Poppins', 'FontSize', 9);
        hold(ax, 'off');
    end
end

function psa_sld_correlation(filepath)
    Raw_Data = readtable(filepath);
    Raw_Data = Raw_Data(:,{'SUBJIDN','SMDIAMBL','PSABL'});
    Raw_Data = Raw_Data(Raw_Data.SMDIAMBL ~= -999,:);
    u_subj = unique(Raw_Data.SUBJIDN);
    psa = nan(1,length(u_subj));
    sld = nan(1, length(u_subj));
    for i = 1:length(u_subj)
        psa(i) = max(table2array(Raw_Data(Raw_Data.SUBJIDN == u_subj(i),'PSABL')));
        sld(i) = max(table2array(Raw_Data(Raw_Data.SUBJIDN == u_subj(i),'SMDIAMBL')));
    end
    mdl = fitlm(sld, psa);
    coeffs = mdl.Coefficients.Estimate;
    ci = coefCI(mdl);
    
    x_line = linspace(min(sld), max(sld), 100);
    y_line = coeffs(1) + coeffs(2) * x_line;
    
    figure; hold on; box on; grid off; set(gca, 'LineWidth', 1.25);
    fontname('Poppins'); fontsize(15,"points");
    scatter(sld, psa, 20, 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(x_line, y_line, '-r', 'LineWidth', 1.5);
    
    txt = sprintf('y = %.3fx + %.3f\nSlope 95%% CI: [%.3f, %.3f]\nIntercept 95%% CI: [%.3f, %.3f]', ...
        coeffs(2), coeffs(1), ci(2,1), ci(2,2), ci(1,1), ci(1,2));
    
    text(0.95, 0.95, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % yscale('log');
    ylabel('PSA Baseline serum conc.'); 
    xlabel('SLD Baseline value');
    title('Correlation between PSA and SLD baseline values');
end


function psa_ntarget_correlation(filepath)
    Raw_Data = readtable(filepath);
    Raw_Data = Raw_Data(:,{'SUBJIDN','NTARGET','PSABL'});
    Raw_Data = Raw_Data(Raw_Data.NTARGET ~= -999,:);
    u_subj = unique(Raw_Data.SUBJIDN);
    psa = nan(1,length(u_subj));
    ntarget = nan(1, length(u_subj));
    for i = 1:length(u_subj)
        psa(i) = max(table2array(Raw_Data(Raw_Data.SUBJIDN == u_subj(i),'PSABL')));
        ntarget(i) = max(table2array(Raw_Data(Raw_Data.SUBJIDN == u_subj(i),'NTARGET')));
    end
    mdl = fitlm(ntarget, psa);
    coeffs = mdl.Coefficients.Estimate;
    ci = coefCI(mdl);
    
    x_line = linspace(min(ntarget), max(ntarget), 100);
    y_line = coeffs(1) + coeffs(2) * x_line;
    
    figure; hold on; box on; grid off; set(gca, 'LineWidth', 1.25);
    fontname('Poppins'); fontsize(15,"points");
    scatter(ntarget, psa, 20, 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(x_line, y_line, '-r', 'LineWidth', 1.5);
    
    txt = sprintf('y = %.3fx + %.3f\nSlope 95%% CI: [%.3f, %.3f]\nIntercept 95%% CI: [%.3f, %.3f]', ...
        coeffs(2), coeffs(1), ci(2,1), ci(2,2), ci(1,1), ci(1,2));
    
    text(0.95, 0.95, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % yscale('log');
    ylabel('PSA Baseline serum conc.'); 
    xlabel('No. of Target lesions');
    title('Correlation between PSA baseline values and count of target lesions');
end


function psa_sld_change_auc_corr(Tumor_Data, PSA_Data)
    
    tumor_auc_tbl = get_normalized_auc(Tumor_Data, 'Tumor_AUC');
    psa_auc_tbl = get_normalized_auc(PSA_Data, 'PSA_AUC');
    
    Final_Table = innerjoin(tumor_auc_tbl, psa_auc_tbl, 'Keys', 'SUBJIDN');
    
    sld = Final_Table.Tumor_AUC;
    psa = Final_Table.PSA_AUC;
    
    % Linear Regression
    mdl = fitlm(sld, psa);
    coeffs = mdl.Coefficients.Estimate;
    ci = coefCI(mdl);
    x_line = linspace(min(sld), max(sld), 100);
    y_line = coeffs(1) + coeffs(2) * x_line;
    
    % Plotting
    figure; hold on; box on; grid off; set(gca, 'LineWidth', 1.25);
    fontname('Poppins'); fontsize(15,"points");
    scatter(sld, psa, 40, 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(x_line, y_line, '-r', 'LineWidth', 1.5);
    txt = sprintf('y = %.3fx + %.3f\nSlope 95%% CI: [%.3f, %.3f]\nIntercept 95%% CI: [%.3f, %.3f]', ...
        coeffs(2), coeffs(1), ci(2,1), ci(2,2), ci(1,1), ci(1,2));
    
    text(0.95, 0.95, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    ylabel('Normalized PSA Change AUC'); 
    xlabel('Normalized Tumor (SLD) Change AUC');
    title('Correlation between Tumor and PSA Normalized AUCs');
end


function auc_tbl = get_normalized_auc(data, var_name)
    u_subj = unique(data.SUBJIDN);
    auc_vals = nan(length(u_subj), 1);
    valid = false(length(u_subj), 1); % To keep track of which subjects has valid data
    
    for i = 1:length(u_subj)
        sub_data = data(data.SUBJIDN == u_subj(i), :);
        
        [afrlt_sorted, sort_idx] = sort(sub_data.AFRLT);
        pchbl_sorted = sub_data.PCHBL(sort_idx);
        
        % Handle duplicate 0 time points and Ignore the first one
        zero_idx = find(afrlt_sorted == 0);
        if length(zero_idx) > 1
            afrlt_sorted(zero_idx(1)) = [];
            pchbl_sorted(zero_idx(1)) = [];
        end
        
        if length(afrlt_sorted) > 1
            time_span = abs(afrlt_sorted(end) - afrlt_sorted(1));
            
            if time_span > 0
                raw_auc = trapz(afrlt_sorted, pchbl_sorted);
                
                % Normalize by total time length
                auc_vals(i) = raw_auc / time_span;
                valid(i) = true;
            end
        end
    end
    
    auc_tbl = table(u_subj(valid), auc_vals(valid), 'VariableNames', {'SUBJIDN', var_name});
end