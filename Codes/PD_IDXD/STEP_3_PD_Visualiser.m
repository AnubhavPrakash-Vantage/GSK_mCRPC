clear; warning off;
load("LesionWaterfall.mat");                                    % Loads the processed Data
load("PSAWaterfall.mat"); 
load('Sim_Params.mat');                                    % Loads the calibrated and pre-assigned parameters for simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load Project model
proj = sbioloadproject("../../GSK227_CRPC.sbproj");     % This path has to be changed as per the location of the model
model = proj.m1; 


seed = 300:319;                                    % This is the RNG Seed value to which you are simulating
Sim_Params.Dose = LesionWaterfall.Dose_mpk;        % in mg/kg 
Sim_Params.N_VP = height(LesionWaterfall.BOR);     % Number of patients to simulate is same as that in the data
Sim_Params.StopTime = 395;                         % Simulating for 1 year
Sim_Params.Time_Sim = 0:49:Sim_Params.StopTime;    % for each of the days, time point where the simulations will be carried out
Sim_Params.Dose_Cycle_length = 21;                 % in days, 21 Days for Q3W
Sim_Params.Doses_days_in_a_cycle = 0;              % On what days you want to administer the dosage, it could be [0,7] if you want d1,d8 Q3W kind of dosing
Sim_Params.Dose_target_name = 'Central_ADC.Cen_ADC_nM'; % Species name in the model where ADC is dosed
Sim_Params.Dar_Dose_target_name = 'Central_ADC.Cen_DummyDAR'; % Species name in the model where DAR is dosed
Sim_Params.Scan_Interval = 49;                      % Time interval between two radiographical scan in simulation, (in days)
Sim_Params.Time_efficacy=0:Sim_Params.Scan_Interval:Sim_Params.StopTime; % response measument days 
Sim_Params.UseParallel_simfcn = true;

%First param must be the total volume of tumor
Sim_Params.StatesToLog = {'TV','PSA', 'Cen_ADC_ugml', 'Cen_tAb_ugml', 'Cen_PL_ngml'}; 
Sim_Params.N_tmr = 1;
Sim_Params.Param_names_for_intrapatient_heterogeneity = {'Tumor_Dia','H_Score','K1','Kmax_p'}; % Parameters that will change across within each patient
Sim_Params.Intrapatient_Heterogeneity = 0.0;       % Variability in parameters across different tumors within patient
Sim_Params.dropout_Switch = true;                  % Switch to turn off/on the dropouts due to New metastats or Non-specific reasons
Sim_Params.NM_Params = [0.01 0.21 15 1.0];         % New-metastates related dropout probability modelling parameters
Sim_Params.NonSpecific_dropouts = 0.21;            % Probability for dropout related to non-radiographical reasons
plot_data_on_waterfall = true;

model_params = { 
    % all the parameters you want to change from what was used in STEP_2
    %param_name          %mu    %sig
    'W'                  70     0;
    'Switch_TumorComp'   1      0;
    'Switch_Tox'         1      0
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SECTION END %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Params = update_params_from_estimates(Sim_Params, model_params, model);

tumor_waterfall_for_all_seed  = nan(length(seed),Sim_Params.N_VP);
psa_waterfall_for_all_seed  = nan(length(seed),height(PSAWaterfall));
pfs_for_all_seed = cell(length(seed));

for s=1:length(seed)
    tic;
    rng(seed(s));
    fprintf('Simulating for seed number %d \n',s);
    Simulation_Results = run_population_simulation(model, Sim_Params, true);
    
    % Storing results for plotting 
    sim_tumor_waterfall = sortrows(Simulation_Results.AllPlots.Sim_Waterfall,-2);
    tumor_waterfall_for_all_seed(s,:) = sim_tumor_waterfall(:,2)';
    psawaterfall = sim_tumor_waterfall(end-height(PSAWaterfall)+1:end,4);
    psa_waterfall_for_all_seed(s,:) = sort(psawaterfall, 'descend')';
    pfs_for_all_seed{s} = Simulation_Results.AllPlots.Sim_PFS;
    
    % PLOTTING
    % plot_waterfall(Simulation_Results, LesionWaterfall, PSAWaterfall); 
    % plot_pfs(Simulation_Results);                             % [OPTIONAL] Lines marked as optional can be commented without any issue
    % plot_swimmer(Simulation_Results);                         % [OPTIONAL]
    plot_spider(Simulation_Results);                            % [OPTIONAL]
    plot_psa_spider(Simulation_Results);                            % [OPTIONAL]
    % plot_observables(Simulation_Results, Sim_Params)          % [OPTIONAL] 
    % summarize_patient_response(Simulation_Results);           % [OPTIONAL]
    drawnow;
    toc;
end
%%
Plot_waterfall_for_all_seed(psa_waterfall_for_all_seed, tumor_waterfall_for_all_seed, plot_data_on_waterfall, LesionWaterfall, PSAWaterfall);
% Plot_PFS_for_all_seed(Sim_Params, pfs_for_all_seed, true);
drawnow;













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTION DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
function Simulation_Results = run_population_simulation(model, Sim_Params, verbose)
    %
    % run_population_simulation: Simulates a virtual population for a given mode and returns all processed results for plotting.
    %
    % INPUTS:
    %   model:              A SimBiology model object (pre-loaded)
    %   Sim_Params:         Information about the crucial parameters to use while simulation
    %   verbose:            If true, prints information while running
    %
    % OUTPUTS:
    %   Simulation_Results:     A struct containing all relevant simulated outputs:
    %                       .AllPlots (Tumor trajectory, Waterfall, PFS, Dropout related informat)
    %                       .VP_Init_Array (Array containing the initial patient related parameters)
    %                       .Dose (the dose amount in mg/Kg)
    %                       .time_series_data_for_all_pts : Time series data for all the logged states(observables)
    %
    
    
    if verbose
        disp('Running Simulation with estimated parameters...');
    end
    

    % Configure SimBiology Model
    csObj = getconfigset(model);
    csObj.SolverType = 'sundials';
    csObj.SolverOptions.RelativeTolerance = 1e-6;
    csObj.SolverOptions.AbsoluteTolerance = 1e-6;
    csObj.SolverOptions.MaxStep = 1;
    csObj.SolverOptions.AbsoluteToleranceScaling = false;
    csObj.StopTime = Sim_Params.StopTime;
    csObj.TimeUnits = 'day';
    csObj.SolverOptions.OutputTimes = Sim_Params.Time_Sim;
    csObj.RuntimeOptions.StatesToLog = Sim_Params.StatesToLog;    

    IHC_Category = nan(Sim_Params.N_VP,1);
                                                                            
    % Initialising VP_Init_Array
    p_names = fieldnames(Sim_Params.param)';
    all_fields = [{'PatientID', 'TumorID'}, p_names];
    
    % Preallocate the full column vector (Total_Rows x 1)
    % This creates the full structure array filled with empty values
    VP_Init_Array = repmat(cell2struct(cell(length(all_fields), 1), all_fields), Sim_Params.N_VP * Sim_Params.N_tmr, 1);
     
    % Initialize Dose-Level Arrays
    n_timepoints = numel(Sim_Params.Time_efficacy);
    SLD = nan(Sim_Params.N_VP * n_timepoints, 6);
    Sim_Vol = cell(Sim_Params.N_VP, 1); 
    PSA = cell(Sim_Params.N_VP, 1);
    Dia = cell(Sim_Params.N_VP, 1); 
    dTS = cell(Sim_Params.N_VP, 1);
    total_indts_rows = Sim_Params.N_VP * Sim_Params.N_tmr * n_timepoints;
    IndTS = nan(total_indts_rows, 5);
    indts_idx = 1; % Counter to track where to write in IndTS
    time_series_data_for_all_pts = struct();
    for i = 1:length(Sim_Params.StatesToLog)
        var_name = sprintf('sim_%s', Sim_Params.StatesToLog{i});
        time_series_data_for_all_pts.(var_name) = NaN(Sim_Params.N_VP, numel(Sim_Params.Time_Sim));
    end  
           
    % Generate VPOP
    Dose_nanomoles_array = nan(Sim_Params.N_VP,1);
    for patient_id = 1:Sim_Params.N_VP               
        [vp_init_Array, ihc_cat] = generate_virtual_patient(Sim_Params);
        [vp_init_Array.PatientID] = deal(patient_id);

        idx_start = (patient_id-1) * Sim_Params.N_tmr + 1;
        idx_end   = patient_id * Sim_Params.N_tmr;
        % Insert directly into the preallocated slots
        VP_Init_Array(idx_start:idx_end) = vp_init_Array;

        IHC_Category(patient_id) = ihc_cat;
        Dose_nanomoles = (Sim_Params.Dose(patient_id) * vp_init_Array(1).(Sim_Params.parameter_name_for_bodyweight) / (1000 * vp_init_Array(1).(Sim_Params.parameter_name_for_mol_weight_ADC))) * 10^9;  % Convert mpk to nanomole
        Dose_nanomoles_array(patient_id) = Dose_nanomoles;
    end
    VPop_for_simfcn = struct2table(VP_Init_Array);
    VPop_for_simfcn = VPop_for_simfcn{:,3:end};

    % Doses
    num_dose_days = length(Sim_Params.Doses_days_in_a_cycle);
    doses_array = cell(Sim_Params.N_VP,num_dose_days*2);% Multiplication with 2 because of 2 targets-- one for ADC dosing and the other for DAR dosing
    for patient_id = 1:Sim_Params.N_VP
        % ADC Doses
        target_index = 0;
        for dose_day = Sim_Params.Doses_days_in_a_cycle            
            dose_table = table(dose_day,...
                               Dose_nanomoles_array(patient_id),...
                               0,...
                               Sim_Params.Dose_Cycle_length,...
                               round(Sim_Params.StopTime/Sim_Params.Dose_Cycle_length) + 1,...
                               'VariableNames',{'StartTime','Amount','Rate','Interval','RepeatCount'});
            dose_table.Properties.VariableUnits = {'day','nanomole','nanomole/day','day',''};
            target_index = target_index + 1;
            doses_array{patient_id,target_index} = dose_table;
        end
        % DAR dosing. Add Dummy Dosing for ensuring proper DAR resetting at each ADC Dose adminisation
        target_index = num_dose_days;
        for dose_day = Sim_Params.Doses_days_in_a_cycle            
            dose_table = table(dose_day,...
                               1,...
                               0,...
                               Sim_Params.Dose_Cycle_length,...
                               round(Sim_Params.StopTime/Sim_Params.Dose_Cycle_length) + 1,...
                               'VariableNames',{'StartTime','Amount','Rate','Interval','RepeatCount'});
            dose_table.Properties.VariableUnits = {'day','mole','mole/day','day',''};
            target_index = target_index + 1;
            doses_array{patient_id,target_index} = dose_table;
        end
    end    

    % Simulate using simfunction
    % Create simfunction
    param_names = fieldnames(vp_init_Array(1)); param_names = param_names(3:end);
    observables_names = Sim_Params.StatesToLog;
    doses_targets = [repmat({Sim_Params.Dose_target_name},[1,num_dose_days]),...
                     repmat({Sim_Params.Dar_Dose_target_name},[1,num_dose_days])];
    F_simfcn = createSimFunction(model,param_names,observables_names,doses_targets,'UseParallel',Sim_Params.UseParallel_simfcn);
    % Simulate
    simData_series = F_simfcn(VPop_for_simfcn,[],doses_array,Sim_Params.Time_Sim);
    
    % Analysis Loop
    simData_series_index = 0;
    for patient_id = 1:Sim_Params.N_VP

        % Preallocate matrices inside the cells for this specific patient
        Sim_Vol{patient_id} = nan(n_timepoints, Sim_Params.N_tmr);
        Dia{patient_id}     = nan(n_timepoints, Sim_Params.N_tmr);
        dTS{patient_id}     = nan(n_timepoints, Sim_Params.N_tmr);  
        PSA{patient_id}     = nan(n_timepoints, Sim_Params.N_tmr);      
        
        for tumor_id = 1:length(vp_init_Array)   
            simData_series_index = simData_series_index+1;
            simData = simData_series(simData_series_index);
            Time_Sim = simData.Time;

            % Extract data
            TV_all = selectbyname(simData, Sim_Params.parameter_name_for_tumor_volume).Data;
            PSA_all = selectbyname(simData, Sim_Params.parameter_name_for_PSA).Data;
                                            
            % Subsample for efficacy
            idx_efficacy = ismember(Time_Sim, Sim_Params.Time_efficacy); 
            TV = TV_all(idx_efficacy);
            PSA{patient_id}(:,tumor_id) = PSA_all(idx_efficacy);

            Sim_Vol{patient_id}(:,tumor_id) = TV;
            Dia{patient_id}(:,tumor_id) = (6 * Sim_Vol{patient_id}(:,tumor_id) / pi).^(1/3);
            dTS{patient_id}(:,tumor_id) = 100 * (Dia{patient_id}(:,tumor_id) - Dia{patient_id}(1,tumor_id)) / Dia{patient_id}(1,tumor_id);
            lesion_details = [repmat(patient_id, n_timepoints, 1), Sim_Params.Time_efficacy', repmat(tumor_id, n_timepoints, 1), Dia{patient_id}(:,tumor_id), dTS{patient_id}(:,tumor_id)];
            IndTS(indts_idx : indts_idx + n_timepoints - 1, :) = lesion_details;
            indts_idx = indts_idx + n_timepoints;
        end    
               
        % Aggregate at patient level
        VP_Diam = sum(Dia{patient_id}, 2);
        Temp_SLD_per_pts = [repmat(patient_id, n_timepoints, 1), Sim_Params.Time_efficacy', VP_Diam];
        Temp_SLD_per_pts(:,4) = 100 * (Temp_SLD_per_pts(:,3) - Temp_SLD_per_pts(1,3)) / Temp_SLD_per_pts(1,3);
        Temp_SLD_per_pts(:,5) = PSA{patient_id};
        Temp_SLD_per_pts(:,6) = 100 * (Temp_SLD_per_pts(:,5) - Temp_SLD_per_pts(1,5)) / Temp_SLD_per_pts(1,5);

        SLD(1+(patient_id-1)*n_timepoints:patient_id*n_timepoints,:) = Temp_SLD_per_pts;
        
        % Store extra logged time_series
        for i = 1:length(Sim_Params.StatesToLog)
            var_name = sprintf('sim_%s', Sim_Params.StatesToLog{i});
            time_series_data_for_all_pts.(var_name)(patient_id,:) = selectbyname(simData, Sim_Params.StatesToLog{i}).Data'; % Take transpose to ensure row vector
        end
    end
    
    % Post-Processing:  BOR, PFS, and Dropout calculation  
    % Apply Dropout Rules
    Sim_SLD = apply_dropouts(Sim_Params, SLD, IndTS);
    PatientList = (1:Sim_Params.N_VP)';
    % Calculate Best Overall Response (for Waterfall)
    Sim_Waterfall = calculate_bor(Sim_SLD, PatientList);
        
    % Calculate PFS
    Sim_PFS = calculate_pfs(Sim_SLD, Sim_Waterfall);
    
    % Calculate Dropout Times 
    DropoutTime = calculate_dropout_times(Sim_SLD, PatientList);
    
    % Calculate Dropout Reasons
    DropoutSummary = calculate_dropout_reasons(Sim_SLD, PatientList);
    
    % Store All Processed Data for Plotting 
    Simulated_Data = struct();
    Simulated_Data.Sim_SLD = Sim_SLD;
    Simulated_Data.IHC_Category = IHC_Category;
    Simulated_Data.DropoutTime = DropoutTime;
    Simulated_Data.DropoutSummary = DropoutSummary;
    Simulated_Data.Sim_PFS = Sim_PFS;
    Simulated_Data.Sim_Waterfall = Sim_Waterfall;
    
    
    % Bundle Final Results
    Simulation_Results = struct();
    Simulation_Results.AllPlots = Simulated_Data;
    Simulation_Results.Dose = Sim_Params.Dose;
    Simulation_Results.time_series_data_for_all_pts = time_series_data_for_all_pts;
    Simulation_Results.VP_Init_Array = struct2table(VP_Init_Array);
end
 

function [vp_init_Array, ihc_cat] = generate_virtual_patient(Sim_Params)
    % generate_virtual_patient
    %
    % Generates parameters for a single virtual patient with N_tmr tumors.
    %
    % INPUTS:
    %   Sim_Params : Struct containing simulation settings
    %
    % OUTPUT:
    %   vp_init_Array : Struct Array of size [1, N_tmr]
    %                   Fields: .PatientID, .TumorID, .Param1, .Param2...
    %   ihc_cat       : Integer (IHC Category)
    
    params = Sim_Params.param;
    % Get List of Parameter Names and Heterogeneity Settings
    all_param_names = fieldnames(params)';
    hetero_list = Sim_Params.Param_names_for_intrapatient_heterogeneity;
    intra_sigma = Sim_Params.Intrapatient_Heterogeneity;
    n_tumors = Sim_Params.N_tmr;

    % Initialize Temporary Storage for Patient-Level Base Values
    patient_base_vals = struct();
    ihc_cat = 3; 

    % Generate the "Base" phenotype for the patient (Population Variability)
    for i = 1:length(all_param_names)
        p_name = all_param_names{i};
        
        % Extract Mu and Sig
        p_mu = params.(p_name).mu;
        p_sig = params.(p_name).sig;
        
        % Generate the patient-specific base value (Log-Normal Distribution)
        if strcmp(p_name,Sim_Params.parameter_name_for_hscore)
            patient_base_vals.(p_name) = betarnd(5, 2.5, [1, 1]) * 300;
        else
            patient_base_vals.(p_name) = p_mu * exp(randn(1,1) * p_sig);
        end
    end
    vp_init_Array = repmat(cell2struct(cell(numel(all_param_names)+2, 1), [{'PatientID', 'TumorID'}, all_param_names]), 1, n_tumors);
    % Loop through Tumors to generate the Struct Array
    for t_id = 1:n_tumors
        
        % 1. Assign Metadata Fields First (to keep them at the left/top of struct)
        vp_init_Array(t_id).PatientID = NaN;
        vp_init_Array(t_id).TumorID = t_id;
        
        % 2. Loop through model parameters
        for i = 1:length(all_param_names)
            p_name = all_param_names{i};
            base_val = patient_base_vals.(p_name);
            
            % Check Heterogeneity Logic
            if ismember(p_name, hetero_list)
                % Apply Intra-patient heterogeneity (Tumor Variability)
                final_val = base_val * exp(randn(1,1) * intra_sigma);
            else
                % Use fixed patient base value
                final_val = base_val;
                
                
            end
            
            % Assign value to the Struct Array
            % This dynamically adds the field 'p_name' to the struct at index 't_id'
            vp_init_Array(t_id).(p_name) = final_val;
        end
    end
    % IHC Category Logic (calculated on the fixed base value)
    mean_h_score = mean([vp_init_Array(:).(Sim_Params.parameter_name_for_hscore)]);
    if mean_h_score >= 1e-4 && mean_h_score <= 75
        ihc_cat = 0;
    elseif mean_h_score > 75 && mean_h_score <= 150
        ihc_cat = 1;
    elseif mean_h_score > 150 && mean_h_score <= 225
        ihc_cat = 2;
    elseif mean_h_score > 225 
        ihc_cat = 3;
    end

    % Transpose to make it a column vector (N_tumors x 1) if preferred, 
    % currently it generates 1 x N_tumors. 
    % Standard practice is often row-based, so this is likely fine.
end
 
function Sim_SLD = apply_dropouts(Sim_Params, SLD, IndTS)
    %
    % Applies dropout logic to simulation results
    % INPUT : 
    %       Sim_Params:         Information about the crucial parameters to use while simulation
    %       IndTS:              Individual Tumor time series data
    %       SLD:                Time-series of the SLD as simulated
    % OUTPUT :
    %       Sim_SLD : The array containing complete Simulated SLD results with dropout informations
    %

    Time_efficacy = Sim_Params.Time_efficacy;
    N_VP= Sim_Params.N_VP;
    NM_Params = Sim_Params.NM_Params;
    NonSpecific_dropouts = Sim_Params.NonSpecific_dropouts;
    dropout_Switch = Sim_Params.dropout_Switch;

    Sim_SLD = zeros(size(SLD, 1), 11); 
    row_pointer = 1;
    PatientList = 1:N_VP; % Get list of all patient IDs

    for p_id = 1:N_VP
        patient_id = PatientList(p_id);
        lesion_level_for_drop_out = IndTS(IndTS(:,1) == patient_id, :);
        current_patient_data = SLD(SLD(:,1) == patient_id, :);
        n_rows = size(current_patient_data, 1);        
        % Create the working matrix: 6 columns of data + 5 columns of zeros
        % This ensures columns 7-11 exist before we loop, preventing resizing
        patient_level_for_drop_out = [current_patient_data, zeros(n_rows, 5)];
        for tumor_id = 2:length(Time_efficacy)
            % Target Lesion PD
            patient_level_for_drop_out(tumor_id,7) = double(patient_level_for_drop_out(tumor_id,4) > 20);
            
            % Non-Target Lesion PD
            Deltas = lesion_level_for_drop_out(lesion_level_for_drop_out(:,2) == Time_efficacy(tumor_id), 5);
            AbsChanges = lesion_level_for_drop_out(lesion_level_for_drop_out(:,2) == Time_efficacy(tumor_id), 4) - ...
                lesion_level_for_drop_out(lesion_level_for_drop_out(:,2) == 0, 4);
            patient_level_for_drop_out(tumor_id,8) = any(Deltas >= 40 & AbsChanges > 5);
            
            % New Metastases
            diam = patient_level_for_drop_out(tumor_id,3);
            
            if dropout_Switch
                probNM = (NM_Params(2) + NM_Params(1)) * (diam^NM_Params(4)) / ...
                ((NM_Params(3)^NM_Params(4)) + (diam^NM_Params(4)));

                patient_level_for_drop_out(tumor_id,9) = rand < probNM;
                
                % Nonspecific dropouts
                if rand < NonSpecific_dropouts
                    patient_level_for_drop_out(tumor_id,10) = 1;
                else
                    patient_level_for_drop_out(tumor_id,10) = 0;
                end
            else
                patient_level_for_drop_out(tumor_id,9) = 0;
                patient_level_for_drop_out(tumor_id,10) = 0;
            end
            
            % Total PD flags
            patient_level_for_drop_out(tumor_id, 11) = sum(patient_level_for_drop_out(tumor_id,7:10));
        end
        % Determine how many rows we are adding
        n_rows = size(patient_level_for_drop_out, 1);
        
        % Insert into the preallocated slots
        Sim_SLD(row_pointer : row_pointer + n_rows - 1, :) = patient_level_for_drop_out;
        
        % Update the pointer for the next iteration
        row_pointer = row_pointer + n_rows;
    end
end
 
function Sim_Waterfall = calculate_bor(Sim_SLD, PatientList)
    %
    % Calculates the Best Overall Response (BOR) for the waterfall plot
    % INPUT : 
    %       Sim_SLD : The array containing complete Simulated SLD results with dropout informations
    %       PatientList : A simple column vector = 1:N_VP  
    %
    % OUTPUT :
    %       Sim_Waterfall : Array containing BOR for each patient along with the corresponding timepoint
    %

    nPatients = numel(PatientList);
    Sim_Waterfall = zeros(nPatients, 5);
    Sim_Waterfall(:,1) = PatientList;
    
    for patient_id = 1:nPatients
        subj = Sim_SLD(Sim_SLD(:,1) == PatientList(patient_id), :);
        pd_idx = find(subj(:,11) >= 1, 1);
        
        if isempty(pd_idx)
            best_change_TV = min(subj(subj(:,2) ~= 0, 4));
            best_time_TV = subj(subj(:,4) == best_change_TV, 2);
            best_change_PSA = min(subj(subj(:,2) ~= 0, 6));
            best_time_PSA = subj(subj(:,6) == best_change_PSA, 2);
        else
            pre_pd_subj = subj(1:pd_idx, :);
            best_change_TV = min(pre_pd_subj(pre_pd_subj(:,2) ~= 0, 4));
            best_time_TV = pre_pd_subj(pre_pd_subj(:,4) == best_change_TV, 2);
            best_change_PSA = min(pre_pd_subj(pre_pd_subj(:,2) ~= 0, 6));
            best_time_PSA = pre_pd_subj(pre_pd_subj(:,6) == best_change_PSA, 2);
        end
        if isempty(best_time_TV)
            best_time_TV = 0;
        end
        if isempty(best_time_PSA)
            best_time_PSA = 0;
        end
        Sim_Waterfall(patient_id,2) = best_change_TV;
        Sim_Waterfall(patient_id,3) = best_time_TV(1);
        Sim_Waterfall(patient_id,4) = best_change_PSA;
        Sim_Waterfall(patient_id,5) = best_time_PSA(1);
    end
end
 
function Sim_PFS = calculate_pfs(Sim_SLD, Sim_Waterfall)
    %
    % Calculates Progression-Free Survival (PFS) data
    % INPUT : 
    %       Sim_SLD : The array containing complete Simulated SLD results with dropout informations
    %       Sim_Waterfall : Array containing BOR for each patient along with the corresponding timepoint
    %
    % OUTPUT :
    %       Sim_PFS : Array containing PFS related data for the simulated population
    %

    Sim_PFS = [Sim_Waterfall(:,1), zeros(size(Sim_Waterfall,1), 5)];
    for i = 1:size(Sim_PFS,1)
        subj = Sim_SLD(Sim_SLD(:,1) == Sim_PFS(i,1), :);
        pd_idx = find(subj(:,11) >= 1, 1);
        
        if isempty(pd_idx)
            Sim_PFS(i,2) = subj(end,2);  % Last time point
            Sim_PFS(i,3) = 0;            % Censored
        else
            Sim_PFS(i,2) = subj(pd_idx,2);  % PD time
            Sim_PFS(i,3) = 1;              % Event occurred
        end
    end
    
    Sim_PFS = sortrows(Sim_PFS, 2);
    Sim_PFS(:,4) = size(Sim_PFS,1):-1:1;
    Sim_PFS(:,5) = 1 - Sim_PFS(:,3) ./ Sim_PFS(:,4);
    Sim_PFS(1,6) = 1;
    for i = 2:size(Sim_PFS,1)
        Sim_PFS(i,6) = Sim_PFS(i-1,6) * Sim_PFS(i,5);
    end
    Sim_PFS = [[0 0 0 0 0 1]; Sim_PFS];  % Add (0,1) starting point
end
 
function DropoutTime = calculate_dropout_times(Sim_SLD, PatientList)
    %
    % Identifies the time of dropout for each patient
    % INPUT : 
    %       Sim_SLD : The array containing complete Simulated SLD results with dropout informations
    %       PatientList : A simple column vector = 1:N_VP  
    %
    % OUTPUT : 
    %       DropoutTime : Dropout time for all patients
    %
    
    nPatients   = numel(PatientList);
    DropoutTime = zeros(nPatients,1);
    for i = 1:nPatients
        pid = PatientList(i);
        data_i = Sim_SLD(Sim_SLD(:,1) == pid, :);
        idx_dropout = find(data_i(:,11) >= 1, 1, 'first');
        
        if ~isempty(idx_dropout)
            DropoutTime(i) = data_i(idx_dropout, 2);
        else
            DropoutTime(i) = max(data_i(:,2));  % No dropout
        end
    end
end
 
function DropoutSummary = calculate_dropout_reasons(Sim_SLD, PatientList)
    % 
    % Identifies the reason for dropout for each patient
    % INPUT : 
    %       Sim_SLD : The array containing complete Simulated SLD results with dropout informations
    %       PatientList : A simple column vector = 1:N_VP  
    %
    % OUTPUT : 
    %       DropoutSummary : Table containing dropout times and dropout reasons for all patients
    %
    
    nPatients   = numel(PatientList);
    DropoutSummary = table(PatientList, zeros(nPatients,1), zeros(nPatients,1), ...
        'VariableNames', {'PatientID','DropoutTime','ReasonCode'});
    
    % ReasonCode: 5=TargetPD, 6=NonTargetPD, 7=NewMet, 8=Nonspecific, 0=No dropout
    for i = 1:nPatients
        pid = PatientList(i);
        data_i = Sim_SLD(Sim_SLD(:,1) == pid,:);
        idx_dropout = find(data_i(:,11) >= 1, 1, 'first');
        
        if ~isempty(idx_dropout)
            DropoutSummary.DropoutTime(i) = data_i(idx_dropout,2);
            reasons = data_i(idx_dropout,7:10);
            if reasons(1)==1
                DropoutSummary.ReasonCode(i) = 5; % Target lesion PD
            elseif reasons(2)==1
                DropoutSummary.ReasonCode(i) = 6; % Non-target lesion PD
            elseif reasons(3)==1
                DropoutSummary.ReasonCode(i) = 7; % New metastases
            elseif reasons(4)==1
                DropoutSummary.ReasonCode(i) = 8; % Nonspecific dropout
            end
        else
            DropoutSummary.DropoutTime(i) = max(data_i(:,2)); % completed
            DropoutSummary.ReasonCode(i) = 0; % no dropout
        end
    end
end
 
function plot_waterfall(Simulation_Results, LesionWaterfall, PSAWaterfall)
    %
    % plot_waterfall: Generates a pooled waterfall plot from simulation results.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %       Waterfall:          Processed data containing information about the waterfall plot from the raw data
    %
    
    gem_hex_codes = rgb2hex(orderedcolors("gem"));

    simulated_waterfall_Data = Simulation_Results.AllPlots.Sim_Waterfall;

    % Sort by Best Overall Response (BOR)
    simulated_waterfall_Data = sortrows(simulated_waterfall_Data, -2);
    N_VP = size(simulated_waterfall_Data,1);   
    sim_waterfall = simulated_waterfall_Data(:,2);    
    psawaterfall = simulated_waterfall_Data(end-height(PSAWaterfall)+1:end,4);
    psawaterfall = sort(psawaterfall, 'descend');
    figure;
    hold on;    
    % Create the bar plot for the Data
    data_prc_change=sort(LesionWaterfall.BOR,'descend');
    bar(data_prc_change, 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    
    % Create a curve for simulated waterfall plot
    plot(1:N_VP, sim_waterfall, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0,'DisplayName', 'Simulated');

    % Add lines to define Stable Disease zone (SD)
    yline(20, '--k','HandleVisibility','off'); yline(-30, '--k','HandleVisibility','off');
    
    % Formatting
    xlabel('Virtual Patients (Sorted by Response)','FontWeight','bold');
    ylabel('BOR (%)','FontWeight','bold');
    global_min = min([min(sim_waterfall(:)), min(data_prc_change(:))]);
    global_max = max([max(sim_waterfall(:)), max(data_prc_change(:))]);
    
    y_min = max(-105, global_min - 10); 
    y_max = global_max + 10;
    ylim([y_min, y_max]);
    xlim([0.5,N_VP+0.5]);
    grid off; box on; set(gca,'LineWidth',1.25);
    pbaspect([3 2 1]);
    fontsize(15,'points');
    legend('Location','southwest','Box','on','FontSize',15);
    title('Calibrating Tumor Waterfall Plot','FontWeight','bold','FontSize',18);
    fontname('Consolas');
    hold off;

    figure;
    hold on;    
    % Create the bar plot for the Data
    data_psa_change=sort(PSAWaterfall.BOR,'descend');
    bar(data_psa_change, 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    
    % Create a curve for simulated waterfall plot
    plot(1:height(data_psa_change), psawaterfall, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0,'DisplayName', 'Simulated');

    % yline(20, '--k','HandleVisibility','off'); yline(-30, '--k','HandleVisibility','off');
    
    % Formatting
    xlabel('Virtual Patients (Sorted by Response)','FontWeight','bold');
    ylabel('BOR (%)','FontWeight','bold');
    global_min = min([min(psawaterfall(:)), min(data_psa_change(:))]);
    global_max = max([max(psawaterfall(:)), max(data_psa_change(:))]);
    
    y_min = max(-105, global_min - 10); 
    y_max = global_max + 10;
    ylim([y_min, y_max]);
    xlim([0.5,height(data_psa_change)+0.5]);
    grid off; box on; set(gca,'LineWidth',1.25);
    pbaspect([3 2 1]);
    fontsize(15,'points');
    legend('Location','southwest','Box','on','FontSize',15);
    title('Calibrating PSA Waterfall Plot','FontWeight','bold','FontSize',18);
    fontname('Consolas');
    hold off;
end
 
function plot_pfs(Simulation_Results)
    %
    % plot_pfs: Generates Kaplan-Meier PFS plots along with 95% confidence interval for each dose.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %

    Sim_PFS = Simulation_Results.AllPlots.Sim_PFS;
    gem_hex_codes = rgb2hex(orderedcolors("gem"));
     
    pfs_times = Sim_PFS(1:end, 2);
    pfs_event_status = Sim_PFS(1:end, 3); % 1 for event, 0 for censored
    
    % Convert event status to logical censoring vector for ecdf
    % ecdf expects 'true' for censored observations
    censoring_flags = ~logical(pfs_event_status);
    
    % Calculate the survivor function and confidence bounds
    [f, x, flo, fup] = ecdf(pfs_times, 'Censoring', censoring_flags, 'Function', 'survivor', 'Alpha', 0.05, 'Bounds', 'on');
    
    % Ensure bounds start at 1 and end consistently
    flo(1) = 1;
    fup(1) = 1;
    if f(end) == 0
        flo(end) = 0;
        fup(end) = 0;
    end
    
    % Create stairstep coordinates for the confidence bounds
    if length(x)>2 %if length of x is 2, then repelem(x(2:end), 2) gives a row vector, leading to an error due to mismatch is array size while vertically concatenating
        x_stairs = [x(1); repelem(x(2:end), 2)];
        lower_bound_stairs = [repelem(flo(1:end-1), 2); flo(end)];
        upper_bound_stairs = [repelem(fup(1:end-1), 2); fup(end)];
    else
        x_stairs = [x(1); repelem(x(2:end), 2)'];
        lower_bound_stairs = [repelem(flo(1:end-1), 2)'; flo(end)];
        upper_bound_stairs = [repelem(fup(1:end-1), 2)'; fup(end)];
    end
    figure;
    hold on;
    % Plot Kaplan-Meier curve
    stairs(x, f * 100, 'Color', gem_hex_codes(7), 'LineWidth', 2,'DisplayName','PFS Curve');

    % Plot the 95% Confidence Limits region
    fill([x_stairs; flipud(x_stairs)], [lower_bound_stairs * 100; flipud(upper_bound_stairs * 100)], hex2rgb(gem_hex_codes(7)), 'FaceAlpha', 0.15, 'EdgeColor', 'none','DisplayName','95% CI');
    
    % Formatting
    xlabel('Time (days)','FontWeight','bold');
    ylabel('Survival Probability(%)','FontWeight','bold');
    ylim([0 102]);
    grid off; box on; set(gca,'LineWidth',1.25);
    pbaspect([3 2 1]);
    fontsize(15,'points'); fontname('Consolas');
    legend('Box','on','Location', 'northeast');
    title('Progression-Free Survival (PFS)','FontWeight','bold','Fontsize', 18);
    hold off;
end

function plot_swimmer(Simulation_Results)
    %
    % plot_swimmer: Generates a swimmer plot for each dose.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %
   
    % Define colors for reasons
    reasonColors = containers.Map({0, 5, 6, 7, 8}, {
        [0.3 0.7 0.3];  % Continue = green
        [0.8 0.5 0.0];  % Target PD = orange
        [0.0 0.4 0.9];  % Non-target PD = blue
        [0.9 0.2 0.2];  % New metastases = red
        [0.5 0.5 0.5]   % Nonspecific = gray
        });
    reasonLabels = {'Continued', 'Target PD', 'Non-target PD', 'New Metastases', 'Nonspecific'};
    reasonCodes  = [0, 5, 6, 7, 8];

    % Extract summary table
    DropoutSummary  = Simulation_Results.AllPlots.DropoutSummary;
    nPatients = size(DropoutSummary, 1);
    
    % Sort patients by dropout time
    DropoutSummary = sortrows(DropoutSummary, "DropoutTime", "descend");
    DropoutSummary.PatientOrder = (1:nPatients)'; % new order for plotting
    
    % Create a new figure for each dose
    figure; 
    hold on;
    
    for i = 1:nPatients
        t_end = DropoutSummary.DropoutTime(i);
        reason = DropoutSummary.ReasonCode(i);
        order = DropoutSummary.PatientOrder(i);
        
        % Plot horizontal bar
        barh(order, t_end, 0.8, 'FaceColor', reasonColors(reason), 'EdgeColor', 'k');
    end
    
    % Formatting
    ylim([0, size(DropoutSummary.DropoutTime,1) + 1]);
    xlim([0, 370]);
    xlabel('Time (days)','FontWeight','bold');
    ylabel('Virtual Patients','FontWeight','bold');
    fontsize(15,'points');
    title('Simulated Swimmer Plot','FontWeight','bold','FontSize',16);
    set(gca, 'YDir', 'reverse'); 
    grid off; box on; set(gca, 'LineWidth', 1.25);
    pbaspect([3 2 1]);
    fontname('Consolas');
    
    % Legend
    h_leg = zeros(1,numel(reasonCodes));
    for k = 1:numel(reasonCodes)
        h_leg(k) = barh(nan, nan, 0.6, 'FaceColor', reasonColors(reasonCodes(k)));
    end
    legend(h_leg, reasonLabels, 'Location', 'southeast');
    hold off;
    
    % Display dropout percentages
    counts = zeros(size(reasonCodes));
    for r = 1:numel(reasonCodes)
        counts(r) = sum(DropoutSummary.ReasonCode == reasonCodes(r));
    end
    percentages = 100 * counts / nPatients;
    disp('Dropout Summary :');
    disp(table(reasonLabels', counts', percentages', 'VariableNames', {'Reason', 'Count', 'Percent'}));
    
end

function plot_spider(Simulation_Results)
    %
    % plot_spider: Generates a spider plot (tumor trajectories) for each dose.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %


    % Define colors for IHC categories
    % IHC_colors = [0.3 0.3 0.3;   % Gray for IHC0
    %               0.0 0.6 0.9;   % Blue for IHC1+
    %               0.3 0.8 0.3;   % Green for IHC2+
    %               0.9 0.2 0.2];  % Red for IHC3+

    AllPlots = Simulation_Results.AllPlots;
    % Extract data for this dose
    Sim_SLD = AllPlots.Sim_SLD;
    DropoutTime = AllPlots.DropoutTime;
    % IHC_Category = AllPlots.IHC_Category;
    
    PatientList = unique(Sim_SLD(:,1));
    nPatients   = numel(PatientList);        
    response_at_1st_visit = numel(nPatients, 1);
    
    figure; 
    hold on; 

    for i = 1:nPatients
        pid = PatientList(i);
        data_i = Sim_SLD(Sim_SLD(:,1) == pid, :);
        t_drop = DropoutTime(i);
        data_cut = data_i(data_i(:,2) <= t_drop, :);
        response_at_1st_visit(i) = data_cut(2,4);

        % Find this patient's IHC category
        % IHC_class = IHC_Category(pid); 

        gem_hex_codes = rgb2hex(orderedcolors("gem"));
        plot(data_cut(:,2), data_cut(:,4), '-', 'Color',gem_hex_codes(1), 'LineWidth', 2);
    end

    % Formatting
    % h_leg = zeros(4,1);
    % for k = 1:4
    %     h_leg(k) = plot(nan, nan, '-', 'Color', IHC_colors(k,:), 'LineWidth', 4); % bogus plots to ensure the IHC Colors are shown only once no matter how many curves are there on the spider plot
    % end
    
    xlabel('Time (days)', 'FontWeight','bold');
    ylabel('% Change in SLD', 'FontWeight','bold');
    fontsize(15, 'points');
    text(10, 80 , sprintf('Mean Depth of Response at 1st Clinical Visit : %.2f', mean(response_at_1st_visit)), 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'left');
    yline(0,'--k');
    title('Simulated Spider Plot');
    ylim([-100 100]);
    xlim([0, 370]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([ 3 2 1]);
    % legend(h_leg, {'IHC0','IHC1+','IHC2+','IHC3+'}, 'Location','northeast');
    fontname('Consolas');
    hold off;

end

function plot_psa_spider(Simulation_Results)
    %
    % plot_psa_spider: Generates a spider plot (psa trajectories) for each dose.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %


    % % Define colors for IHC categories
    % IHC_colors = [0.3 0.3 0.3;   % Gray for IHC0
    %               0.0 0.6 0.9;   % Blue for IHC1+
    %               0.3 0.8 0.3;   % Green for IHC2+
    %               0.9 0.2 0.2];  % Red for IHC3+

    AllPlots = Simulation_Results.AllPlots;
    % Extract data for this dose
    Sim_SLD = AllPlots.Sim_SLD;
    DropoutTime = AllPlots.DropoutTime;
    % IHC_Category = AllPlots.IHC_Category;
    
    PatientList = unique(Sim_SLD(:,1));
    nPatients   = numel(PatientList);        
    % response_at_1st_visit = numel(nPatients, 1);
    
    figure; 
    hold on; 

    for i = 1:nPatients
        pid = PatientList(i);
        data_i = Sim_SLD(Sim_SLD(:,1) == pid, :);
        t_drop = DropoutTime(i);
        data_cut = data_i(data_i(:,2) <= t_drop, :);
        % response_at_1st_visit(i) = data_cut(2,4);
        data_cut(:,6) =  min(100,data_cut(:,6));
        % Find this patient's IHC category
        % IHC_class = IHC_Category(pid); 
        
        gem_hex_codes = rgb2hex(orderedcolors("gem"));
        if min(data_cut(2:end,4))<=0
            plot(data_cut(:,2), data_cut(:,6), '-', 'Color', gem_hex_codes(1), 'LineWidth', 2);
        end
        % plot(data_cut(:,2), data_cut(:,6), '-', 'Color', IHC_colors(IHC_class+1,:), 'LineWidth', 2);
    end

    % % Formatting
    % h_leg = zeros(4,1);
    % for k = 1:4
    %     h_leg(k) = plot(nan, nan, '-', 'Color', IHC_colors(k,:), 'LineWidth', 4); % bogus plots to ensure the IHC Colors are shown only once no matter how many curves are there on the spider plot
    % end
    
    xlabel('Time (days)', 'FontWeight','bold');
    ylabel('% Change in PSA', 'FontWeight','bold');
    fontsize(15, 'points');
    yline(0,'--k');
    title('Simulated PSA trajectories');
    ylim([-100 100]);
    xlim([0, 370]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([ 3 2 1]);
    % legend(h_leg, {'IHC0','IHC1+','IHC2+','IHC3+'}, 'Location','northeast');
    fontname('Consolas');
    hold off;

end

function plot_observables(Simulation_Results,Sim_Params)
    %
    % plot__observables: Generates plots for all the states_to_log as
    % specified at the top of the script, except Tumor Volume. You can log
    % whatever states you want to log and it will plot it
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated output
    %       Sim_Params:         Information about the crucial parameters to use while simulation
    %
    
    
    
    Time_Sim = Sim_Params.Time_Sim(1:(Sim_Params.Scan_Interval+1));
    % Extract data from the results struct
    ts = Simulation_Results.time_series_data_for_all_pts;

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    nAnalytes  = length(fields(ts));
    figure;
    hold on;

    for i = 1:length(Sim_Params.StatesToLog) 

        subplot(2, ceil(nAnalytes/2), i); 
        var_name = sprintf('sim_%s', Sim_Params.StatesToLog{i});
        working_ts = ts.(var_name);
        working_ts = working_ts(:,1:Sim_Params.Scan_Interval+1);
        
        % Plot multiple simulated trajectory (N_VP curves)
        plot(Time_Sim, working_ts, 'Color', gem_hex_codes(i), 'LineWidth', 2); %if this line gives error due to gem_hex_codes, just remove the color property from here. It gives error because you are logging more than 7 states.
        
        title(sprintf('Plotting Time Series %d',i-1));

        xlabel('Time (days)'); ylabel(Sim_Params.StatesToLog{i},'Interpreter','none');
        xlim([-1 43]);
        ylim([0.1*min(working_ts(:)), 10*max(working_ts(:))]);
        box on; set(gca, 'YScale', 'log');
        
        fontname('Consolas');
        pbaspect([3 2 1]);
        hold off;
        fontsize(12,'points');
    end
end
 
function PatientTable = summarize_patient_response(Simulation_Results)
    %
    % summarize_patient_response: Collects and prints patient-level summary (PID, Dose, H-score, BestResponse)
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %
    % OUTPUTS:
    %       PatientTable:       A table with columns {PID, Dose_mpk, Hscore, BestResponse}. This is printed and returned
    %
        
    % Extract data from results struct
    Dose = Simulation_Results.Dose;
    AllPlots = Simulation_Results.AllPlots;
    VP_Init_Array = Simulation_Results.VP_Init_Array;
    Sim_SLD      = AllPlots.Sim_SLD;
    DropoutTime  = AllPlots.DropoutTime;
        
    PatientList = unique(Sim_SLD(:,1));
    nPatients   = numel(PatientList);
    PatientSummary = zeros(nPatients,4);  % Columns will be [PID, Dose, Hscore, BestResp]

    for i = 1:nPatients
        pid = PatientList(i);
        
        % corresponding antigen % positivity / H-score
        lesions = VP_Init_Array(VP_Init_Array.PatientID == pid,:);
        Hscore = mean(lesions.H_Score); 
        
        % trajectory truncated at dropout
        data_i  = Sim_SLD(Sim_SLD(:,1) == pid, :);
        t_drop  = DropoutTime(i);
        data_cut = data_i(data_i(:,2) <= t_drop, :);
        
        % Best response = minimum % change in SLD
        % Ensure data_cut(2:end,4) is not empty
        if size(data_cut, 1) > 1
            best_resp = min(data_cut(2:end, 4));
        else
            best_resp = data_cut(1, 4); % Use baseline if no follow-up
        end
        
        % add row
        PatientSummary(i,:) = [pid, Dose, Hscore, best_resp];
    end
    
    % Convert to table
    PatientTable = array2table(PatientSummary, ...
        'VariableNames', {'PID','Dose_mpk','Hscore','BestResponse'});
    disp(PatientTable);
end

function Plot_waterfall_for_all_seed(psa_waterfall_for_all_seed, tumor_waterfall_for_all_seed, plot_data, LesionWaterfall, PSAWaterfall)
    figure; 
    hold on;
    box on; set(gca,'LineWidth', 1.25);

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    if plot_data
        bar(sort(LesionWaterfall.BOR,'descend'), 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    end
    
    meanVals = mean(tumor_waterfall_for_all_seed, 1);
    lower_bound_95 = prctile(tumor_waterfall_for_all_seed, 2.5, 1); 
    upper_bound_95 = prctile(tumor_waterfall_for_all_seed, 97.5, 1);
   
    x = 1:numel(meanVals);
    
    % Shaded region (min–max)
    fill([x fliplr(x)], [lower_bound_95 fliplr(upper_bound_95)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(7), ...
        'HandleVisibility','off');
    
    % Overlay lines
    plot(x, meanVals, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Mean Simulated Waterfall');
    
    % Calculate % of patients with response ≤ -30% 
    ORR_val = sum(meanVals <= -30) / numel(meanVals) * 100;
    % Mention ORR on the plot
    text_loc_x = max(x);   
    text_loc_y = max(upper_bound_95)-1; 
    text_str = sprintf('Mean ORR is: %.1f%%', ORR_val);
    text(text_loc_x, text_loc_y, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');

    % Cosmetics
    ylabel('% Change in SLD',FontWeight="bold");
    xlabel('Patient No.',FontWeight="bold");
    fontsize(15,'points');
    title('I-Dxd PD Calibration',FontSize = 18,FontWeight="bold");
    legend('Box','off', 'Location', 'southwest','FontSize',13);
    % ylim([-100,100]);
    % xlim([0.5, Sim_Params.N_VP+0.5]);
            
    pbaspect([3 2 1]);
    fontname('Consolas');


    % PSA Waterfall plotting
    figure; 
    hold on;
    box on; set(gca,'LineWidth', 1.25);

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    if plot_data
        bar(sort(PSAWaterfall.BOR,'descend'), 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    end
    
    psa_waterfall_for_all_seed = min(100, psa_waterfall_for_all_seed);
    meanVals = mean(psa_waterfall_for_all_seed, 1);
    lower_bound_95 = prctile(psa_waterfall_for_all_seed, 2.5, 1); 
    upper_bound_95 = prctile(psa_waterfall_for_all_seed, 97.5, 1);
   
    x = 1:numel(meanVals);
    
    % Shaded region (min–max)
    fill([x fliplr(x)], [lower_bound_95 fliplr(upper_bound_95)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(7), ...
        'HandleVisibility','off');
    
    % Overlay lines
    plot(x, meanVals, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Mean PSA Waterfall');
   
    % Cosmetics
    ylabel('% Change in Serum PSA',FontWeight="bold");
    xlabel('Patient No.',FontWeight="bold");
    fontsize(15,'points');
    title('I-Dxd PSA Calibration',FontSize = 18,FontWeight="bold");
    legend('Box','off', 'Location', 'southwest','FontSize',13);
    ylim([-100,100]);
    % xlim([0.5, Sim_Params.N_VP+0.5]);
            
    pbaspect([3 2 1]);
    fontname('Consolas');
end

function Plot_PFS_for_all_seed(Sim_Params, pfs_for_all_seed, plot_95_CI, plot_choice, seed_indices)
    
    % Set default plot_choice 
    if nargin < 5 || isempty(plot_choice)
        plot_choice = "all seeds";
    end
    if strcmpi(plot_choice,'certain seeds') & isempty(seed_indices)
        error('If you specify the plot choice as "certain seeds", then you must specify the seed_indices as the last argument in the function Plot_PFS_for_all_seed.');
    end
    % Get the color palette
    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    nColors = size(gem_hex_codes, 1);
    if strcmp(plot_choice, "all seeds")
        % PLOT ALL SEEDS
        
        % Get dimensions from the main cell array
        nSeeds = length(pfs_for_all_seed);
        
        % Loop over each DOSE - one figure per dose
    
        figure; 
        hold on;
        grid off;
        box on; set(gca, 'LineWidth', 1.25);
        
        % Loop over each SEED for the current dose
        for sIndex = 1:nSeeds
            
            Sim_PFS = pfs_for_all_seed{sIndex};
            
            % Extract data for ecdf
            pfs_times = Sim_PFS(1:end, 2);
            pfs_event_status = Sim_PFS(1:end, 3);
            censoring_flags = ~logical(pfs_event_status);
            
            % Calculate survivor function and bounds
            [f, x, flo, fup] = ecdf(pfs_times, ...
                                    'Censoring', censoring_flags, ...
                                    'Function', 'survivor', 'Alpha', 0.05, 'Bounds', 'on');
            
            % The following manual adjustment part is necessary to ensure the 95% range is plotted properly
            flo(1) = 1;
            fup(1) = 1;
            if f(end) == 0
                flo(end) = 0;
                fup(end) = 0;
            end

            x(end+1) = pfs_times(end); %#ok<AGROW>
            f(end+1) = f(end); %#ok<AGROW>
            
            % Create stairstep coordinates for confidence bounds
            if length(x)>2 %if length of x is 2, then repelem(x(2:end), 2) gives a row vector, leading to an error due to mismatch is array size while vertically concatenating
                x_stairs = [x(1); repelem(x(2:end), 2)];
                lower_bound_stairs = [repelem(flo(1:end), 2); flo(end)];
                upper_bound_stairs = [repelem(fup(1:end), 2); fup(end)];
            else
                x_stairs = [x(1); repelem(x(2:end), 2)'];
                lower_bound_stairs = [repelem(flo(1:end), 2)'; flo(end)];
                upper_bound_stairs = [repelem(fup(1:end), 2)'; fup(end)];
            end
            
            
            % Get color for this seed (wraps around if nSeeds > nColors)
            color_index = mod(sIndex - 1, nColors) + 1;
            color_sim = gem_hex_codes(color_index); % Use string array indexing

            if plot_95_CI
                % Plot the 95% Confidence Limits region
                fill([x_stairs; flipud(x_stairs)], [lower_bound_stairs * 100; flipud(upper_bound_stairs * 100)], ...
                     hex2rgb(color_sim), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
            end

            % Plot the Kaplan-Meier estimate
            stairs(x, f * 100, 'Color', color_sim, 'LineWidth', 2.5, 'DisplayName', sprintf('Seed %d PFS', sIndex));
            
        end % End of seed loop
        
        % Cosmetics for the DOSE figure
        xlabel('Time (days)',FontWeight="bold");
        ylabel('PFS Probability (%)',FontWeight="bold");
        fontsize(15,'points'); 
        title('Simulated PFS for I-Dxd Preliminary Calibration', FontSize=18, FontWeight="bold");
        ylim([0 102]); xlim([0,366]);
        % legend('Box','on', 'Location', 'northeast', FontSize= 15);
        pbaspect([3 2 1]);
        fontname('Consolas');
        hold off;
        
    elseif strcmp(plot_choice, "certain seeds")
        % PLOT SINGLE SEED
        figure; 
        hold on;
        grid off;
        box on; set(gca, 'LineWidth', 1.25);
        nseeds = numel(seed_indices);
       
        for sIndex=1:nseeds
        
            Sim_PFS = pfs_for_all_seed{sIndex};
            pfs_times = Sim_PFS(1:end, 2);
            pfs_event_status = Sim_PFS(1:end, 3); % 1 for event, 0 for censored
            
            % Convert event status to logical censoring vector for ecdf
            censoring_flags = ~logical(pfs_event_status);
            
            % Calculate the survivor function and confidence bounds
            [f, x, flo, fup] = ecdf(pfs_times, ...
                                    'Censoring', censoring_flags, ...
                                    'Function', 'survivor', 'Alpha', 0.05, 'Bounds', 'on');
            
            % Important manual adjustment to ensure proper plotting of 95%ile range
            flo(1) = 1;
            fup(1) = 1;
            if f(end) == 0
                flo(end) = 0;
                fup(end) = 0;
            else
                % If curve doesn't end at 0, hold last bound
                if length(f) > 1
                    flo(end) = flo(end-1);  
                    fup(end) = fup(end-1);
                end
            end
            
            % Create stairstep coordinates for the confidence bounds
            x_stairs = [x(1); repelem(x(2:end), 2)];
            lower_bound_stairs = [repelem(flo(1:end-1), 2); flo(end)];
            upper_bound_stairs = [repelem(fup(1:end-1), 2); fup(end)];        
            color_sim = gem_hex_codes(sIndex);  %if this line gives error then remove the color input from the below fill and stairs function
            if plot_95_CI
                % Plot the 95% Confidence Limits region
                fill([x_stairs; flipud(x_stairs)], [lower_bound_stairs * 100; flipud(upper_bound_stairs * 100)], ...
                     hex2rgb(color_sim), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
            end
    
            % Plot the Kaplan-Meier estimate (the survival curve)
            stairs(x, f * 100, 'Color', color_sim, 'LineWidth', 2.5, 'DisplayName', 'Simulated PFS');
        end    
        % Cosmetics
        xlabel('Time (days)',FontWeight="bold");
        ylabel('PFS Probability (%)',FontWeight="bold");
        fontsize(15,'points'); 
        title('Simulated PFS for I-Dxd Preliminary Calibration', 'FontSize', 16, FontWeight="bold");
        ylim([0 100]); xlim([0,Sim_Params.StopTime+1]);
        legend('Box','off', 'Location', 'best', FontSize= 12);
        pbaspect([3 2 1]);
        fontname('Consolas');
        hold off;
          
    else
        error('Invalid plot_choice. Must be "all seeds" or "certain seeds".');
    end

end

function Plot_dose_response(Dose, Sim_Params, All_Dose_Waterfall)
    %
    % plot_dose_response: Plots box plot indicating the dose vs response (Best % change in SLD) trend
    %
    % INPUTS:
    %   Dose:                  Array containing the doses for which the script simulates
    %   All_Dose_Waterfall:    The cell array containing responses for all seeds across all doses
    %
    % OUTPUTs:
    %   A box plot showing dose vs response trend
    %
    n_doses = length(Dose);
    n_seeds = size(All_Dose_Waterfall{1},1);
    n_vp = Sim_Params.N_VP;
    allDoseVals = nan(n_vp*n_seeds*n_doses,1);
    allRespVals = nan(n_vp*n_seeds*n_doses,1);
    
    for dose = 1:n_doses
        
        RespVals = All_Dose_Waterfall{dose}(:);
        RespVals = RespVals(:); % flatten to column vector
        
        allRespVals((dose-1)*n_vp*n_seeds+1:dose*n_vp*n_seeds,1) = RespVals;
        allDoseVals((dose-1)*n_vp*n_seeds+1:dose*n_vp*n_seeds,1) = repmat(Dose(dose), numel(RespVals), 1);

    end
    figure;
    boxplot(allRespVals, allDoseVals, 'Symbol', 'o', 'Colors', hex2rgb("#DD5400"), 'Symbol', 'o', 'PlotStyle', 'traditional', 'OutlierSize', 5, 'Positions', unique(Dose), 'Widths', 1);
    set(findobj(gca, 'Type', 'line', '-regexp', 'Tag', 'Whisker|Adjacent Value|Box|Median'), 'LineWidth', 1.5, 'LineStyle', '-');
    xlabel('Dose (mg/kg)', 'FontWeight', 'bold');
    ylabel('Best % Change in SLD', 'FontWeight', 'bold');
    grid off; box on;
    set(gca, 'LineWidth', 1.25, 'FontName', 'Consolas', 'FontSize', 15);
    title('Dose vs Best Response', 'FontSize', 16, 'FontWeight', 'bold');
end

function Sim_Params = update_params_from_estimates(Sim_Params, model_params_e, model)
    % update_tox_params_from_estimates
    %
    % Updates Tox_Params.param based on fitted values provided in a simplified 
    % format: {ParamName, Mu_Value, Sig_Value}.
    %
    % INPUTS:
    %   Tox_Params     : Structure containing .param.(param_name).mu/sig
    %   model_params_e : Nx3 Cell array {'ParamName', mu_val, sig_val}
    %                    Use [] or NaN to skip updating a specific value.
    %   model          : SimBiology model object (needed for default values)
    %
    % OUTPUT:
    %   Tox_Params     : Updated structure

    n_params = size(model_params_e, 1);

    for i = 1:n_params
        p_name  = model_params_e{i, 1};
        mu_in   = model_params_e{i, 2};
        sig_in  = model_params_e{i, 3};
        
        % Check if inputs are valid (not empty and not NaN)
        has_mu  = ~isempty(mu_in) && isnumeric(mu_in) && ~isnan(mu_in);
        has_sig = ~isempty(sig_in) && isnumeric(sig_in) && ~isnan(sig_in);

        % Check if parameter already exists in Tox_Params
        if isfield(Sim_Params.param, p_name)
            % EXISTING PARAMETER: Update provided values
            if has_mu
                Sim_Params.param.(p_name).mu = mu_in;
            end
            if has_sig
                Sim_Params.param.(p_name).sig = sig_in;
            end
        else
            % NEW PARAMETER
            final_mu = []; %#ok<NASGU>
            final_sig = []; %#ok<NASGU>

            % Determine Mu Value
            if has_mu
                final_mu = mu_in;
            else
                % Use Model Value
                pObj = sbioselect(model, 'Name', p_name);
                if isempty(pObj)
                     error('Parameter "%s" is not in Tox_Params and not found in the Model. Cannot determine Mu.', p_name);
                end
                final_mu = pObj.Value;
                fprintf('WARNING : New parameter "%s" added. Mu value missing in input; using model default: %g.', p_name, final_mu);
            end

            % Determine Sig Value
            if has_sig
                final_sig = sig_in;
            else
                % Set Sig to 0
                final_sig = 0;
                fprintf('WARNING : New parameter "%s" added. Sig value missing in input; setting Sig to 0.', p_name);
            end

            % Add to Structure
            Sim_Params.param.(p_name).mu  = final_mu;
            Sim_Params.param.(p_name).sig = final_sig;
        end
    end
end