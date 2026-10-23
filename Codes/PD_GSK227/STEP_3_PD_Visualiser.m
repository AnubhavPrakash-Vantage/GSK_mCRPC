clear; warning off;
load("New Analysis-sept/Lesion_Waterfall_8mpk_all.mat"); 
load("New Analysis-sept/PSA_Waterfall_8mpk_all.mat"); 
load("New Analysis-sept/Mean_spider_clinical_8mpk_all.mat");
load("New Analysis-sept/Median_PSA_spider_clinical_8mpk_all.mat");
load('Sim_Params.mat');                                    % Loads the calibrated and pre-assigned parameters for simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load Project model
proj = sbioloadproject("../../GSK227_CRPC.sbproj");     % This path has to be changed as per the location of the model
model = proj.m1; 

seed = 1:100;                                    % This is the RNG Seed value to which you are simulating
Sim_Params.Dose = [6,8,10,12];        % in mg/kg 
Sim_Params.N_VP = 30;     % Number of patients to simulate is same as that in the data
Sim_Params.variability = 1; 
Sim_Params.StopTime = 340;                       % Simulating for 1 year
Sim_Params.Time_Sim = 0:1:Sim_Params.StopTime;    % for each of the days, time point where the simulations will be carried out
Sim_Params.Dose_Cycle_length = 21;                 % in days, 21 Days for Q3W
Sim_Params.Doses_days_in_a_cycle = 0;              % On what days you want to administer the dosage, it could be [0,7] if you want d1,d8 Q3W kind of dosing
Sim_Params.Dose_target_name = 'Central_ADC.Cen_ADC_nM'; % Species name in the model where ADC is dosed
Sim_Params.Dar_Dose_target_name = 'Central_ADC.Cen_DummyDAR'; % Species name in the model where DAR is dosed
Sim_Params.Scan_Interval = 42;                      % Time interval between two radiographical scan in simulation, (in days)
Sim_Params.Time_efficacy=0:Sim_Params.Scan_Interval:Sim_Params.StopTime; % response measument days 
Sim_Params.UseParallel_simfcn = true;

%First param must be the total volume of tumor
Sim_Params.StatesToLog = {'TV','PSA', 'Cen_ADC_ugml', 'Cen_tAb_ugml', 'Cen_PL_ngml','TME_p_PLLysoFree_nM'}; 
Sim_Params.N_tmr = 1;
Sim_Params.Param_names_for_intrapatient_heterogeneity = {'Tumor_Dia','H_Score','K1','Kmax_p'}; % Parameters that will change across within each patient
Sim_Params.Intrapatient_Heterogeneity = 0.0;       % Variability in parameters across different tumors within patient
Sim_Params.dropout_Switch = true;                  % Switch to turn off/on the dropouts due to New metastats or Non-specific reasons
plot_data_on_waterfall = true;

model_params = { 
    % all the parameters you want to change from what was used in STEP_2
    %param_name          %mu    %sig
    'W'                  70     0;
    'Switch_TumorComp'   1      0;
    'Switch_Tox'         1      0;
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SECTION END %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Params = update_params_from_estimates(Sim_Params, model_params, model);

% Configure SimBiology Model
csObj = getconfigset(model);
csObj.SolverType = 'sundials';
csObj.SolverOptions.RelativeTolerance = 1e-4;
csObj.SolverOptions.AbsoluteTolerance = 1e-4;
csObj.SolverOptions.MaxStep = 1;
csObj.SolverOptions.AbsoluteToleranceScaling = false;
csObj.StopTime = Sim_Params.StopTime;
csObj.TimeUnits = 'day';
csObj.SolverOptions.OutputTimes = Sim_Params.Time_Sim;
csObj.RuntimeOptions.StatesToLog = Sim_Params.StatesToLog;

% Create simfunction
param_names = fieldnames(Sim_Params.param);
observables_names = Sim_Params.StatesToLog;
num_dose_days = length(Sim_Params.Doses_days_in_a_cycle);
doses_targets = [repmat({Sim_Params.Dose_target_name},[1,num_dose_days]),...
                 repmat({Sim_Params.Dar_Dose_target_name},[1,num_dose_days])];
F_simfcn = createSimFunction(model,param_names,observables_names,doses_targets,'UseParallel',Sim_Params.UseParallel_simfcn);




tumor_waterfall_for_all_seed  = nan(length(seed),Sim_Params.N_VP);
psa_waterfall_for_all_seed  = nan(length(seed),Sim_Params.N_VP);
ttr_for_all_seed = nan(length(seed),Sim_Params.N_VP);
orr_for_all_seed = nan(length(seed),1);
psa50_for_all_seed = nan(length(seed),1);

mean_spider = nan(9, length(seed)+1); % 7 is the max no. of visits possible with the current stoptime and scan interval
mean_spider = array2table(mean_spider, 'VariableNames',['Unique_visits', string(1:length(seed))]);
mean_spider.Unique_visits = [0:8]';
mean_psa_spider = nan(9, length(seed)+1); % 7 is the max no. of visits possible with the current stoptime and scan interval
mean_psa_spider = array2table(mean_psa_spider, 'VariableNames',['Unique_visits', string(1:length(seed))]);
mean_psa_spider.Unique_visits = [0:8]';

min_spider = nan(9, length(seed)+1); % 7 is the max no. of visits possible with the current stoptime and scan interval
min_spider = array2table(min_spider, 'VariableNames',['Unique_visits', string(1:length(seed))]);
min_spider.Unique_visits = [0:8]';
min_psa_spider = nan(9, length(seed)+1); % 7 is the max no. of visits possible with the current stoptime and scan interval
min_psa_spider = array2table(min_psa_spider, 'VariableNames',['Unique_visits', string(1:length(seed))]);
min_psa_spider.Unique_visits = [0:8]';

max_spider = nan(9, length(seed)+1); % 7 is the max no. of visits possible with the current stoptime and scan interval
max_spider = array2table(max_spider, 'VariableNames',['Unique_visits', string(1:length(seed))]);
max_spider.Unique_visits = [0:8]';
max_psa_spider = nan(9, length(seed)+1); % 7 is the max no. of visits possible with the current stoptime and scan interval
max_psa_spider = array2table(max_psa_spider, 'VariableNames',['Unique_visits', string(1:length(seed))]);
max_psa_spider.Unique_visits = [0:8]';

pfs_for_all_seed = cell(length(seed));
Tumor_response_alldose = nan([Sim_Params.N_VP*length(Sim_Params.Dose),2]);
PSA_response_alldose = [];

Dose = Sim_Params.Dose;
for d = 1:length(Dose)
    
    for s=1:length(seed)
        Sim_Params.Dose_index = d;
        tic;
        rng(seed(s), 'twister');
        % fprintf('Simulating for seed number %d \n',s);
        fprintf('Simulating for dose = %d and seed = %d \n',Dose(d),seed(s));
        Simulation_Results = run_population_simulation(Sim_Params, F_simfcn, true);
        
        % Storing results for plotting 
        % disp(sum(Simulation_Results.AllPlots.DropoutSummary.ReasonCode == 8)*100/Sim_Params.N_VP); % 8 is the reason code for Bone_metastasis
        [~, TTR] = calculate_TTR(Simulation_Results.AllPlots.Sim_SLD(:,[1, 2, 4]));
        ttr_for_all_seed(s,:) = sort(TTR);
        orr_for_all_seed(s) = Simulation_Results.AllPlots.ORR;
        psa50_for_all_seed(s) = Simulation_Results.AllPlots.PSA50;
    
        [Mean_spider_simulated, Min_spider_simulated, Max_spider_simulated] = calculate_mean_resp_at_clinical_visits(Simulation_Results.AllPlots.Sim_SLD(:,[1,2,4,12]));
        [Mean_psa_spider_simulated, Min_psa_spider_simulated, Max_psa_spider_simulated] = calculate_mean_resp_at_clinical_visits(Simulation_Results.AllPlots.Sim_SLD(:,[1,2,6,12]));
       
        % the following is to resolve an error regarding dimension mismatch when no patients make it till 7th visit
        padded_col = nan(9,1); padded_col(1:height(Mean_spider_simulated)) = Mean_spider_simulated.Mean_Resp_Sim;
        mean_spider.(s+1) =  padded_col;
        padded_col = nan(9,1); padded_col(1:height(Mean_psa_spider_simulated)) = Mean_psa_spider_simulated.Mean_Resp_Sim;
        mean_psa_spider.(s+1) =  padded_col;
    
        padded_col = nan(9,1); padded_col(1:height(Min_spider_simulated)) = Min_spider_simulated.Min_Resp_Sim;
        min_spider.(s+1) =  padded_col;
        padded_col = nan(9,1); padded_col(1:height(Min_psa_spider_simulated)) = Min_psa_spider_simulated.Min_Resp_Sim;
        min_psa_spider.(s+1) =  padded_col;
    
        padded_col = nan(9,1); padded_col(1:height(Max_spider_simulated)) = Max_spider_simulated.Max_Resp_Sim;
        max_spider.(s+1) =  padded_col;
        padded_col = nan(9,1); padded_col(1:height(Max_psa_spider_simulated)) = Max_psa_spider_simulated.Max_Resp_Sim;
        max_psa_spider.(s+1) =  padded_col;
    
        sim_tumor_waterfall = sortrows(Simulation_Results.AllPlots.Sim_Waterfall,-2);
        tumor_waterfall_for_all_seed(s,:) = sim_tumor_waterfall(:,2)';
        psawaterfall = sim_tumor_waterfall(:,4);
        % psawaterfall = sim_tumor_waterfall(:,4);%sim_tumor_waterfall(sim_tumor_waterfall(:,2)<0,4);
        psalen = height(psawaterfall);
        psawaterfall = sort(psawaterfall, 'descend');
        psa_waterfall_for_all_seed(s,:) = sort(psawaterfall, 'descend')';
        pfs_for_all_seed{s} = Simulation_Results.AllPlots.Sim_PFS;
        
        % PLOTTING
        % plot_waterfall(Simulation_Results, 'Data', Lesion_Waterfall, PSA_Waterfall); 
        % plot_pfs(Simulation_Results);                             % [OPTIONAL] Lines marked as optional can be commented without any issue
        % plot_swimmer(Simulation_Results);                         % [OPTIONAL]
        % plot_spider(Simulation_Results, Sim_Params, median(Lesion_Waterfall.TTR, 'omitnan'));  
        % plot_psa_spider(Simulation_Results);                            % [OPTIONAL]
        % plot_observables(Simulation_Results, Sim_Params)          % [OPTIONAL] 
        % summarize_patient_response(Simulation_Results);           % [OPTIONAL]
        drawnow;
    
    
        % % For creating dose reponse plots
        Tumor_response_alldose(Sim_Params.N_VP*(d-1)+1:Sim_Params.N_VP*d, 2) = sim_tumor_waterfall(:,2);
        Tumor_response_alldose(Sim_Params.N_VP*(d-1)+1:Sim_Params.N_VP*d, 1) = repmat(Dose(d),[Sim_Params.N_VP,1]);

        PSA_response_alldose(end+1:end+psalen,2) = min(100,psawaterfall);
        PSA_response_alldose(end-psalen+1:end,1) = repmat(Dose(d),[psalen,1]);
        toc;
    end

    % plot_mean_spider_for_all_seed(Mean_spider_clinical, mean_spider, min_spider, max_spider, median(Lesion_Waterfall.TTR, 'omitnan'), ttr_for_all_seed);
    % %%
    % plot_mean_psa_spider_for_all_seed(Median_PSA_spider_clinical, mean_psa_spider, min_psa_spider, max_psa_spider);
    %%
    % Plot_waterfall_for_all_seed(psa_waterfall_for_all_seed, tumor_waterfall_for_all_seed, plot_data_on_waterfall, Lesion_Waterfall, PSA_Waterfall, orr_for_all_seed, psa50_for_all_seed, Sim_Params);
    %%
    % Plot_waterfall_for_all_seed1(psa_waterfall_for_all_seed, tumor_waterfall_for_all_seed, plot_data_on_waterfall, Lesion_Waterfall, PSA_Waterfall);
    % Plot_PFS_for_all_seed(Sim_Params, Dose(d), pfs_for_all_seed, true);
    % drawnow;
    % analyze_binomial_p(Sim_Params.N_VP,psa50_for_all_seed);

end


Plot_dose_response(Tumor_response_alldose);
Plot_dose_response(PSA_response_alldose);



























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTION DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
function Simulation_Results = run_population_simulation(Sim_Params, F_simfcn, verbose)
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
        Dose_nanomoles = (Sim_Params.Dose(Sim_Params.Dose_index) * vp_init_Array(1).(Sim_Params.parameter_name_for_bodyweight) / (1000 * vp_init_Array(1).(Sim_Params.parameter_name_for_mol_weight_ADC))) * 10^9;  % Convert mpk to nanomole
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
    [ORR, PSA50, Sim_Waterfall] = calculate_bor(Sim_SLD, PatientList);

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
    Simulated_Data.ORR = ORR;
    Simulated_Data.PSA50 = PSA50;
    % Bundle Final Results
    Simulation_Results = struct();
    Simulation_Results.AllPlots = Simulated_Data;
    Simulation_Results.Dose = Sim_Params.Dose(Sim_Params.Dose_index);
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
        p_sig = params.(p_name).sig*Sim_Params.variability;
        
        % Generate the patient-specific base value (Log-Normal Distribution)
        if strcmp(p_name,Sim_Params.parameter_name_for_hscore)
            patient_base_vals.(p_name) = betarnd(5, 2.5, [1, 1]) * 300*Sim_Params.variability + (1-Sim_Params.variability)*p_mu;
        else
            patient_base_vals.(p_name) = p_mu * exp(randn(1,1) * p_sig);
        end
    end
    vp_init_Array = repmat(cell2struct(cell(numel(all_param_names)+2, 1), [{'PatientID', 'TumorID'}, all_param_names]), 1, n_tumors);
    % Loop through Tumors to generate the Struct Array
    for t_id = 1:n_tumors
        
        % Assign Metadata Fields First (to keep them at the left/top of struct)
        vp_init_Array(t_id).PatientID = NaN;
        vp_init_Array(t_id).TumorID = t_id;
        
        % Loop through model parameters
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
    BLDrop_Params = Sim_Params.BLDrop_Params;
    NonSpecific_dropouts = Sim_Params.NonSpecific_dropouts;
    dropout_Switch = Sim_Params.dropout_Switch;

    Sim_SLD = zeros(size(SLD, 1), 13); 
    row_pointer = 1;
    PatientList = 1:N_VP; % Get list of all patient IDs
    for p_id = 1:N_VP
        patient_id = PatientList(p_id);
        baseline_BL_count = round(3.4928 * exp(randn(1,1)*1.9736));
        lesion_level_for_drop_out = IndTS(IndTS(:,1) == patient_id, :);
        current_patient_data = SLD(SLD(:,1) == patient_id, :);
        n_rows = size(current_patient_data, 1);        
        patient_level_for_drop_out = [current_patient_data, zeros(n_rows, 6)];
        for time_point = 2:length(Time_efficacy)
            % Target Lesion PD
            patient_level_for_drop_out(time_point,7) = double(patient_level_for_drop_out(time_point,4) > 20);
            
            % Non-Target Lesion PD
            Deltas = lesion_level_for_drop_out(lesion_level_for_drop_out(:,2) == Time_efficacy(time_point), 5);
            AbsChanges = lesion_level_for_drop_out(lesion_level_for_drop_out(:,2) == Time_efficacy(time_point), 4) - ...
                lesion_level_for_drop_out(lesion_level_for_drop_out(:,2) == 0, 4);
            patient_level_for_drop_out(time_point,8) = any(Deltas >= 40 & AbsChanges > 5);
            
            % New Metastases
            diam = patient_level_for_drop_out(time_point,3);
            
            if dropout_Switch
                
                probNM = NM_Params(1) + (NM_Params(2) * diam^NM_Params(4)) / ((NM_Params(3)^NM_Params(4)) + (diam^NM_Params(4)));
                patient_level_for_drop_out(time_point,9) = rand < probNM;

                
                if time_point ~= 2
                    probBLDrop_SLD1 = (BLDrop_Params(1) * diam) / (BLDrop_Params(2) + diam);
                    probBLDrop_NBL = BLDrop_Params(3) * baseline_BL_count / (BLDrop_Params(4) + baseline_BL_count);
                    diam_prev = patient_level_for_drop_out(time_point-1,3);
                    probBLDrop_SLD2 = (BLDrop_Params(1) * diam_prev) / (BLDrop_Params(2) + diam_prev);
                    patient_level_for_drop_out(time_point,10) = rand < (probBLDrop_SLD1 + probBLDrop_NBL)*(probBLDrop_SLD2 + probBLDrop_NBL);
                end
                
                % Nonspecific dropouts
                if rand < NonSpecific_dropouts
                    patient_level_for_drop_out(time_point,11) = 1;
                else
                    patient_level_for_drop_out(time_point,11) = 0;
                end
            else
                patient_level_for_drop_out(time_point,9:11) = 0;
            end
            
            % Total PD flags
            patient_level_for_drop_out(time_point, 12) = sum(patient_level_for_drop_out(time_point,7:11));
            patient_level_for_drop_out(time_point, 13) = baseline_BL_count; 
        end
        % Determine how many rows we are adding
        n_rows = size(patient_level_for_drop_out, 1);
        
        % Insert into the preallocated slots
        Sim_SLD(row_pointer : row_pointer + n_rows - 1, :) = patient_level_for_drop_out;
        
        % Update the pointer for the next iteration
        row_pointer = row_pointer + n_rows;
    end
end
 
function [ORR, PSA50, Sim_Waterfall] = calculate_bor(Sim_SLD, PatientList)
    %
    % Calculates the Best Overall Response (BOR) for the waterfall plot
    % INPUT : 
    %       Sim_SLD : The array containing complete Simulated SLD results with dropout informations
    %       PatientList : A simple column vector = 1:N_VP  
    %
    % OUTPUT :
    %       Sim_Waterfall : Array containing BOR for each patient along with the corresponding timepoint
    %       ORR : Objective Overall Response in % as per RECIST 1.1 Criteria

    nPatients = numel(PatientList);
    Sim_Waterfall = zeros(nPatients, 5);
    Sim_Waterfall(:,1) = PatientList;
    n_ORR = 0;
    n_PSA50 = 0;
    for patient_id = 1:nPatients
        subj = Sim_SLD(Sim_SLD(:,1) == PatientList(patient_id), :);
        drop_idx = find(subj(:,12) >= 1, 1);
        is_progressive = evaluate_psa_progression(subj); % PSA progression flag
        
        if isempty(drop_idx)
            best_change_TV = min(subj(subj(:,2) ~= 0, 4));
            best_time_TV = subj(subj(:,4) == best_change_TV, 2);
            best_change_PSA = min(subj(subj(:,2) ~= 0, 6));
            best_time_PSA = subj(subj(:,6) == best_change_PSA, 2);
            best_change_ORR = best_change_TV;
        else

            % for waterfall calculation
            till_drop_subj = subj(1:drop_idx, :);
            best_change_TV = min(till_drop_subj(till_drop_subj(:,2) ~= 0, 4));
            best_time_TV = till_drop_subj(till_drop_subj(:,4) == best_change_TV, 2);
            best_change_PSA = min(till_drop_subj(till_drop_subj(:,2) ~= 0, 6));
            best_time_PSA = till_drop_subj(till_drop_subj(:,6) == best_change_PSA, 2);
            
            % for calculation ORR numerator
            reason_row = subj(drop_idx,:);
            non_spec_drp_flag = (reason_row(:,11)==1) && (reason_row(:,12)==1); % (Non-specific dropout is true) && (only one dropout criteria is satisfied)
            if non_spec_drp_flag
                best_change_ORR = best_change_TV;
            else
                if drop_idx == 2
                    best_change_ORR = 0;
                else
                    pre_drop_subj = subj(1:drop_idx-1, :);
                    best_change_ORR = min(pre_drop_subj(pre_drop_subj(:,2) ~= 0, 4));
                end
            end
        end
        if isempty(best_time_TV)
            best_time_TV = 0;
        end

        Sim_Waterfall(patient_id,2) = best_change_TV;
        Sim_Waterfall(patient_id,3) = best_time_TV(1);

        Sim_Waterfall(patient_id,4) = best_change_PSA;
        if ~isempty(best_time_PSA)
            Sim_Waterfall(patient_id,5) = best_time_PSA(1);
        else
            Sim_Waterfall(patient_id,5) = 0;
        end

        if (best_change_ORR <= -30) && ~is_progressive
            n_ORR = n_ORR + 1;
        end
        if evaluate_psa50_response(subj) && (best_change_ORR < 0) % PSA50 response flag
            n_PSA50 = n_PSA50 + 1;
        end
    end
    ORR = n_ORR/nPatients;
    PSA50 = n_PSA50/nPatients; % count of patients
end
 
function is_responder = evaluate_psa50_response(subj)
    % EVALUATE_PSA50_RESPONSE Determines if an mCRPC patient achieves a 
    % confirmed PSA50 response based on PCWG3 criteria.
    %
    % Input: 
    %   subj - Matrix where:
    %          Col 2: Time (days)
    %          Col 5: PSA concentration (ng/mL)
    %          Col 12: Dropout flag (> 0 indicates dropout)
    % Output:
    %   is_responder - Boolean flag (true/false)

    % Default output
    is_responder = false;
    
    dropout_idx = find(subj(:, 12) > 0, 1, 'first');
    if ~isempty(dropout_idx)
        subj = subj(1:dropout_idx, :);
    end
    
    num_rows = size(subj, 1);
    time_days = subj(:, 2);
    psa = subj(:, 5);
    baseline_psa = psa(1);
    
    % Evaluate PCWG3 PSA50 Criteria
    for i = 2:num_rows
        % Check if the current point shows a >= 50% drop from baseline
        has_initial_drop = (psa(i) <= 0.5 * baseline_psa);
        
        if has_initial_drop
            % Look ahead for a consecutive confirmatory time point
            for k = (i+1):num_rows
                
                % If any intermediate point rises back above the -50% threshold, the streak is broken.
                if psa(k) > 0.5 * baseline_psa
                    break; 
                end
                
                time_diff = time_days(k) - time_days(i);
                % Confirmatory test must be >= 3 weeks (21 days) apart
                if time_diff >= 21
                    is_responder = true;
                    return; % Confirmed PSA50 responder
                end
            end
        end
    end
end
   
function is_progressive = evaluate_psa_progression(subj)
    % EVALUATE_PSA_PROGRESSION Determines if an mCRPC patient has progressed
    % based on PCWG3 PSA criteria.
    %
    % Input: 
    %   subj - Matrix where:
    %          Col 2: Time (days)
    %          Col 5: PSA concentration (ng/mL)
    %          Col 12: Dropout flag (> 0 indicates dropout)
    % Output:
    %   is_progressive - Boolean flag (true/false)

    is_progressive = false; % Default output
    
    dropout_idx = find(subj(:, 12) > 0, 1, 'first');
    if ~isempty(dropout_idx)
        subj = subj(1:dropout_idx, :);
    end
    
    time_days = subj(:, 2);
    psa = subj(:, 5);
    
    % Evaluate PCWG3 Criteria 
    num_rows = size(subj, 1);
    for i = 2:num_rows
        nadir = min(psa(1:i));
        
        % Check if current point satisfies initial progression thresholds:
        % 1. Absolute increase of >= 2 ng/mL from nadir
        % 2. Relative increase of >= 25% from nadir
        has_initial_rise = (psa(i) >= nadir + 2) && (psa(i) >= 1.25 * nadir);
        
        if has_initial_rise
            % Look ahead for a consecutive confirmatory time point
            for k = (i+1):num_rows
                time_diff = time_days(k) - time_days(i);
                
                % Confirmatory test must be >= 3 weeks (21 days) later
                if time_diff < 21
                    continue; 
                end                
                
                % If the initial rise happens before Day 84 (12 weeks), 
                % confirmation must occur at or after Day 84.
                if (time_days(i) < 84) && (time_days(k) < 84)
                    continue;
                end
                
                % The confirmation point must also sustain the threshold
                is_sustained = (psa(k) >= nadir + 2) && (psa(k) >= 1.25 * nadir);
                
                if is_sustained
                    is_progressive = true;
                    return; % Progression confirmed; exit function early
                end
            end
        end
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

    Sim_PFS = [Sim_Waterfall(:,1), zeros(size(Sim_Waterfall,1), 2)];
    for i = 1:size(Sim_PFS,1)
        subj = Sim_SLD(Sim_SLD(:,1) == Sim_PFS(i,1), :);
        pd_idx = find(subj(:,12) >= 1, 1);
        
        if isempty(pd_idx)
            Sim_PFS(i,2) = subj(end,2);  % Last time point
            Sim_PFS(i,3) = 0;            % Censored
        else
            Sim_PFS(i,2) = subj(pd_idx,2);  % PD time
            Sim_PFS(i,3) = 1;              % Event occurred
        end
    end
    
    Sim_PFS = sortrows(Sim_PFS, 2);
    % Sim_PFS(:,4) = size(Sim_PFS,1):-1:1;
    % Sim_PFS(:,5) = 1 - Sim_PFS(:,3) ./ Sim_PFS(:,4);
    % Sim_PFS(1,6) = 1;
    % for i = 2:size(Sim_PFS,1)
    %     Sim_PFS(i,6) = Sim_PFS(i-1,6) * Sim_PFS(i,5);
    % end
    Sim_PFS = [[0 0 0]; Sim_PFS];  % Add (0,1) starting point

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
        idx_dropout = find(data_i(:,12) >= 1, 1, 'first');
        
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
    
       % ReasonCode: 5=TargetPD, 6=NonTargetPD, 7=NewMet, 8=Bone mets, 9=Nonspecific, 0=No dropout
    for i = 1:nPatients
        pid = PatientList(i);
        data_i = Sim_SLD(Sim_SLD(:,1) == pid,:);
        idx_dropout = find(data_i(:,12) >= 1, 1, 'first');
        
        if ~isempty(idx_dropout)
            DropoutSummary.DropoutTime(i) = data_i(idx_dropout,2);
            reasons = data_i(idx_dropout,7:11);
            if reasons(4)==1
                DropoutSummary.ReasonCode(i) = 8; % Bone metastases
            elseif reasons(1)==1
                DropoutSummary.ReasonCode(i) = 5; % Target lesion PD
            elseif reasons(2)==1
                DropoutSummary.ReasonCode(i) = 6; % Non-target lesion PD
            elseif reasons(3)==1
                DropoutSummary.ReasonCode(i) = 7; % New metastases
            elseif reasons(5)==1
                DropoutSummary.ReasonCode(i) = 9; % Nonspecific dropout
            end
        else
            DropoutSummary.DropoutTime(i) = max(data_i(:,2)); % completed
            DropoutSummary.ReasonCode(i) = 0; % no dropout
        end
    end
end
 
function plot_waterfall(Simulation_Results, plot_mode, Lesion_Waterfall, PSA_Waterfall)
    %
    % plot_waterfall: Generates a pooled waterfall plot from simulation results.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %       plot_mode :         'Data' or 'Sim' Choice between plotting a comparison with the data, or plot only the simulation as bars
    %       Waterfall:          Processed data containing information about the waterfall plot from the raw data
    %

    gem_hex_codes = rgb2hex(orderedcolors("gem"));

    simulated_waterfall_Data = Simulation_Results.AllPlots.Sim_Waterfall;

    % Sort by Best Overall Response (BOR)
    simulated_waterfall_Data = sortrows(simulated_waterfall_Data, -2);
    N_VP = size(simulated_waterfall_Data,1);   
    sim_waterfall = simulated_waterfall_Data(:,2);    
    if strcmp(plot_mode,'Data')
        psawaterfall = min(100,simulated_waterfall_Data(:,4));%end-height(PSA_Waterfall)+1:end,4));
        psawaterfall = sort(psawaterfall, 'descend');
        figure;
        hold on;    
        % Create the bar plot for the Data
        data_prc_change=sort(Lesion_Waterfall.BOR,'descend');
        if height(data_prc_change) == N_VP
            bar(data_prc_change, 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
        end

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
        fontname('Poppins');
        hold off;
    
        figure;
        hold on;    
        % Create the bar plot for the Data
        data_psa_change=sort(PSA_Waterfall.BOR,'descend');
        if height(psawaterfall) == height(data_psa_change)
            bar(data_psa_change, 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
        end

        % Create a curve for simulated waterfall plot
        plot(1:height(psawaterfall), psawaterfall, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0,'DisplayName', 'Simulated');
    
        % yline(20, '--k','HandleVisibility','off'); yline(-30, '--k','HandleVisibility','off');
        
        % Formatting
        xlabel('Virtual Patients (Sorted by Response)','FontWeight','bold');
        ylabel('BOR (%)','FontWeight','bold');
        global_min = min([min(psawaterfall(:)), min(data_psa_change(:))]);
        global_max = max([max(psawaterfall(:)), max(data_psa_change(:))]);
        
        y_min = max(-105, global_min - 10); 
        y_max = global_max + 10;
        ylim([y_min, y_max]);
        xlim([0.5,height(psawaterfall)+0.5]);
        grid off; box on; set(gca,'LineWidth',1.25);
        pbaspect([3 2 1]);
        fontsize(15,'points');
        legend('Location','southwest','Box','on','FontSize',15);
        title('Calibrating PSA Waterfall Plot','FontWeight','bold','FontSize',18);
        fontname('Poppins');
        hold off;
    else
        psawaterfall = min(100,simulated_waterfall_Data(:,4));%min(100,simulated_waterfall_Data(simulated_waterfall_Data(:,2)<0,4));
        psawaterfall = sort(psawaterfall, 'descend');
        figure;
        hold on;    
        % Create a curve for simulated waterfall plot
        bar(sim_waterfall, 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Sim');
    
        % Add lines to define Stable Disease zone (SD)
        yline(20, '--k','HandleVisibility','off'); yline(-30, '--k','HandleVisibility','off');
        
        % Formatting
        xlabel('Virtual Patients (Sorted by Response)','FontWeight','bold');
        ylabel('BOR (%)','FontWeight','bold');
       
        ylim([-105, 105]);
        xlim([0.5,N_VP+0.5]);
        grid off; box on; set(gca,'LineWidth',1.25);
        pbaspect([3 2 1]);
        fontsize(15,'points');
        legend('Location','southwest','Box','on','FontSize',15);
        title('Simulated Tumor Waterfall Plot','FontWeight','bold','FontSize',18);
        fontname('Poppins');
        hold off;
    
        figure;
        hold on;    
        
        % Create a curve for simulated waterfall plot
        bar(psawaterfall, 'FaceColor', gem_hex_codes(1), 'BarWidth', 0.6, 'EdgeColor',gem_hex_codes(1), 'DisplayName','PSA Sim');
        
        % Formatting
        xlabel('Virtual Patients (Sorted by Response)','FontWeight','bold');
        ylabel('BOR (%)','FontWeight','bold');
        ylim([-105, 105]);
        xlim([0.5,height(psawaterfall)+0.5]);
        grid off; box on; set(gca,'LineWidth',1.25);
        pbaspect([3 2 1]);
        fontsize(15,'points');
        legend('Location','southwest','Box','on','FontSize',15);
        title('Simulated PSA Waterfall Plot','FontWeight','bold','FontSize',18);
        fontname('Poppins');
        hold off
    end
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
    fontsize(15,'points'); fontname('Poppins');
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
    reasonColors = containers.Map({0, 5, 6, 7, 8, 9}, {
        [0.3 0.7 0.3];  % Continue = green
        [0.8 0.5 0.0];  % Target PD = orange
        [0.0 0.4 0.9];  % Non-target PD = blue
        [0.9 0.2 0.2];  % New metastases = red
        [0.4 0.4 0.4];  % Bone  mets = gray
        [1 1 1]   % Nonspecific = gray
        });
    reasonLabels = {'Continued', 'Target PD', 'Non-target PD', 'New Metastases', 'Bone Metastases', 'Nonspecific'};
    reasonCodes  = [0, 5, 6, 7, 8, 9];

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
        barh(order, t_end, 0.8, 'FaceColor', reasonColors(reason), 'EdgeColor', 'none');
    end
    
    % Formatting
    ylim([0, size(DropoutSummary.DropoutTime,1) + 1]);
    xlim([0, 300]);
    xlabel('Time (days)','FontWeight','bold');
    ylabel('Virtual Patients','FontWeight','bold');
    fontsize(15,'points');
    title('Simulated Swimmer Plot','FontWeight','bold','FontSize',16);
    set(gca, 'YDir', 'reverse'); 
    grid off; box on; set(gca, 'LineWidth', 1.25);
    pbaspect([3 2 1]);
    fontname('Poppins');
    
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
 
function plot_spider(Simulation_Results, Sim_Params, median_TTR)
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
    [median_TTR_Sim, ~] = calculate_TTR(Sim_SLD(:,[1, 2, 4])); % this subset of Sim_SLD contains the patient_id, Time, % change in SLD
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
        plot(data_cut(:,2), data_cut(:,4), '-', 'Color', gem_hex_codes(2), 'LineWidth', 1,'Marker','.','MarkerSize',8);
    end

    
    
    xlabel('Time (days)', 'FontWeight','bold');
    ylabel('% Change in SLD', 'FontWeight','bold');
    fontsize(15, 'points');
    text(10, 100 , sprintf('Mean Depth of Response at 1st Clinical Visit : %.2f', mean(response_at_1st_visit, 'omitnan')), 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'left');
    text(10, 90, '-----------------------------------------------------------', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'left');
    text(10, 80 , sprintf('Median TTR (Clinical) : %.1f', median_TTR), 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'left');
    text(10, 70 , sprintf('Median TTR (Simulated) : %.1f', median_TTR_Sim), 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'left');
    
    yline(0,'--k');
    title('Simulated Tumor Trajectories');
    ylim([-110 110]);
    xlim([0, Sim_Params.StopTime + 5]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([ 3 2 1]);
    % Formatting
    % h_leg = zeros(4,1);
    % for k = 1:4
    %     h_leg(k) = plot(nan, nan, '-', 'Color', IHC_colors(k,:), 'LineWidth', 4); % bogus plots to ensure the IHC Colors are shown only once no matter how many curves are there on the spider plot
    % end
    % legend(h_leg, {'IHC0','IHC1+','IHC2+','IHC3+'}, 'Location','northeast');
    fontname('Poppins');
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
        % if min(data_cut(2:end,4))<=0
        %     plot(data_cut(:,2), data_cut(:,6), '-', 'Color', gem_hex_codes(1), 'LineWidth', 2);
        % end
        plot(data_cut(:,2), data_cut(:,6), '-', 'Color', gem_hex_codes(1), 'LineWidth', 1);
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
    xlim([0, 735]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([1 1 1]);
    % legend(h_leg, {'IHC0','IHC1+','IHC2+','IHC3+'}, 'Location','northeast');
    fontname('Poppins');
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
        
        fontname('Poppins');
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
 
function [median_TTR_Sim, TTR] = calculate_TTR(Sim_SLD_subset)
    
    subjects = unique(Sim_SLD_subset(:,1));
    TTR = NaN(length(subjects), 1);
    for i = 1:length(subjects)
        subject_id = subjects(i);
        
        subj_data = Sim_SLD_subset(Sim_SLD_subset(:,1) == subject_id, :);
                
        t = subj_data(:,2); y = subj_data(:,3);
        
        idx = find(y <= -30, 1, 'first');
        
        if ~isempty(idx)
            % Get the coordinates for the points immediately before and after crossing
            t0 = t(idx-1);
            y0 = y(idx-1);
            t1 = t(idx);
            y1 = y(idx);
            
            % Safeguard against division by zero in case of identical y values
            if y1 == y0
                TTR(i) = t1;
            else
                % Linear interpolation 
                TTR(i) = t0 + ((-30 - y0) / (y1 - y0)) * (t1 - t0);
            end
        end
        % If idx is empty (trajectory never crossed -30), the TTR value remains NaN
    end
    median_TTR_Sim = median(TTR, 'omitnan');
end
 
function [Mean_response_table, Min_response_table, Max_response_table] = calculate_mean_resp_at_clinical_visits(Sim_SLD_subset)
    
    num_rows = size(Sim_SLD_subset,1); 

    keep_rows = true(num_rows, 1); % Default to keeping all rows
    has_dropped_out = false;       % Tracker for the current patient
    current_subj = Sim_SLD_subset(1, 1); 
    
    for i = 1:num_rows
        % If we reach a new patient, reset the tracker
        if Sim_SLD_subset(i, 1) ~= current_subj
            current_subj = Sim_SLD_subset(i, 1);
            has_dropped_out = false; 
        end
        
        if has_dropped_out
            % If patient dropped out in a prior visit, discard this row
            keep_rows(i) = false;
        else
            % If patient drops out AT this visit, keep this row, but mark 
            % 'has_dropped_out' as true so subsequent rows are discarded.
            if Sim_SLD_subset(i, 4) > 0
                has_dropped_out = true;
            end
        end
    end
    
    Clean_Sim_SLD = Sim_SLD_subset(keep_rows, :);
    
    % Assign Visit Numbers
    num_clean_rows = size(Clean_Sim_SLD, 1);
    Visit_Num = zeros(num_clean_rows, 1);
    current_subj = Clean_Sim_SLD(1, 1); 
    visit_count = 0; 
    
    for i = 1:num_clean_rows
        if Clean_Sim_SLD(i, 1) ~= current_subj
            % New patient found, reset tracker
            current_subj = Clean_Sim_SLD(i, 1);
            visit_count = 0;
        end
        Visit_Num(i) = visit_count;
        visit_count = visit_count + 1;
    end
    
    Clean_Sim_SLD(:, 5) = Visit_Num; 
    
    unique_visits = unique(Clean_Sim_SLD(:, 5));
    mean_resp_list = NaN(length(unique_visits), 1);
    max_resp_list = NaN(length(unique_visits), 1);
    min_resp_list = NaN(length(unique_visits), 1);
    
    for k = 1:length(unique_visits)
        vis = unique_visits(k);
        responses = min(100, Clean_Sim_SLD(Clean_Sim_SLD(:, 5) == vis, 3),'includemissing'); 
        mean_resp_list(k) = mean(responses, 'omitnan');
        max_resp_list(k) = max(responses);
        min_resp_list(k) = min(responses);
        
    end
    
    Mean_response_table = table(unique_visits, mean_resp_list, 'VariableNames', {'Sim_visits', 'Mean_Resp_Sim'});
    Max_response_table = table(unique_visits, max_resp_list, 'VariableNames', {'Sim_visits', 'Max_Resp_Sim'});
    Min_response_table = table(unique_visits, min_resp_list, 'VariableNames', {'Sim_visits', 'Min_Resp_Sim'});
end
  
function plot_mean_spider_for_all_seed(Clinical_table, Simulated_table, min_table, max_table, median_TTR, ttr_for_all_seed)
    
    trial_medians = median(ttr_for_all_seed, 2, 'omitnan');
    overall_median_ttr = median(trial_medians, 'omitnan');
    ci_95_ttr = prctile(trial_medians, [25, 75]);
    lower_ci = ci_95_ttr(1);
    upper_ci = ci_95_ttr(2);

    sim_data_matrix = table2array(Simulated_table(:, 2:end));
    min_data_matrix = table2array(min_table(:, 2:end));
    max_data_matrix = table2array(max_table(:, 2:end));
    mean_sim_across_seeds = mean(sim_data_matrix, 2, 'omitnan'); 
    upper_bound = max(max_data_matrix')';
    lower_bound = min(min_data_matrix')';

    figure; hold on;
    box on; grid off; set(gca, 'LineWidth', 1.25);
    x_fill = [Simulated_table.Unique_visits; flipud(Simulated_table.Unique_visits)];
    y_fill = [upper_bound; flipud(lower_bound)];
    fill(x_fill, y_fill, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Simulated 95% CI');
    
    errorbar(Clinical_table.Clinical_visits, Clinical_table.Mean_Resp_Clinical, ...
        Clinical_table.Mean_Resp_Clinical - Clinical_table.Min_Resp_Clinical, ...
        Clinical_table.Max_Resp_Clinical - Clinical_table.Mean_Resp_Clinical, ...
        'o','Color','k','LineWidth',2,'MarkerFaceColor','k','DisplayName','Clinical(Range)');
    % plot(Clinical_table.Clinical_visits, Clinical_table.Mean_Resp_Clinical, 'DisplayName','Clinical');
    plot(Simulated_table.Unique_visits, mean_sim_across_seeds, '-', 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName','Simulated (Mean)');
    xlabel('Patient Visit No.');
    ylabel('Mean response in SLD (%)');
    legend('Location','southeast','Box','on');
    text(7, 95 , sprintf('Median TTR (Clinical) : %.1f', median_TTR), 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');
    text(7, 85 , sprintf('Median TTR[IQR] (Simulated) : %.1f[%.1f-%.1f]', overall_median_ttr,lower_ci, upper_ci), 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');
    
    ylim([-105, 105]); xlim([-0.1, 7.2])
    fontsize(15,'points');
    title('Comparing Mean response for clinical and simulated data','FontWeight','bold','FontSize',18);
    fontname('Poppins');
    hold off;
end
   
function plot_mean_psa_spider_for_all_seed(Clinical_table, Simulated_table,  min_table, max_table)
    
    sim_data_matrix = table2array(Simulated_table(:, 2:end));
    min_data_matrix = table2array(min_table(:, 2:end));
    max_data_matrix = table2array(max_table(:, 2:end));
    mean_sim_across_seeds = mean(sim_data_matrix, 2, 'omitnan'); 
    upper_bound = max(max_data_matrix')';
    lower_bound = min(min_data_matrix')';

    figure; hold on;
    box on; grid off; set(gca, 'LineWidth', 1.25);
    x_fill = [Simulated_table.Unique_visits; flipud(Simulated_table.Unique_visits)];
    y_fill = [upper_bound; flipud(lower_bound)];
    fill(x_fill, y_fill, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Simulated 95% CI');
    
    errorbar(Clinical_table.Clinical_visits, Clinical_table.Mean_Resp_Clinical, ...
        Clinical_table.Mean_Resp_Clinical - Clinical_table.Min_Resp_Clinical, ...
        Clinical_table.Max_Resp_Clinical - Clinical_table.Mean_Resp_Clinical, ...
        'o','Color','k','LineWidth',2,'MarkerFaceColor','k','DisplayName','Clinical(Range)');
    % plot(Clinical_table.Clinical_visits, Clinical_table.Mean_Resp_Clinical, 'DisplayName','Clinical');
    plot(Simulated_table.Unique_visits, mean_sim_across_seeds, '-', 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName','Simulated (Mean)');
    xlabel('Patient Visit No.');
    ylabel('Mean response in SLD (%)');
    legend('Location','southeast','Box','on');

    ylim([-105, 105]); xlim([-0.1, 7.2])
    fontsize(15,'points');
    title('Comparing Mean PSA response for clinical and simulated data','FontWeight','bold','FontSize',17);
    fontname('Poppins');
    hold off;
end
 
function Plot_waterfall_for_all_seed1(psa_waterfall_for_all_seed, tumor_waterfall_for_all_seed, plot_data, Lesion_Waterfall, PSA_Waterfall)
    figure; 
    hold on;
    box on; set(gca,'LineWidth', 1.25);

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    if plot_data
        bar(sort(Lesion_Waterfall.BOR,'descend'), 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    end

    meanVals = mean(tumor_waterfall_for_all_seed, 1,'omitmissing');
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
    title('PD Calibration',FontSize = 18,FontWeight="bold");
    legend('Box','off', 'Location', 'southwest','FontSize',13);
    % ylim([-100,100]);
    % xlim([0.5, Sim_Params.N_VP+0.5]);

    pbaspect([3 2 1]);
    fontname('Poppins');


    % PSA Waterfall plotting
    figure; 
    hold on;
    box on; set(gca,'LineWidth', 1.25);

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    if plot_data
        bar(sort(PSA_Waterfall.BOR,'descend'), 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
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
    title('PSA Calibration',FontSize = 18,FontWeight="bold");
    legend('Box','off', 'Location', 'southwest','FontSize',13);
    ylim([-100,100]);
    % xlim([0.5, Sim_Params.N_VP+0.5]);

    pbaspect([3 2 1]);
    fontname('Poppins');
end
 
function Plot_waterfall_for_all_seed(psa_waterfall_for_all_seed, tumor_waterfall_for_all_seed, plot_data, Lesion_Waterfall, PSA_Waterfall, orr_for_all_seed, psa50_for_all_seed, Sim_Params)
    figure; 
    hold on;
    box on; set(gca,'LineWidth', 1.25);

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    if plot_data && height(Lesion_Waterfall) == Sim_Params.N_VP
        bar(sort(Lesion_Waterfall.BOR,'descend'), 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    end

    medianVals = median(tumor_waterfall_for_all_seed, 1, 'omitnan');
    prc5 = prctile(tumor_waterfall_for_all_seed, 5, 1); 
    prc95 = prctile(tumor_waterfall_for_all_seed, 95, 1);
    n=length(medianVals);
    medianCI_upper= nan([1,n]);
    medianCI_lower= nan([1,n]);
    prc5CI_upper= nan([1,n]);
    prc5CI_lower= nan([1,n]);
    prc95CI_upper= nan([1,n]);
    prc95CI_lower= nan([1,n]);

    for i = 1:length(medianVals)
        CI95=boot95CI(tumor_waterfall_for_all_seed(:,i),0.5);
        medianCI_lower(i)=CI95(1);
        medianCI_upper(i)= CI95(2);

        CI95=boot95CI(tumor_waterfall_for_all_seed(:,i),0.05);
        prc5CI_lower(i)=CI95(1);
        prc5CI_upper(i)= CI95(2);

        CI95=boot95CI(tumor_waterfall_for_all_seed(:,i),0.95);
        prc95CI_lower(i)=CI95(1);
        prc95CI_upper(i)= CI95(2);
    end
    x = 1:numel(medianVals);

    % Shaded region (min–max)
    fill([x fliplr(x)], [prc5CI_lower fliplr(prc5CI_upper)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(5), ...
        'HandleVisibility','off');
    fill([x fliplr(x)], [prc95CI_lower fliplr(prc95CI_upper)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(3), ...
        'HandleVisibility','off');
    fill([x fliplr(x)], [medianCI_lower fliplr(medianCI_upper)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(7), ...
        'HandleVisibility','off');

    % Overlay lines
    plot(x, medianVals, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Median Simulated Waterfall');
    plot(x, prc95, '-', 'Color', gem_hex_codes(3), 'LineWidth', 2.0, 'DisplayName', '95th percentile');
    plot(x, prc5, '-', 'Color', gem_hex_codes(5), 'LineWidth', 2.0, 'DisplayName', '5th percentile');
    yline(20, '--k','HandleVisibility','off'); yline(-30, '--k','HandleVisibility','off');
    yline(-30, '--k','HandleVisibility','off'); yline(-30, '--k','HandleVisibility','off');

    median_orr = median(orr_for_all_seed);
    ORR_CI_95 = prctile(orr_for_all_seed, [25, 75]);

    % Mention ORR on the plot
    text_loc_x = max(x)-1;   
    text_loc_y = 90; 
    text_str = sprintf('Simulated ORR [IQR] : %.1f%% [%.1f%% - %.1f%%]', 100*median_orr,100*ORR_CI_95(1),100*ORR_CI_95(2));
    text_str_1 = sprintf('Clinical ORR : %.1f%%', sum(Lesion_Waterfall.OVLRSPN < 2)/height(Lesion_Waterfall) * 100);

    text(text_loc_x, text_loc_y-10, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');
    text(text_loc_x, text_loc_y, text_str_1, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');

    % Cosmetics
    ylabel('% Change in SLD',FontWeight="bold");
    xlabel('Patient No.',FontWeight="bold");
    fontsize(15,'points');
    title('Multi-seed PD Simulation',FontSize = 18,FontWeight="bold");
    legend('Box','off', 'Location', 'southwest','FontSize',13);
    ylim([-105,105]);
    xlim([0.5, Sim_Params.N_VP + 0.5]);

    pbaspect([3 2 1]);
    fontname('Poppins');
    hold off;

    % PSA Waterfall plotting
    figure; 
    hold on;
    box on; set(gca,'LineWidth', 1.25);

    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    if plot_data && height(PSA_Waterfall) == Sim_Params.N_VP
        bar(sort(PSA_Waterfall.BOR,'descend'), 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    end

    psa_waterfall_for_all_seed = min(100, psa_waterfall_for_all_seed);
    medianVals = median(psa_waterfall_for_all_seed, 1);
    prc5 = prctile(psa_waterfall_for_all_seed, 5, 1); 
    prc95 = prctile(psa_waterfall_for_all_seed, 95, 1);

    n=length(medianVals);
    medianCI_upper= nan([1,n]);
    medianCI_lower= nan([1,n]);
    prc5CI_upper= nan([1,n]);
    prc5CI_lower= nan([1,n]);
    prc95CI_upper= nan([1,n]);
    prc95CI_lower= nan([1,n]);

    for i = 1:length(medianVals)
        CI95=boot95CI(psa_waterfall_for_all_seed(:,i),0.5);
        medianCI_lower(i)=CI95(1);
        medianCI_upper(i)= CI95(2);

        CI95=boot95CI(psa_waterfall_for_all_seed(:,i),0.05);
        prc5CI_lower(i)=CI95(1);
        prc5CI_upper(i)= CI95(2);

        CI95=boot95CI(psa_waterfall_for_all_seed(:,i),0.95);
        prc95CI_lower(i)=CI95(1);
        prc95CI_upper(i)= CI95(2);
    end
    x = 1:numel(medianVals);

    % Shaded region (min–max)
    fill([x fliplr(x)], [prc5CI_lower fliplr(prc5CI_upper)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(3), ...
        'HandleVisibility','off');
    fill([x fliplr(x)], [prc95CI_lower fliplr(prc95CI_upper)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(5), ...
        'HandleVisibility','off');
    fill([x fliplr(x)], [medianCI_lower fliplr(medianCI_upper)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',gem_hex_codes(7), ...
        'HandleVisibility','off');

    % Overlay lines
    plot(x, medianVals, '-', 'Color', gem_hex_codes(7), 'LineWidth', 2.0, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Median Simulated Waterfall');
    plot(x, prc95, '-', 'Color', gem_hex_codes(5), 'LineWidth', 2.0, 'DisplayName', '95th percentile');
    plot(x, prc5, '-', 'Color', gem_hex_codes(3), 'LineWidth', 2.0, 'DisplayName', '5th percentile');

    mean_psa50 = mean(psa50_for_all_seed);
    PSA50_CI_95 = prctile(psa50_for_all_seed, [25, 75]);

    % Mention ORR on the plot
    text_loc_x = max(x)-0.5;   
    text_loc_y = 90; 
    text_str = sprintf('Simulated PSA50 [IQR] : %.1f%% [%.1f%% - %.1f%%]', 100*mean_psa50/height(PSA_Waterfall),100*PSA50_CI_95(1)/height(PSA_Waterfall),100*PSA50_CI_95(2)/height(PSA_Waterfall));
    text_str_1 = sprintf('Clinical PSA50 : %.1f%%', sum(PSA_Waterfall.BOR<=-50)/height(PSA_Waterfall) * 100);

    text(text_loc_x, text_loc_y-10, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');
    text(text_loc_x, text_loc_y, text_str_1, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');

    % Cosmetics
    ylabel('% Change in Serum PSA',FontWeight="bold");
    xlabel('Patient No.',FontWeight="bold");
    fontsize(15,'points');
    title('Multi-seed PSA Simulation',FontSize = 18,FontWeight="bold");
    legend('Box','off', 'Location', 'southwest','FontSize',13);
    ylim([-100,100]);
    % xlim([0.5, Sim_Params.N_VP+0.5]);

    pbaspect([3 2 1]);
    fontname('Poppins');
    hold off;
end


function Plot_PFS_for_all_seed(Sim_Params, Dose, pfs_for_all_seed, plot_95_CI, plot_choice, seed_indices)
    %
    % Plot_PFS_for_all_seed: Uses the output PFS curves of the simulated
    %  virtual population for the bunch of seeds and plots the mean and 95%
    %  CI for all plotted PFS curved
    %
    % INPUTS:
    %   Dose:               The dose(s) used for simulated
    %   Sim_Params:         Information about the crucial parameters used during simulation
    %   pfs_for_all_seed:   Stored array containing PFS curve for all seeds
    %   plot_95_CI :        Boolean flag to decide whether or not to plot the 95% CI interval for each PFS curve
    %   plot_choice :       Allowed values are 'mean of seeds'(default value) or 'all seeds' or
    %                       'certain seeds'. To decide whether to plot the mean PFS with 95%tile band or simply
    %                       plot PFS for all simulated seeds or only for few of them, 
    %                       might be useful for eyeballing some particular seed curves
    %   seed_indices :      To be provided when plot_choice is 'certain
    %                       seeds'. It provides the seed for which the PFS 
    %                       data needs to be filtered and plotted.
    %
    % OUTPUTS:     A set of PFS curves for each seeds on a sample graph.
    %
    % Set default plot_choice 
    if nargin < 5 || isempty(plot_choice)
        plot_choice = "mean of seeds";
    end
    if strcmpi(plot_choice,'certain seeds') & isempty(seed_indices)
        error('If you specify the plot choice as "certain seeds", then you must specify the seed_indices as the last argument in the function Plot_PFS_for_all_seed.');
    end
    % Get the color palette
    gem_hex_codes = rgb2hex(orderedcolors("gem"));
    nColors = size(gem_hex_codes, 1);
    if strcmp(plot_choice, "mean of seeds")
        % PLOT ALL SEEDS
        
        % Get dimensions from the main cell array
        nSeeds = length(pfs_for_all_seed);
        
        % Loop over each DOSE - one figure per dose
    
        figure; 
        hold on;
        grid off;
        box on; set(gca, 'LineWidth', 1.25);
        time_array = Sim_Params.Time_efficacy;
        varNames = "T_" + string(time_array);
        f_seeds = array2table(zeros(0, length(time_array)), 'VariableNames',cellstr(varNames));
        % Loop over each SEED for the current dose
        for sIndex = 1:nSeeds
            
            Sim_PFS = pfs_for_all_seed{sIndex};
            
            % Extract data for ecdf
            pfs_times = Sim_PFS(1:end, 2);
            pfs_event_status = Sim_PFS(1:end, 3);
            censoring_flags = ~logical(pfs_event_status);
            
            % Calculate survivor function and bounds
            [f, x, ~, ~] = ecdf(pfs_times, 'Censoring', censoring_flags, ...
                                    'Function', 'survivor', 'Alpha', 0.05, ...
                                    'Bounds', 'on');
            % Dynamically find the maximum follow-up time
            max_time = max(pfs_times);
            
            % Append if the max time is greater than the last event time
            if x(end) < max_time
                x = [x; max_time];
                f = [f; f(end)];
            end

            new_row = interp1(x, f, time_array, 'previous', 'extrap');
            new_row_table = array2table(new_row, 'VariableNames', f_seeds.Properties.VariableNames);
            f_seeds = [f_seeds; new_row_table];    %#ok<AGROW>
        end % End of seed loop
        
        % Plot the Kaplan-Meier estimate
        f_seeds = f_seeds .* 100;
        fseed_matrix = f_seeds{:, :}; 
        
        mean_vals = mean(fseed_matrix, 1);
        p95_vals  = prctile(fseed_matrix, 95, 1);
        p5_vals   = prctile(fseed_matrix, 5, 1);
        [x_stairs, p95_stairs] = stairs(time_array, p95_vals);
        [~, p5_stairs] = stairs(time_array, p5_vals);
        x_fill = [x_stairs; flipud(x_stairs)];
        y_fill = [p95_stairs; flipud(p5_stairs)];
        fill(x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
            'DisplayName', '5th - 95th Percentile');
        
        % Plot the mean PFS
        stairs(time_array, mean_vals, 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Mean PFS');       

        % Cosmetics for the DOSE figure
        xlabel('Time (days)',FontWeight="bold");
        ylabel('PFS Probability (%)',FontWeight="bold");
        fontsize(15,'points'); 
        title(sprintf('Simulated PFS with %d seeds (Dose = %d mg/kg)', nSeeds, Dose), FontSize=18, FontWeight="bold");
        ylim([0 102]); xlim([0,366]);
        legend('Box','on', 'Location', 'northeast', FontSize= 15);
        pbaspect([3 2 1]);
        fontname('Poppins');
        hold off;

    elseif strcmp(plot_choice, "all seeds")
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
        title(sprintf('Simulated PFS with %d seeds (Dose = %d mg/kg)', nSeeds, Dose), FontSize=18, FontWeight="bold");
        ylim([0 102]); xlim([0,366]);
        legend('Box','on', 'Location', 'northeast', FontSize= 15);
        pbaspect([3 2 1]);
        fontname('Poppins');
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
        title(sprintf('Simulated PFS (Dose = %d mg/kg; Seed Indices = %s)', Dose, num2str(seed_indices)), FontSize=16, FontWeight="bold");
        ylim([0 100]); xlim([0,Sim_Params.StopTime+1]);
        legend('Box','off', 'Location', 'best', FontSize= 12);
        pbaspect([3 2 1]);
        fontname('Poppins');
        hold off;
          
    else
        error('Invalid plot_choice. Must be "mean of seeds" or "all seeds" or "certain seeds".');
    end

end

function Plot_dose_response(All_Dose_resp)
    %
    % plot_dose_response: Plots box plot indicating the dose vs response trend
    %
    % INPUTS:
    %   All_Dose_resp:    The array containing responses across all doses
    %                     (Column 1: Doses, Column 2: Responses)
    %
    % OUTPUTs:
    %   A box chart showing dose vs response trend
    %
    figure; hold on;
    % mean6 = mean(PSA_Waterfall.BOR(PSA_Waterfall.Dose_mpk==6));
    % mean8 = mean(PSA_Waterfall.BOR(PSA_Waterfall.Dose_mpk==8));
    % plot([6,8],[mean6,mean8],'--d', 'Color', 'k', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6,'DisplayName','Clinical mean BOR');
    % 
    % a = boxchart(PSA_Waterfall.Dose_mpk-0.2,PSA_Waterfall.BOR, ...
    %     'MarkerStyle', 'o', ...
    %     'MarkerSize', 5, ...
    %     'BoxFaceColor', 'none', 'BoxEdgeColor','k', ...
    %     'BoxWidth', 0.15,'DisplayName','Clinical Data'); 
    % 
    % a.LineWidth = 1.5;
    
    b = boxchart(All_Dose_resp(:,1), All_Dose_resp(:,2), ...
        'MarkerStyle', 'o', 'MarkerSize', 5, ...
        'BoxFaceColor', '#DD5400', ...
        'BoxWidth', 0.15,'DisplayName','Simulated Data', 'WhiskerLineColor','#DD5400');
    
    b.LineWidth = 1.5;
    
    doses_x = unique(All_Dose_resp(:,1));
    medians_y = arrayfun(@(d) mean(All_Dose_resp(All_Dose_resp(:,1)==d, 2), 'omitnan'), doses_x);
    plot(doses_x, medians_y, '--s', 'Color', '#DD5400', 'LineWidth', 1, 'MarkerFaceColor', '#DD5400', 'MarkerSize', 6,'DisplayName','Simulated mean BOR');
    
    xlabel('Dose (mg/kg)', 'FontWeight', 'bold');
    ylabel('Best % Change in SLD', 'FontWeight', 'bold');
    grid on; box on;
    set(gca, 'LineWidth', 1.25, 'FontName', 'Poppins', 'FontSize', 15);
    title('Dose vs Best Response', 'FontSize', 16, 'FontWeight', 'bold');
    legend();
    hold off;
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

function CI95 = boot95CI(array_to_sample_from,prctile)
    arr = nan([3000,1]);
    for i=1:3000
        arr(i,1) = quantile(randsample(array_to_sample_from,length(array_to_sample_from),true),prctile);
    end
    CI95 = [quantile(arr,0.025);quantile(arr,0.975)];
end
 
function [mean_p, median_p, ci_mean, ci_median, pi_p] = analyze_binomial_p(N, p_array)
    % Ensure the array is a column vector
    p_array = p_array(:);
    n_samples = length(p_array); % Should be 100 based on your prompt
    
    % 1. Mean
    mean_p = mean(p_array);
    
    % 2. Median
    median_p = median(p_array);
    
    % --- Confidence Intervals ---
    % Calculate critical t-value for 95% CI (alpha = 0.05, two-tailed)
    alpha = 0.05;
    t_crit = tinv(1 - alpha/2, n_samples - 1);
    z_crit = norminv(1 - alpha/2); % ~1.96 for large N / theoretical formulas
    
    % 3. 95% CI of the Mean
    % Formula: Mean +/- t * (Standard Deviation / sqrt(n))
    std_p = std(p_array);
    standard_error_mean = std_p / sqrt(n_samples);
    margin_error_mean = t_crit * standard_error_mean;
    ci_mean = [mean_p - margin_error_mean, mean_p + margin_error_mean];
    
    % 4. 95% CI of the Median
    % Using the asymptotic standard error of the median for normal-ish data:
    % SE_median ~= 1.2533 * SE_mean
    standard_error_median = 1.2533 * standard_error_mean;
    margin_error_median = t_crit * standard_error_median;
    ci_median = [median_p - margin_error_median, median_p + margin_error_median];
    
    % 5. 95% Prediction Interval of a new p-value
    % Since each p is from a binomial distribution out of N, the variance 
    % of a new p-value (X/N) is expected to be p(1-p)/N. 
    % We use the mean_p as our best estimate for the true probability.
    expected_variance = (mean_p * (1 - mean_p)) / N;
    
    % The prediction interval combines the uncertainty of our mean estimate 
    % AND the natural binomial variance of a new observation.
    % PI = Mean +/- z * sqrt( SE_mean^2 + Binomial_Variance )
    pi_margin = z_crit * sqrt(standard_error_mean^2 + expected_variance);
    
    % Bound the PI between 0 and 1 since p cannot physically exceed these bounds
    lower_pi = max(0, mean_p - pi_margin);
    upper_pi = min(1, mean_p + pi_margin);
    pi_p = [lower_pi, upper_pi];
    
    % --- Display Results (Optional, for convenience) ---
    fprintf('Mean of p: %.4f\n', mean_p);
    fprintf('Median of p: %.4f\n', median_p);
    fprintf('95%% CI for Mean: [%.4f, %.4f]\n', ci_mean(1), ci_mean(2));
    fprintf('95%% CI for Median: [%.4f, %.4f]\n', ci_median(1), ci_median(2));
    fprintf('95%% Prediction Interval for a new p: [%.4f, %.4f]\n', pi_p(1), pi_p(2));
end
 