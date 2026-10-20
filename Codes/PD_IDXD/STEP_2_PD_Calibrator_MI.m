clear; warning off; parfevalOnAll(@warning,0,'off','all');
load("LesionWaterfall.mat"); 
load("PSAWaterfall.mat"); 
Sim_Params = struct(); % struct to pass all crucial choice related parameters into the simulating function

%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Automated_calibration_switch = true;                   % If true, the calibration will be done by the solver, you will need to provide the parameter initial values and their bounds.
Calibrating_what = 'LESION';                           % Valid inputs are 'LESION' or 'PSA'.

% Load Project model
proj = sbioloadproject("../../GSK227_CRPC.sbproj");     % This path has to be changed as per the location of the model
model = proj.m1; 

Sim_Params.seed = 0;                                   % This is the RNG Seed value to which you are fitting
Sim_Params.Dose = LesionWaterfall.Dose_mpk;             % in mg/kg , should be matching with the dose used in Step_1_Data_processor script
Sim_Params.N_VP = height(LesionWaterfall.BOR);          % Number of patients to simulate is same as that in the data
Sim_Params.variability = 1;                             
Sim_Params.StopTime = 395;                              % Simulating for 1 year
Sim_Params.Time_Sim = 0:49:Sim_Params.StopTime;         % for each of the days, time point where the simulations will be carried out
Sim_Params.Dose_Cycle_length = 21;                      % in days, 21 Days for Q3W
Sim_Params.Doses_days_in_a_cycle = 0;                   % On what days you want to administer the dosage, it could be [0,7] if you want d1,d8 Q3W kind of dosing
Sim_Params.Dose_target_name = 'Central_ADC.Cen_ADC_nM'; % Species name in the model where ADC is dosed
Sim_Params.Dar_Dose_target_name = 'Central_ADC.Cen_DummyDAR'; % Species name in the model where DAR is dosed
Sim_Params.Scan_Interval = 49;                          % Time interval between two radiographical scan in simulation, (in days)
Sim_Params.Time_efficacy=0:Sim_Params.Scan_Interval:Sim_Params.StopTime; % response measument days 

Sim_Params.UseParallel_simfcn = true;

Sim_Params.StatesToLog = {'TV','PSA', 'Cen_ADC_ugml', 'Cen_tAb_ugml', 'Cen_PL_ngml'};
Sim_Params.parameter_name_for_PSA = 'PSA';
Sim_Params.parameter_name_for_tumor_volume = 'TV';
Sim_Params.parameter_name_for_hscore = 'H_Score';
Sim_Params.parameter_name_for_bodyweight = 'W';
Sim_Params.parameter_name_for_mol_weight_ADC = 'MW_ADC';
Sim_Params.N_tmr = 1;
Sim_Params.Param_names_for_intrapatient_heterogeneity = {'Tumor_Dia','H_Score','K1','Kmax_p'}; % Parameters that will change across within each patient
Sim_Params.Intrapatient_Heterogeneity = 0.0;       % Variability in parameters across different tumors within patient
Sim_Params.dropout_Switch = true;                  % Switch to turn off/on the dropouts due to New metastats or Non-specific reasons
Sim_Params.NM_Params = [0.01 0.21 15 1.0];         % New-metastates related dropout probability modelling parameters
Sim_Params.NonSpecific_dropouts = 0.21;            % Probability for dropout related to non-radiographical reasons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-exploratory params, that are different from model values
model_params_drug_specific = {
    %param_name    %mu        %sig
    'MW_ADC'       147000     0;
    'MW_PL'        494        0;
    'DAR_t0'       4          0;
    'Kd_ADC'       0.8        0.1;
    'Kon_ADC'      8.64       0.12;
    };

model_params_indication_specific = {
    %param_name            %mu        %sig
    'Tumor_Dia'            49         1.2;
    'H_Score'              150        0.475;
    'Free_fraction_PL'     1.2        0;
    'Ag_p'                 220000     0.2;
    'Ag_n'                 0          0.0;
    'Single_Cell_Vol_p'    0.27e-12   0.0;
    'K_Deg_Lyso'           0.69       0;
    'Kin_PL_p'             273        0;
    'Kin_PL_n'             273        0;
    'Kint_ADC'             1.944      0;
    'KC50_p'               0.49       0.15;
    'K1'                   0.0498      0.15;
    'Kmax_p'               0.0615      0.17;
    'Cp'                   1.9          0.3;
    % 'n'                    2.5        0.2;
    'CL_PSA'               0.7e-4     0.8;
    'Initial_PSA_ngml'     76         1.1;
    'k_res'                0.015      0.2;
    };

model_params_patient_specific = {
    %param_name    %mu         %sig
    'V1_ADC'       0.039035    0.2;
    'V2_ADC'       0.019844    0.2;
    'CL_ADC'       0.00524     0.2;
    'CLD_ADC'      0.010561    0.2;
    'K_P_dec'      0.03        0.0;
    'V1_PL'        0.04        0.2;
    'V2_PL'        0.07        0.2;
    'CL_PL'        3.5         0.2;
    'CLD_PL'       0.003       0.2
    };

model_params_ne = { % all the parameters that you haven't specified above and you are not estimating/or are not properly set in the model as default value
    %param_name          %mu    %sig
    'W'                  70     0;
    'Switch_TumorComp'   1      0;
    'Switch_Tox'         1      0
    };

model_params_ne = [
    model_params_ne;
    model_params_drug_specific;
    model_params_indication_specific;
    model_params_patient_specific
    ];

model_params_ne = cell2table(model_params_ne,'VariableNames',{'param_name','mu','sig'});

% Exploratory params (can also be called parameters to estimate)
model_params_e = { 
    %param_name     %statistic   %initial_value    %lower_bound  %upper_bound  
    % 'CL_PSA'          'mu'           0.6e-4             1e-5       1.5e-3;
    % 'Cp'              'mu'           3.5                1          6;
    % 'Cp'              'sig'          0.3               0.001         1;  
    % 'CL_PSA'          'sig'          0.8               0.001         1;
    % 'k_res'             'mu'         0.015             1e-4          0.03;
    % 'k_res'            'sig'         0.2                0.001         1;  
    'K1'              'sig'        0.15               0.01          1;
    'Kmax_p'          'mu'         0.06094               0.001         1;
    'Kmax_p'          'sig'        0.17               0.01          1;
    'Free_fraction_PL'  'mu'         1.2               0.001         10;
    };
if ~isempty(model_params_e)
    model_params_e = cell2table(model_params_e,'VariableNames',{'param_name','statistic','initial_value','lower_bound','upper_bound'});
else
    model_params_e = array2table(zeros(0,5),'VariableNames',{'param_name','statistic','initial_value','lower_bound','upper_bound'});
end
%%%Note : 
% 1) If a parameter's sig is specified in model_params_e
%    but mu is neither specified in model_params_e nor in model_params_ne,
%    then the MODEL VALUE is taken as the parameter's mu
% 2) If a parameter's mu is specified in model_params_e
%    but sig is neither specified in model_params_e nor in model_params_ne,
%    then the parameter's sig is taken as ZERO

% Combine non-exploratory and exploratory parameters
[model_params,indices_e] = combine_params(model_params_e,model_params_ne,model);

initial_values = model_params_e.initial_value;    

if Automated_calibration_switch
    
    lower_bound    = model_params_e.lower_bound; %#ok<*UNRCH>
    upper_bound    = model_params_e.upper_bound;

    objFun=@(fit_Params) objfun(fit_Params, model_params,indices_e, model, Sim_Params, LesionWaterfall, PSAWaterfall, Calibrating_what);
    fprintf('Initiating automated Calibration with %d parameters...\n\n', length(initial_values));
    fminconopts=optimoptions('fmincon','Display', 'iter-detailed', 'MaxIterations', 20, 'Algorithm', 'interior-point', 'UseParallel',true,'PlotFcn','optimplotfval');
    tic;
    [bestParams,bestError, exitflag, output] = fmincon(objFun, initial_values, [], [], [], [], lower_bound, upper_bound, [], fminconopts);
    time_taken = toc;
    fprintf('Took %.1f seconds to calibrate.\n', time_taken);

else
    %%%%%%%%%%%%%%%%% USE THIS PART FOR MANUAL FITTING %%%%%%%%%%%%%%%%%%%
    bestParams = initial_values;   
end

Sim_Params.param = update_model_params_from_e(model_params,bestParams,indices_e);

fprintf('================================================\n');
fprintf('Calibrated the parameters:\n');
fields = model_params_e.param_name; % Get list of parameter names
for i = 1:numel(fields)
    fprintf('%s : %.4g\n', fields{i}, bestParams(i));
end
fprintf('================================================\n');
fprintf('\nWill use these parameters to simulate and visualise the waterfall calibration\n\n');
save('Sim_Params.mat','Sim_Params');
fprintf('Saved the parameter struct for further visualisation.\n\n');


%%%%%%%%%%%%%%%%%%% SIMULATE WITH ESTIMATED PARAMETERS %%%%%%%%%%%%%%%%%%%%
rng(Sim_Params.seed, 'twister');    % Set RNG seed for reproducibility
tic;
Simulation_Results = run_population_simulation(model, Sim_Params, true);
time_taken = toc;
fprintf('Took %.1f seconds to calibrate.', time_taken);

%%
plot_waterfall(Simulation_Results, LesionWaterfall, PSAWaterfall); 
% plot_pfs(Simulation_Results);                             % [OPTIONAL] Lines marked as optional can be commented without any issue
% plot_swimmer(Simulation_Results);                         % [OPTIONAL]
% plot_spider(Simulation_Results, Sim_Params);                % [OPTIONAL]
plot_psa_spider(Simulation_Results, Sim_Params);            % [OPTIONAL]
% plot_observables(Simulation_Results, Sim_Params)          % [OPTIONAL] 
% summarize_patient_response(Simulation_Results);           % [OPTIONAL]








































%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective Function 

function err = objfun(fit_Params, model_params,indices_e, model, Sim_Params, LesionWaterfall, PSAWaterfall, Calibrating_what) %#ok<DEFNU>
    %
    % Defines the objective function that will be minimized for automated calibration process.
    %
    % INPUTS:
    %       fit_Params:         The value for the parameters (to be estimated) at which they will be evaluated in an iteration
    %       model_params:       Set of values for all the parameters to be set and simulated in the model
    %       indices_e:          Indices of the parameters(to be estimated) in the model_params, so that they can replaced
    %       model:              A SimBiology model object (to be simulated)
    %       Sim_Params:         Information about the crucial parameters to use while simulation
    %       LesionWaterfall:    Processed data containing information about the lesion waterfall plot from the raw data
    %       PSAWaterfall:       Processed data containing information about the PSA waterfall plot from the raw data
    %       Calibrating_what:   'PSA' or 'LESION'. Will use the corresponding waterfalls to calibrate
    %       
    %
    % OUTPUTS:
    %       err:                The square root of sum of square of differences between the Clinical waterfall and the simulated waterfall
    %                 
    
    
    rng(Sim_Params.seed, "twister");    % Set RNG seed for reproducibility

    % Pass the parameters
    Sim_Params.param = update_model_params_from_e(model_params,fit_Params,indices_e);
 
    Simulation_Results = run_population_simulation(model, Sim_Params, false); % Run simulation 
    simulated_waterfall_Data = Simulation_Results.AllPlots.Sim_Waterfall;      % Extract the Simulated Waterfall data
    simulated_waterfall_Data = sortrows(simulated_waterfall_Data, -2);

    if strcmp(Calibrating_what,'PSA')
        psawaterfall = simulated_waterfall_Data(end-height(Waterfall)+1:end,4);
        y = sort(psawaterfall, 'descend');
        % y = simulated_waterfall_Data(:,2); 
        obs=sort(PSAWaterfall.BOR,'descend');     % The clinical waterfall data
    elseif strcmp(Calibrating_what, 'LESION')
        obs=sort(LesionWaterfall.BOR,'descend');     % The clinical waterfall data
        y = simulated_waterfall_Data(:,2);
    else
        error('Wrong input for the parameter "Calibrating_what" in objfun. Use either "PSA" or "LESION".')
    end

    Diff= y-obs;
    err=sum(Diff.^2); % Calculating the error metric to minimise
end

% 
% function Simulation_Results = run_population_simulation(model, Sim_Params, verbose)
%     %
%     % run_population_simulation: Simulates a virtual population for a given mode and returns all processed results for plotting.
%     %
%     % INPUTS:
%     %   model:              A SimBiology model object (pre-loaded)
%     %   Sim_Params:         Information about the crucial parameters to use while simulation
%     %   verbose:            If true, prints information while running
%     %
%     % OUTPUTS:
%     %   Simulation_Results:     A struct containing all relevant simulated outputs:
%     %                       .AllPlots (Tumor trajectory, Waterfall, PFS, Dropout related informat)
%     %                       .All_PK/ADC/TAB (array for PK output)
%     %                       .VP_Init_Array (Array containing the initial patient related parameters)
%     %                       .Dose (the dose amount in mg/Kg)
% 
% 
%     if verbose
%         disp('Running Simulation with estimated parameters...');
%     end
% 
% 
%     % Configure SimBiology Model
%     csObj = getconfigset(model);
%     csObj.SolverType = 'sundials';
%     csObj.SolverOptions.RelativeTolerance = 1e-5;
%     csObj.SolverOptions.AbsoluteTolerance = 1e-5;
%     csObj.SolverOptions.MaxStep = 1;
%     csObj.SolverOptions.AbsoluteToleranceScaling = false;
%     csObj.StopTime = Sim_Params.StopTime;
%     csObj.TimeUnits = 'day';
%     csObj.SolverOptions.OutputTimes = Sim_Params.Time_Sim;
%     csObj.RuntimeOptions.StatesToLog = Sim_Params.StatesToLog;
% 
%     % Add Dummy Dosing for ensuring proper DAR resetting at each ADC Dose adminisation
%     d2 = sbiodose( 'DAR_Dose', 'repeat');
%     d2.Amount = 1;
%     d2.AmountUnits = 'mole';
%     d2.TimeUnits = 'day';
%     d2.StartTime = 0;
%     d2.Interval = Sim_Params.Dose_Cycle_length;
%     d2.RepeatCount = round(Sim_Params.StopTime/Sim_Params.Dose_Cycle_length) + 1;
%     d2.TargetName = Sim_Params.Dar_Dose_target_name;
%     DAR_doses = [];
%     for dose_day = Sim_Params.Doses_days_in_a_cycle
%         DAR_dose = copyobj(d2);
%         DAR_dose.Name = sprintf('DAR_dose_day_%d', dose_day);
%         DAR_dose.StartTime = dose_day;
%         DAR_doses = [DAR_doses, DAR_dose]; %#ok<AGROW>
%     end
%     d2 = DAR_doses; clear DAR_dose DAR_doses
% 
%     IHC_Category = nan(Sim_Params.N_VP,1);
% 
%     % Initialising VP_Init_Array
%     p_names = fieldnames(Sim_Params.param)';
%     all_fields = [{'PatientID', 'TumorID'}, p_names];
% 
%     % Preallocate the full column vector (Total_Rows x 1)
%     % This creates the full structure array filled with empty values
%     VP_Init_Array = repmat(cell2struct(cell(length(all_fields), 1), all_fields), Sim_Params.N_VP * Sim_Params.N_tmr, 1);
% 
% 
% 
%     % ADC Dosing Setup
%     d1 = sbiodose('ADC_dose', 'repeat');
%     d1.AmountUnits = 'nanomole';
%     d1.TimeUnits = 'day';
%     d1.StartTime = 0;
%     d1.Interval = Sim_Params.Dose_Cycle_length;
%     d1.RepeatCount = round(Sim_Params.StopTime/Sim_Params.Dose_Cycle_length) + 1;
%     d1.TargetName = Sim_Params.Dose_target_name;
%     ADC_doses = [];
%     for dose_day = Sim_Params.Doses_days_in_a_cycle
%         ADC_dose = copyobj(d1);
%         ADC_dose.Name = sprintf('ADC_dose_day_%d', dose_day);
%         ADC_dose.StartTime = dose_day;
%         ADC_doses = [ADC_doses , ADC_dose]; %#ok<AGROW>
%     end
%     d1 = ADC_doses; clear ADC_dose ADC_doses
% 
% 
% 
%     % Initialize Dose-Level Arrays
%     n_timepoints = numel(Sim_Params.Time_efficacy);
%     SLD = nan(Sim_Params.N_VP * n_timepoints, 6);
%     Sim_Vol = cell(Sim_Params.N_VP, 1); 
%     PSA = cell(Sim_Params.N_VP, 1); 
%     Dia = cell(Sim_Params.N_VP, 1); 
%     dTS = cell(Sim_Params.N_VP, 1);
%     total_indts_rows = Sim_Params.N_VP * Sim_Params.N_tmr * n_timepoints;
%     IndTS = nan(total_indts_rows, 5);
%     indts_idx = 1; % Counter to track where to write in IndTS
%     time_series_data_for_all_pts = struct();
%     for i = 1:length(Sim_Params.StatesToLog)
%         var_name = sprintf('sim_%s', Sim_Params.StatesToLog{i});
%         time_series_data_for_all_pts.(var_name) = NaN(Sim_Params.N_VP, numel(Sim_Params.Time_Sim));
%     end  
% 
% 
%     % Patient Simulation Loop
%     for patient_id = 1:Sim_Params.N_VP
% 
%         if verbose
%             fprintf('Simulating for patient no. %d\n',patient_id);
%         end
% 
%         % Preallocate matrices inside the cells for this specific patient
%         Sim_Vol{patient_id} = nan(n_timepoints, Sim_Params.N_tmr);
%         Dia{patient_id}     = nan(n_timepoints, Sim_Params.N_tmr);
%         dTS{patient_id}     = nan(n_timepoints, Sim_Params.N_tmr);
%         PSA{patient_id}     = nan(n_timepoints, Sim_Params.N_tmr);
% 
%         [vp_init_Array, ihc_cat] = generate_virtual_patient(Sim_Params);
%         [vp_init_Array.PatientID] = deal(patient_id);
%         vp_init_Array.Tumor_Dia = min(150, vp_init_Array.Tumor_Dia);
%         idx_start = (patient_id-1) * Sim_Params.N_tmr + 1;
%         idx_end   = patient_id * Sim_Params.N_tmr;
%         % Insert directly into the preallocated slots
%         VP_Init_Array(idx_start:idx_end) = vp_init_Array;
% 
%         IHC_Category(patient_id) = ihc_cat;
%         Dose_nanomoles = (Sim_Params.Dose(patient_id) * vp_init_Array(1).(Sim_Params.parameter_name_for_bodyweight) / (1000 * vp_init_Array(1).(Sim_Params.parameter_name_for_mol_weight_ADC))) * 10^9;  % Convert mpk to nanomole
%         for dose_num = 1:length(d1)
%             d1(dose_num).Amount = Dose_nanomoles;
%         end
% 
%         param_names = fieldnames(vp_init_Array(1));
%         for tumor_id = 1:length(vp_init_Array)
%             % Set patient-specific (variable) parameters
%             for x = 3:length(param_names)
%                 pObj = sbioselect(model, 'Type', 'parameter', 'Name', param_names{x});
%                 pObj.Value = vp_init_Array(tumor_id).(param_names{x});
%             end        
%             % Run simulation
%             simData = sbiosimulate(model, csObj, [], [d1 d2]);
%             Time_Sim = simData.Time;
% 
%             % Extract data
%             TV_all = selectbyname(simData, Sim_Params.parameter_name_for_tumor_volume).Data;  
%             PSA_all = selectbyname(simData, Sim_Params.parameter_name_for_PSA).Data;
%             % Subsample for efficacy
%             idx_efficacy = ismember(Time_Sim, Sim_Params.Time_efficacy); 
%             TV = TV_all(idx_efficacy);
%             PSA{patient_id}(:,tumor_id) = PSA_all(idx_efficacy);
% 
%             Sim_Vol{patient_id}(:,tumor_id) = TV;
%             Dia{patient_id}(:,tumor_id) = (6 * Sim_Vol{patient_id}(:,tumor_id) / pi).^(1/3);
%             dTS{patient_id}(:,tumor_id) = 100 * (Dia{patient_id}(:,tumor_id) - Dia{patient_id}(1,tumor_id)) / Dia{patient_id}(1,tumor_id);
%             lesion_details = [repmat(patient_id, n_timepoints, 1), Sim_Params.Time_efficacy', repmat(tumor_id, n_timepoints, 1), Dia{patient_id}(:,tumor_id), dTS{patient_id}(:,tumor_id)];
%             IndTS(indts_idx : indts_idx + n_timepoints - 1, :) = lesion_details;
%             indts_idx = indts_idx + n_timepoints;
%         end    
% 
%         % Aggregate at patient level
%         VP_Diam = sum(Dia{patient_id}, 2);
%         Temp_SLD_per_pts = [repmat(patient_id, n_timepoints, 1), Sim_Params.Time_efficacy', VP_Diam];
%         Temp_SLD_per_pts(:,4) = 100 * (Temp_SLD_per_pts(:,3) - Temp_SLD_per_pts(1,3)) / Temp_SLD_per_pts(1,3);
%         Temp_SLD_per_pts(:,5) = PSA{patient_id};
%         Temp_SLD_per_pts(:,6) =  100 * (Temp_SLD_per_pts(:,5) - Temp_SLD_per_pts(1,5)) / Temp_SLD_per_pts(1,5);
% 
%         SLD(1+(patient_id-1)*n_timepoints:patient_id*n_timepoints,:) = Temp_SLD_per_pts;
% 
%         % Store extra logged time_series
%         for i = 1:length(Sim_Params.StatesToLog)
%             var_name = sprintf('sim_%s', Sim_Params.StatesToLog{i});
%             time_series_data_for_all_pts.(var_name)(patient_id,:) = selectbyname(simData, Sim_Params.StatesToLog{i}).Data'; % Take transpose to ensure row vector
%         end
%     end
% 
%     % Post-Processing:  BOR, PFS, and Dropout calculation  
%     % Apply Dropout Rules
%     Sim_SLD = apply_dropouts(Sim_Params, SLD, IndTS);
%     PatientList = (1:Sim_Params.N_VP)';
%     % Calculate Best Overall Response (for Waterfall)
%     Sim_Waterfall = calculate_bor(Sim_SLD, PatientList);
% 
%     % Calculate PFS
%     Sim_PFS = calculate_pfs(Sim_SLD, Sim_Waterfall);
% 
%     % Calculate Dropout Times 
%     DropoutTime = calculate_dropout_times(Sim_SLD, PatientList);
% 
%     % Calculate Dropout Reasons
%     DropoutSummary = calculate_dropout_reasons(Sim_SLD, PatientList);
% 
%     % Store All Processed Data for Plotting 
%     Simulated_Data = struct();
%     Simulated_Data.Sim_SLD = Sim_SLD;
%     Simulated_Data.IHC_Category = IHC_Category;
%     Simulated_Data.DropoutTime = DropoutTime;
%     Simulated_Data.DropoutSummary = DropoutSummary;
%     Simulated_Data.Sim_PFS = Sim_PFS;
%     Simulated_Data.Sim_Waterfall = Sim_Waterfall;
% 
% 
%     % Bundle Final Results
%     Simulation_Results = struct();
%     Simulation_Results.AllPlots = Simulated_Data;
%     Simulation_Results.Dose = Sim_Params.Dose;
%     Simulation_Results.time_series_data_for_all_pts = time_series_data_for_all_pts;
%     Simulation_Results.VP_Init_Array = struct2table(VP_Init_Array);
% end
% 

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
        Temp_SLD_per_pts(:,6) =  100 * (Temp_SLD_per_pts(:,5) - Temp_SLD_per_pts(1,5)) / Temp_SLD_per_pts(1,5);

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
        p_sig = params.(p_name).sig*Sim_Params.variability;
        
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
        
        Sim_Waterfall(patient_id,2) = best_change_TV;
        if isempty(best_time_TV)
            disp(best_time_TV)
            best_time_TV = 0;
        end
        if isempty(best_time_PSA)
            disp(best_time_PSA)
            best_time_PSA = 0;
        end
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
    psawaterfall = min(100, psawaterfall);
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

function plot_spider(Simulation_Results, Sim_Params)
    %
    % plot_spider: Generates a spider plot (tumor trajectories) for each dose.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %


    % Define colors for IHC categories
    IHC_colors = [0.3 0.3 0.3;   % Gray for IHC0
                  0.0 0.6 0.9;   % Blue for IHC1+
                  0.3 0.8 0.3;   % Green for IHC2+
                  0.9 0.2 0.2];  % Red for IHC3+

    AllPlots = Simulation_Results.AllPlots;
    % Extract data for this dose
    Sim_SLD = AllPlots.Sim_SLD;
    DropoutTime = AllPlots.DropoutTime;
    IHC_Category = AllPlots.IHC_Category;
    
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
        IHC_class = IHC_Category(pid); 
        plot(data_cut(:,2), data_cut(:,4), '-', 'Color', IHC_colors(IHC_class+1,:), 'LineWidth', 2);
    end

    % Formatting
    h_leg = zeros(4,1);
    for k = 1:4
        h_leg(k) = plot(nan, nan, '-', 'Color', IHC_colors(k,:), 'LineWidth', 4); % bogus plots to ensure the IHC Colors are shown only once no matter how many curves are there on the spider plot
    end
    
    xlabel('Time (days)', 'FontWeight','bold');
    ylabel('% Change in SLD', 'FontWeight','bold');
    fontsize(15, 'points');
    text(10, 80 , sprintf('Mean Depth of Response at 1st Clinical Visit : %.2f', mean(response_at_1st_visit)), 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'left');
    yline(0,'--k');
    title('Simulated Tumor Trajectories');
    ylim([-100 100]);
    xlim([0, Sim_Params.StopTime + 5]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([ 3 2 1]);
    legend(h_leg, {'IHC0','IHC1+','IHC2+','IHC3+'}, 'Location','northeast');
    fontname('Consolas');
    hold off;

end

function plot_psa_spider(Simulation_Results, Sim_Params)
    %
    % plot_psa_spider: Generates a spider plot (psa trajectories) for each dose.
    %
    % INPUTS:
    %       Simulation_Results :  A struct containing all relevant simulated outputs
    %


    % Define colors for IHC categories
    IHC_colors = [0.3 0.3 0.3;   % Gray for IHC0
                  0.0 0.6 0.9;   % Blue for IHC1+
                  0.3 0.8 0.3;   % Green for IHC2+
                  0.9 0.2 0.2];  % Red for IHC3+

    AllPlots = Simulation_Results.AllPlots;
    % Extract data for this dose
    Sim_SLD = AllPlots.Sim_SLD;
    DropoutTime = AllPlots.DropoutTime;
    IHC_Category = AllPlots.IHC_Category;
    
    PatientList = unique(Sim_SLD(:,1));
    nPatients   = numel(PatientList);        
    
    figure; 
    hold on; 

    for i = 1:nPatients
        pid = PatientList(i);
        data_i = Sim_SLD(Sim_SLD(:,1) == pid, :);
        t_drop = DropoutTime(i);
        data_cut = data_i(data_i(:,2) <= t_drop, :);
        % response_at_1st_visit(i) = data_cut(2,4);

        % Find this patient's IHC category
        IHC_class = IHC_Category(pid); 
        if min(data_cut(2:end,4))<=0
            plot(data_cut(:,2), data_cut(:,6), '-', 'Color', IHC_colors(IHC_class+1,:), 'LineWidth', 2);
        end
    end

    % Formatting
    h_leg = zeros(4,1);
    for k = 1:4
        h_leg(k) = plot(nan, nan, '-', 'Color', IHC_colors(k,:), 'LineWidth', 4); % bogus plots to ensure the IHC Colors are shown only once no matter how many curves are there on the spider plot
    end
    
    xlabel('Time (days)', 'FontWeight','bold');
    ylabel('% Change in PSA', 'FontWeight','bold');
    fontsize(15, 'points');
    yline(0,'--k');
    title('Simulated PSA trajectories');
    ylim([-100 100]);
    xlim([0, Sim_Params.StopTime]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([ 3 2 1]);
    legend(h_leg, {'IHC0','IHC1+','IHC2+','IHC3+'}, 'Location','northeast');
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

function [model_params,indices_e] = combine_params(model_params_e,model_params_ne,model)
      
    % Iterate through the parameters in model_params_e, and add it to model_params one-by-one
    % First split the parameters based on mu and sig
    mu_model_params_e = model_params_e(strcmp(model_params_e.statistic,'mu'),:);
    sig_model_params_e = model_params_e(strcmp(model_params_e.statistic,'sig'),:);
    
    n_new_rows = height(mu_model_params_e) + height(sig_model_params_e);
    % Create a cell buffer to hold: {ParamName, Mu, Sig}
    new_params_buffer = cell(n_new_rows, 3); 
    buffer_idx = 1;

    % Iterate over the mu values of exploratory parameters
    for param_num = 1:height(mu_model_params_e)
        current_param_name = mu_model_params_e.param_name{param_num};
        current_mu_val = mu_model_params_e.initial_value(param_num);
    
        % Check if sig value is provided for this param in model_params_e
        req_row_number = strcmp(sig_model_params_e.param_name,current_param_name);
        if sum(req_row_number)~=0
            current_sig_val = sig_model_params_e.initial_value(req_row_number);       
        else
            % Check if sig value is provided for this param in model_params_ne
            req_row_number = strcmp(model_params_ne.param_name,current_param_name);
            if sum(req_row_number)~=0
                current_sig_val = model_params_ne.sig(req_row_number);
            else
                current_sig_val = 0; % Set the sigma to zero in case it is not provided either in exploratory or non-exploratory
            end
        end
        new_params_buffer(buffer_idx, :) = {current_param_name, current_mu_val, current_sig_val};
        buffer_idx = buffer_idx + 1;    
    end
    
    % Iterate over the sig values of exploratory parameters
    for param_num = 1:height(sig_model_params_e)
        current_param_name = sig_model_params_e.param_name{param_num};
        current_sig_val = sig_model_params_e.initial_value(param_num);  
    
        % Check if mu value is provided for this param in model_params_e
        req_row_number = strcmp(mu_model_params_e.param_name,current_param_name);
        if sum(req_row_number)~=0
            current_mu_val = mu_model_params_e.initial_value(req_row_number); 
        else
            % Check if sig value is provided for this param in model_params_ne
            req_row_number = strcmp(model_params_ne.param_name,current_param_name);
            if sum(req_row_number)~=0
                current_mu_val = model_params_ne.mu(req_row_number);
            else
                current_mu_val = sbioselect(model,'Name',current_param_name);
                if isempty(current_mu_val)
                    error('Parameter %s, for which sigma is provided in exploratory params, does not exist in model',...
                        current_param_name);
                end
                warning('The mean value for Parameter %s is being set to the model default value!!!',current_param_name);
                current_mu_val = current_mu_val.Value;
            end
            

        end

    
        new_params_buffer(buffer_idx, :) = {current_param_name, current_mu_val, current_sig_val};
        buffer_idx = buffer_idx + 1;   
    end
    new_params_table = cell2table(new_params_buffer, 'VariableNames', model_params_ne.Properties.VariableNames);
    model_params = [model_params_ne; new_params_table];
    [~,ix] = unique(model_params.param_name,'stable','last');
    model_params = model_params(ix,:);
    
    % Get the locations (i.e. indices) of model_params_e in the complete parameters set, i.e., in model_params
    indices_e = nan(height(model_params_e),2);
    for param_num = 1:height(model_params_e)
        current_param_name = model_params_e.param_name{param_num};
        current_statistic = model_params_e.statistic{param_num};
        req_row_number = find(strcmp(model_params.param_name,current_param_name));
        if strcmp(current_statistic,'mu')
            indices_e(param_num,:) = [req_row_number,2];
        elseif strcmp(current_statistic,'sig')
            indices_e(param_num,:) = [req_row_number,3];
        end
    end
end

function model_params = update_model_params_from_e(model_params_table, params_values_e, indices_e)
    % update_model_params_from_e
    % 
    % INPUTS:
    %   model_params_table : The 30x3 table containing param_name, mu, sig
    %   params_values_e    : Array of new estimated values
    %   indices_e          : Nx2 matrix of [row, col] indices mapping values to the table
    %
    % OUTPUT:
    %   model_params      : The struct with format params.varname.mu / params.varname.sig
    
    % Update the Table values first
    % We use curly brace indexing table{row, col} to assign values
    for i = 1:length(params_values_e)
        r = indices_e(i,1);
        c = indices_e(i,2);
        model_params_table{r, c} = params_values_e(i);
    end
    
    % Convert the updated Table into the "param_vary" struct format
    model_params = struct();
    
    for i = 1:height(model_params_table)
        p_name = model_params_table.param_name{i};
        p_mu   = model_params_table.mu(i);
        p_sig  = model_params_table.sig(i);
        
        % Create the nested structure: params.Name.mu and params.Name.sig
        model_params.(p_name).mu = p_mu;
        model_params.(p_name).sig = p_sig;
    end
end
 