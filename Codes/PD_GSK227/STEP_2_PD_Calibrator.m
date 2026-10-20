clear; 
warning off; parfevalOnAll(@warning,0,'off','all');
load("New Analysis-sept/Lesion_Waterfall_8mpk_all.mat"); 
load("New Analysis-sept/PSA_Waterfall_8mpk_all.mat"); 
load("New Analysis-sept/Mean_spider_clinical_8mpk_all.mat");
load('New Analysis-sept/PFS_Table_8mpk_all.mat');

Sim_Params = struct(); % struct to pass all crucial choice related parameters into the simulating function

%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Automated_calibration_switch = false;                  % If true, the calibration will be done by the solver, you will need to provide the parameter initial values and their bounds.
Calibrating_what = 'PSA';                           % Valid inputs are 'LESION' or 'PSA'.
Sim_Params.UseParallel_simfcn = ~Automated_calibration_switch;
Sim_Params.Solvertype = 'global';                      % 'local' or 'global' are the valid inputs. 

% Load Project model
proj = sbioloadproject("../../GSK227_CRPC.sbproj");    % This path has to be changed as per the location of the model
model = proj.m1; 

Sim_Params.seed = 0;                                  % This is the RNG Seed value to which you are fitting
Sim_Params.N_VP = height(Lesion_Waterfall.BOR);        % Number of patients to simulate is same as that in the data
Sim_Params.Dose = PSA_Waterfall.Dose_mpk;       % in mg/kg , should be matching with the dose used in Step_1_Data_processor script
Sim_Params.variability = 1;                             
Sim_Params.StopTime = 340;                              % Simulating for 1 year
Sim_Params.Time_Sim = 0:1:Sim_Params.StopTime;          % for each of the days, time point where the simulations will be carried out
Sim_Params.Dose_Cycle_length = 21;                      % in days, 21 Days for Q3W
Sim_Params.Doses_days_in_a_cycle = 0;                   % On what days you want to administer the dosage, it could be [0,7] if you want d1,d8 Q3W kind of dosing
Sim_Params.Dose_target_name = 'Central_ADC.Cen_ADC_nM'; % Species name in the model where ADC is dosed
Sim_Params.Dar_Dose_target_name = 'Central_ADC.Cen_DummyDAR'; % Species name in the model where DAR is dosed
Sim_Params.Scan_Interval = 42;                          % Time interval between two radiographical scan in simulation, (in days)
Sim_Params.Time_efficacy=0:Sim_Params.Scan_Interval:Sim_Params.StopTime; % response measument days 



Sim_Params.StatesToLog = {'TV','PSA'};
Sim_Params.parameter_name_for_PSA = 'PSA';
Sim_Params.parameter_name_for_tumor_volume = 'TV';
Sim_Params.parameter_name_for_hscore = 'H_Score';
Sim_Params.parameter_name_for_bodyweight = 'W';
Sim_Params.parameter_name_for_mol_weight_ADC = 'MW_ADC';
Sim_Params.N_tmr = 1;
Sim_Params.Param_names_for_intrapatient_heterogeneity = {'Tumor_Dia','H_Score','K1','Kmax_p'}; % Parameters that will change across within each patient
Sim_Params.Intrapatient_Heterogeneity = 0.0;               % Variability in parameters across different tumors within patient
Sim_Params.dropout_Switch = true;                          % Switch to turn off/on the dropouts due to New metastats or Non-specific reasons
Sim_Params.NM_Params = [0.028404, 0.14067, 56.603, 1.2763];% New-metastates related dropout probability modelling parameters
Sim_Params.BLDrop_Params = [0.18 15  0.19 3.5];            % Drop out parameters related to bone lesions/metastases
Sim_Params.NonSpecific_dropouts = 0.06245;                 % Probability for dropout related to non-radiographical reasons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-exploratory params, that are different from model values
model_params_drug_specific = {
    %param_name    %mu        %sig
    'MW_ADC'       148000     0;
    'MW_PL'        533.55        0;
    'DAR_t0'       4          0;
    'Kd_ADC'       35.4        0.1;
    'Kon_ADC'      5.2       0.12;
    };

model_params_indication_specific = {
    %param_name            %mu        %sig
    'Tumor_Dia'            49         1.2;
    'H_Score'              225        0.2;
    'Free_fraction_PL'     0.325707     0;
    'Ag_p'                 220000     0.2;
    'Ag_n'                 0          0.0;
    'Single_Cell_Vol_p'    0.27e-12   0.0;
    'K_Deg_Lyso'           4.2        0;
    'Kin_PL_p'             2982       0;
    'Kin_PL_n'             2982       0;
    'Kint_ADC'             5.9047     0;
    'KC50_p'               8.78       0.15;
    'K1'                   0.0498     0.081215;
    'Kmax_p'               0.725674   0.075566;
    'Cp'                   0.057523  0.498943;
    'CL_PSA'               0.083857  1.360775;
    'Initial_PSA_ngml'     76         1.1;
    'K_res'                2.997727  0.860893;
    'saturation_exp'       1/3          0;
    'Kmax_multiplier'      0.275657   0.392483;
    };

model_params_patient_specific = {
    %param_name    %mu         %sig
    'V1_ADC'       3.554,       0.1578;
    'V2_ADC'       1.567,       0.3326;
    'CL0_ADC'      24*0.03305,     0.216;
    'CLD_ADC'      24*0.0331,      0.10;
    'V1_PL'        54.53,       0.4251;
    'V2_PL'        1787,     0.10;
    'CL_PL'        24*16.62,       0.3647;
    'CLD_PL'       24*1.703,       1.12;
    'alpha_ADC'    0.6611,        0.10;
    'beta_ADC'     0.033744,     0.10;
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
    % 'CL_PSA'          'mu'           0.083857            1e-2      0.5;
    % 'Cp'              'mu'           0.057523         0.01      0.2;
    % 'Cp'              'sig'          0.498943            0.4       0.5;  
    % 'CL_PSA'          'sig'          1.360775            0.5       2;
    % 'K_res'           'mu'           2.997727            2         3;
    % 'K_res'           'sig'          0.860893            0.5       1;
    % 'Kmax_multiplier'   'mu'           0.2757           0.25    0.75;
    % 'Kmax_multiplier'   'sig'          0.3925           0.05    0.5;
    % 'Kmax_p'            'mu'           0.7257           0.01    1;
    % 'Kmax_p'            'sig'          0.07557          0.05    0.5;
    % 'Free_fraction_PL'  'mu'           0.3257             0.1     10;
    % 'K1'                'sig'          0.08122          0.05     0.5;   
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
[model_params, indices_e] = combine_params(model_params_e,model_params_ne,model);
initial_values = model_params_e.initial_value;    

% Configure SimBiology Model
csObj = getconfigset(model);
csObj.SolverType = 'sundials';
csObj.SolverOptions.RelativeTolerance = 1e-4;
csObj.SolverOptions.AbsoluteTolerance = 1e-6;
csObj.SolverOptions.MaxStep = 1;
csObj.SolverOptions.AbsoluteToleranceScaling = false;
csObj.StopTime = Sim_Params.StopTime;
csObj.TimeUnits = 'day';
csObj.SolverOptions.OutputTimes = Sim_Params.Time_Sim;
csObj.RuntimeOptions.StatesToLog = Sim_Params.StatesToLog;    

% Create simfunction
param_names =  model_params{:,1};
observables_names = Sim_Params.StatesToLog;
num_dose_days = length(Sim_Params.Doses_days_in_a_cycle);
doses_targets = [repmat({Sim_Params.Dose_target_name},[1,num_dose_days]),...
                 repmat({Sim_Params.Dar_Dose_target_name},[1,num_dose_days])];
F_simfcn = createSimFunction(model,param_names,observables_names,doses_targets,'UseParallel',Sim_Params.UseParallel_simfcn);

if Automated_calibration_switch
    
    lower_bound    = model_params_e.lower_bound; %#ok<*UNRCH>
    upper_bound    = model_params_e.upper_bound;

    objFun=@(fit_Params) objfun(fit_Params, model_params,indices_e, Sim_Params, F_simfcn, Lesion_Waterfall, PSA_Waterfall, Calibrating_what);
    fprintf('Initiating automated Calibration with %d parameters...\n\n', length(initial_values));
    

    if strcmp(Sim_Params.Solvertype, 'global')

        % ga options
        safeObjFun = @(x) real(objFun(x));
        nvars = numel(initial_values);
    
        %Setup GA options with the LOADED population
        % r=load('output_gen_30.mat');
        gaopts = optimoptions('ga', ...  
            'PopulationSize', 160, ... 
            'MaxGenerations', 30, ...           
            'Display', 'iter', ...              
            'UseParallel', true, ...            
            'PlotFcn', {@gaplotbestf, @gaplotscores, @gaplotdistance}, ... 
            'OutputFcn', @monitor_ga);%, ...
            % 'InitialPopulationMatrix', r.current_pop); % CRITICAL: This connects the old data
                      
        
        % Execute GA
        fprintf('Initiating GA Calibration ...\n');
        tic;
        
        % Running the GA
        [bestParams, bestError, exitflag, output, population] = ...
            ga(safeObjFun, nvars, [], [], [], [], lower_bound, upper_bound, [], gaopts);
        
        time_taken = toc;
        fprintf('==============================================\n');
        fprintf('Took %.1f seconds to calibrate.\n', time_taken);
        fprintf('==============================================\n');

    elseif strcmp(Sim_Params.Solvertype, 'local')
        fminconopts=optimoptions('fmincon','Display', 'iter-detailed', 'MaxIterations', 2, 'Algorithm', 'interior-point', 'UseParallel',true,'PlotFcn','optimplotfval');
        tic;
        [bestParams,bestError, exitflag, output] = fmincon(objFun, initial_values, [], [], [], [], lower_bound, upper_bound, [], fminconopts);
        time_taken = toc;
        
        fprintf('==============================================\n');
        fprintf('Took %.1f seconds to calibrate.\n', time_taken);
        fprintf('==============================================\n');

    else
        error('Check the definition of Sim_Params.Solvertype. It should either be "local" or "gloabl".');
    end
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
Sim_Params.UseParallel_simfcn = true;
Simulation_Results = run_population_simulation( Sim_Params, F_simfcn, true);
time_taken = toc;
fprintf('Took %.1f seconds to simulate.\n', time_taken);


plot_waterfall(Simulation_Results, Lesion_Waterfall, PSA_Waterfall); 
plot_pfs(Simulation_Results, PFS_Table);                             % [OPTIONAL] Lines marked as optional can be commented without any issue
plot_spider(Simulation_Results, Sim_Params, median(Lesion_Waterfall.TTR, 'omitnan'));              % [OPTIONAL]
plot_psa_spider(Simulation_Results, Sim_Params);            % [OPTIONAL]
% plot_swimmer(Simulation_Results);                         % [OPTIONAL]
% plot_observables(Simulation_Results, Sim_Params)          % [OPTIONAL] 
% summarize_patient_response(Simulation_Results);           % [OPTIONAL]
% Mean_spider_simulated = calculate_mean_resp_at_clinical_visits(Simulation_Results.AllPlots.Sim_SLD(:,[1,2,4,12]));
% plot_mean_response_at_visits(Mean_spider_clinical, Mean_spider_simulated);

































%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective Function 

function err = objfun(fit_Params, model_params,indices_e, Sim_Params, F_simfcn, Lesion_Waterfall, PSA_Waterfall, Calibrating_what) 
    %
    % Defines the objective function that will be minimized for automated calibration process.
    %
    % INPUTS:
    %       fit_Params:         The value for the parameters (to be estimated) at which they will be evaluated in an iteration
    %       model_params:       Set of values for all the parameters to be set and simulated in the model
    %       indices_e:          Indices of the parameters(to be estimated) in the model_params, so that they can replaced
    %       model:              A SimBiology model object (to be simulated)
    %       Sim_Params:         Information about the crucial parameters to use while simulation
    %       Lesion_Waterfall:    Processed data containing information about the lesion waterfall plot from the raw data
    %       PSA_Waterfall:       Processed data containing information about the PSA waterfall plot from the raw data
    %       Calibrating_what:   'PSA' or 'LESION'. Will use the corresponding waterfalls to calibrate
    %       
    %
    % OUTPUTS:
    %       err:                The square root of sum of square of differences between the Clinical waterfall and the simulated waterfall
    %                 
    
    
    rng(Sim_Params.seed, "twister");    % Set RNG seed for reproducibility

    % Pass the parameters
    Sim_Params.param = update_model_params_from_e(model_params,fit_Params,indices_e);
 
    Simulation_Results = run_population_simulation(Sim_Params, F_simfcn,  false); % Run simulation 
    simulated_waterfall_Data = Simulation_Results.AllPlots.Sim_Waterfall;      % Extract the Simulated Waterfall data
    simulated_waterfall_Data = sortrows(simulated_waterfall_Data, -2);

    if strcmp(Calibrating_what,'PSA')
        psawaterfall = min(100,simulated_waterfall_Data(end-height(PSA_Waterfall)+1:end,4));
        y = sort(psawaterfall, 'descend');
        % y = simulated_waterfall_Data(:,2); 
        obs=sort(PSA_Waterfall.BOR,'descend');     % The clinical waterfall data
    elseif strcmp(Calibrating_what, 'LESION')
        obs=sort(Lesion_Waterfall.BOR,'descend');     % The clinical waterfall data
        y = simulated_waterfall_Data(:,2);
    else
        error('Wrong input for the parameter "Calibrating_what" in objfun. Use either "PSA" or "LESION".')
    end

    Diff= y-obs;
    err=sqrt(mean(Diff.^2)); % Calculating the error metric to minimise
end
 
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
        if evaluate_psa50_response(subj) %&& (best_change_ORR < 0) % PSA50 response flag
            n_PSA50 = n_PSA50 + 1;
        end
    end
    ORR = n_ORR/nPatients;
    PSA50 = n_PSA50; % count of patients
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
 
function plot_waterfall(Simulation_Results, Lesion_Waterfall, PSA_Waterfall)
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
    psawaterfall = simulated_waterfall_Data(max(1,end-height(PSA_Waterfall)+1):end,4);
    psawaterfall = sort(psawaterfall, 'descend');
    psawaterfall = min(100, psawaterfall);
    figure;
    hold on;    
    % Create the bar plot for the Data
    data_prc_change=sort(Lesion_Waterfall.BOR,'descend');
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
    % Mentioning ORR on the plot
    text_loc_x = max(N_VP - 1);   
    text_loc_y = y_max - 5; 
    text_str = sprintf('Clinical ORR is: %.1f%%', sum(Lesion_Waterfall.OVLRSPN < 2)/height(Lesion_Waterfall) * 100);
    text(text_loc_x, text_loc_y, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');
    
    
    text_str = sprintf('Simulated ORR is: %.1f%%', 100*Simulation_Results.AllPlots.ORR);
    text(text_loc_x, text_loc_y - 10, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');

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
    bar(data_psa_change, 'FaceColor', gem_hex_codes(1), 'BarWidth', 1, 'EdgeColor',gem_hex_codes(1), 'DisplayName','Data');
    
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
    xlim([0.5,height(data_psa_change)+0.5]);
    text_loc_x = height(data_psa_change);   
    text_loc_y = y_max - 15; 
    text_str = sprintf('Clinical PSA50 is: %.1f%%', sum(PSA_Waterfall.BOR<=-50)/height(PSA_Waterfall) * 100);
    text(text_loc_x, text_loc_y, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');
    
    
    text_str = sprintf('Simulated PSA50 is: %.1f%%', 100*Simulation_Results.AllPlots.PSA50/height(PSA_Waterfall));
    text(text_loc_x, text_loc_y - 20, text_str, 'FontSize', 9, 'FontWeight', 'normal', 'Color', 'k', 'Margin', 2,  'HorizontalAlignment', 'right');

    grid off; box on; set(gca,'LineWidth',1.25);
    pbaspect([3 2 1]);
    fontsize(15,'points');
    legend('Location','southwest','Box','on','FontSize',15);
    title('Calibrating PSA Waterfall Plot','FontWeight','bold','FontSize',18);
    fontname('Poppins');
    hold off;
end
 
function plot_pfs(Simulation_Results, PFS_Table)
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
    
    time = PFS_Table.PFSD;
    censor = PFS_Table.PFSCNSRN; % 0 = PD (Event), 1 = Censored
    
    % Calculate Kaplan-Meier survival curve & 95% CIs using ecdf
    [f, x, ~, ~] = ecdf(time, 'Censoring', censor, 'Function', 'survivor', 'Alpha', 0.05);   
    
    % Plot PFS Kaplan-Meier Curve
    stairs(x, f*100, 'k-', 'LineWidth', 2, 'DisplayName', 'Clinical PFS');
    
    % Mark censored observations with tick marks (+)
    cens_idx = (censor == 1);
    if any(cens_idx)
        cens_times = time(cens_idx);
        [x_unique, u_idx] = unique(x, 'last');
        cens_f = interp1(x_unique, f(u_idx), cens_times, 'previous');
        plot(cens_times, 100*cens_f, 'r+', 'LineWidth', 2, 'MarkerSize', 8, 'HandleVisibility','off');
    end
    

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
    median_TTR_Sim = calculate_TTR(Sim_SLD(:,[1, 2, 4])); % this subset of Sim_SLD contains the patient_id, Time, % change in SLD
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
    
        data_cut(:,6) =  min(100,data_cut(:,6));
        % Find this patient's IHC category
        IHC_class = IHC_Category(pid); 
        % if min(data_cut(2:end,4))<=0
        plot(data_cut(:,2), data_cut(:,6), '-', 'Color', 'k', 'LineWidth', 2);
        % end
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
    xlim([0, Sim_Params.StopTime]);
    grid off; box on; set(gca, 'LineWidth',1.25);
    pbaspect([ 3 2 1]);
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
    
    
    Time_Sim = Sim_Params.Time_Sim;%(1:min(length(Sim_Params.Time_Sim), Sim_Params.Scan_Interval+1));
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
        working_ts = working_ts(:,1:length(Sim_Params.Time_Sim));
        
        % Plot multiple simulated trajectory (N_VP curves)
        plot(Time_Sim, working_ts, 'Color', gem_hex_codes(i), 'LineWidth', 1); %if this line gives error due to gem_hex_codes, just remove the color property from here. It gives error because you are logging more than 7 states.
        
        title(sprintf('Plotting Time Series %d',i-1));

        xlabel('Time (days)'); ylabel(Sim_Params.StatesToLog{i},'Interpreter','none');
        xlim([-1 Sim_Params.StopTime]);
        % ylim([0.1*min(working_ts(:)), 10*max(working_ts(:))]);
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
    AllPlots = Simulation_Results.AllPlots;
    VP_Init_Array = Simulation_Results.VP_Init_Array;
    Sim_SLD      = AllPlots.Sim_SLD;
    DropoutTime  = AllPlots.DropoutTime;
        
    PatientList = unique(Sim_SLD(:,1));
    nPatients   = numel(PatientList);
    PatientSummary = zeros(nPatients,4);  % Columns will be [PID, Dose, Hscore, BestResp]

    for i = 1:nPatients
        pid = PatientList(i);
        Dose = Simulation_Results.Dose(pid);        
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

function [model_params, indices_e] = combine_params(model_params_e,model_params_ne,model)
      
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
    model_params = sortrows(model_params); % To ensure alphabetical ordering.

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
  
function median_TTR_Sim = calculate_TTR(Sim_SLD_subset)
    
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
 
function Mean_response_table = calculate_mean_resp_at_clinical_visits(Sim_SLD_subset)
    
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
    
    for k = 1:length(unique_visits)
        vis = unique_visits(k);
        responses = Clean_Sim_SLD(Clean_Sim_SLD(:, 5) == vis, 3); 
        
        mean_resp_list(k) = mean(responses, 'omitnan');
    end
    
    Mean_response_table = table(unique_visits, mean_resp_list, 'VariableNames', {'Sim_visits', 'Mean_Resp_Sim'});
end
 
function plot_mean_response_at_visits(Clinical_table, Simulated_table)
    figure; hold on;
    box on; grid off; set(gca, 'LineWidth', 1.25);
    errorbar(Clinical_table.Clinical_visits, Clinical_table.Mean_Resp_Clinical, ...
        Clinical_table.Mean_Resp_Clinical - Clinical_table.Min_Resp_Clinical, ...
        Clinical_table.Max_Resp_Clinical - Clinical_table.Mean_Resp_Clinical, ...
        'o','LineWidth',2,'MarkerFaceColor','auto','DisplayName','Clinical(Range)');
    % plot(Clinical_table.Clinical_visits, Clinical_table.Mean_Resp_Clinical, 'DisplayName','Clinical');
    plot(Simulated_table.Sim_visits, Simulated_table.Mean_Resp_Sim, 'DisplayName','Simulated');
    xlabel('Patient Visit No.');
    ylabel('Mean response in SLD (%)');
    legend('Location','southwest','Box','on');
    ylim([-105, 105])
    fontsize(15,'points');
    title('Comparing Mean response for clinical and simulated data','FontWeight','bold','FontSize',18);
    fontname('Poppins');
    hold off;
end
 
% Function to monitor GA optimizer params at each generation
function [state, options, optchanged] = monitor_ga(options, state, flag)
    optchanged = false;
    
    switch flag
        case 'iter'
            % 
            [minScore, idx] = min(state.Score);
            best_p = state.Population(idx, :);
            current_score = state.Score;
            current_pop = state.Population;
            current_gen = state.Generation;
            
            % Print Summary to Command Window
            fprintf('\n>>> Generation %d Completed\n', current_gen);
            fprintf('Current Best Error: %.6f\n', minScore);
            
            % Print Best Parameters
            fprintf('Best Parameters: ');
            fprintf('%0.6f  ', best_p); 
            fprintf('\n'); 
            % ------------------------------------
            
            % INCREMENTAL SAVE
            fname = sprintf('output_gen_%d.mat', current_gen);
            save(fname, 'best_p', 'minScore', 'current_pop', 'current_gen',"current_score");
            
            
            fprintf('Progress saved to %s\n', fname);
            fprintf('--------------------------------------\n');
    end
end
 