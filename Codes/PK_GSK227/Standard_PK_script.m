% This script is for PK fitting and PK simulation both.
% There are total 3 ways to run this script using the option
% User_Inputs.run_type. In case User_Inputs.run_type = 'auto_fitting', the
% script will do automated fitting using a MATLAB in-built optimizer and plot the single mean curve. In
% case User_Inputs.run_type = 'manual_fitting', the script will simply run
% single patient simulation with the given params and give single fit plots. In case
% User_Inputs.run_type = 'simulation', the script will use the user given
% parameters and run VPoP and plot the results.

clear; warning off;
%% USER INPUTS

User_Inputs = struct(); % This struct contains all the user inputs that go into functions

% Set the seed for reproducibility
User_Inputs.seed = 0; % This seed is used for both calibration as well as for optimizer

% model names
User_Inputs.projfile = "../../GSK227_CRPC.sbproj";
User_Inputs.model_id = 1; % access the first model in the project file
User_Inputs.StatesToLog = {'Cen_ADC_ugml', 'Cen_tAb_ugml', 'Cen_PL_ugml'};
User_Inputs.ADC_Model_Species = 'Cen_ADC_ugml';
User_Inputs.TAB_Model_Species = 'Cen_tAb_ugml';
User_Inputs.PL_Model_Species = 'Cen_PL_ugml';

% PK Data file
User_Inputs.PK_file = "PK_Data_for_Calibration.xlsx"; % For calibrations
User_Inputs.PK_Summary_file = 'PK_Summary.xlsx'; % For visualization


% Specify columns
User_Inputs.Time_column_name = 'Time';
User_Inputs.Dose_column_name = 'Dose_mpk';
User_Inputs.Group_column_name = 'DoseID';
User_Inputs.Observables_column_names = {'Mean_SERUM_GSK5757810_CONCENTRATION', 'Mean_SERUM_GSK5764227_CONJUGATED_CONCENTRATION','Mean_SERUM_GSK5764227_TOTAL_CONCENTRATION'};
User_Inputs.groupdata_var_names = {'DoseID','Dose_mpk','Time','PL_ugml','ADC_ugml','tAb_ugml','Dose_nmole'};
User_Inputs.groupdata_var_units = {'','milligram/kilogram','hour','microgram/milliliter','microgram/milliliter','microgram/milliliter','nanomole'};

% Response Map
User_Inputs.responseMap =  {'Cen_ADC_ugml = ADC_ugml','Cen_tAb_ugml = tAb_ugml','Cen_PL_ugml = PL_ugml'};



% Give user the choice of type of run
User_Inputs.run_type = 'manual'; % 'auto_fitting' or 'manual_fitting' or 'simulation'


% ------ Common Options for both fitting and Simulation -----
% Parameters that are not estimated
Non_Exloratory_Params = {
     % Fixed Parameter Name,          mu,           sigma;    
     'W',                             70,         0.0;
     'DAR_t0',                        4,            0.0;
     'MW_ADC'                         151000,         0;
     'MW_PL'                           533.55,         0;
     'Switch_TumorComp',              0,            0.0;
     'Switch_TMDD',                   1,            0.0;
     'Switch_Tox',                    1,            0.0;
     'V1_ADC'       3.554,       0.1578;
    'V2_ADC'       1.567,       0.3326;
    'CL0_ADC'      0.03305,     0.216;
    'CLD_ADC'      0.0331,      0.10;
    'V1_PL'        54.11,       0.4231;
    'V2_PL'        1750,     0.10;
    'CL_PL'        16.59,       0.3675;
    'CLD_PL'       1.714,       1.148;
    'alpha_ADC'    0.6611,        0.10;
    'beta_ADC'     0.033744,     0.10;
     'K_P_dec',                       0,            0.00;
};
User_Inputs.Non_Exploratory_Params_Table = cell2table(Non_Exloratory_Params, ...
    'VariableNames', {'Name','mu', 'sigma'});

% Parameters that are estimated or explored
Exploratory_Params = {
    %'Parameter', 'init_guess',    'lower_lim', 'upper_lim',      'Statistic'
    % 'V1_ADC',      0.054073825,     0.001,       0.5,              'mu';
    % 'V2_ADC',      0.009850577,     0.001,       0.5,              'mu';
    % 'CL_ADC',      0.01,     0.001,       0.5,              'mu';
    % 'CLD_ADC',     0.007991744,     0.0005,      0.01,             'mu';
    % 'V1_PL',       0.070672592,     0.04,        0.1,              'mu';
    % 'V2_PL',       0.035035142,     0.03,        0.1,              'mu';
    % 'CL_PL',       4.008860814,     4.6,         6.0,              'mu';
    % 'CLD_PL',      0.117373022,     0.21,        0.43,             'mu';
    % 'K_P_dec',     0.010502780,     0.007,       0.012,            'mu';
     
    % Below parameters are to be uncommented when TMDD is turned on
    'Kd_ADC',      35.4,         32.0,        36.0,            'mu';
    %'Kon_ADC',     5.2,          4.0,         6.0,             'mu';
    %'Initial_PT',  2.0,          1.0,         3.0,             'mu';
    %'K_syn_PT',    10.0          8.0,         11.0,            'mu';
    };

User_Inputs.Exploratory_Params_Table = cell2table(Exploratory_Params, ...
    'VariableNames', {'Name','InitialValue','LowerBound','UpperBound', 'Statistic'});


% ------------------------------ Options for auto_fitting ----------------------

User_Inputs.fit_pooled = false; % true for pooled and false for unpooled
User_Inputs.plot_separate_unpooled = true; % Keep this true to plot with unpooled parameter set for different doses, false for plotting with mean params
User_Inputs.solvername = 'ga'; % lsqnonlin, fminunc, lsqcurvefit, fminsearch, fmincon, scattersearch, patternsearch, ga, particleswarm
User_Inputs.errormodel = 'exponential'; % constant, proportional, combined, exponential


% ------------------------------Options for Simulation ----------------------

User_Inputs.N_VP = 1;                             % Number of virtual patients
User_Inputs.Sampling_time = 0:1:1600;                % Time grid for PK profile (hours)



% --------------------------- Dose Related Inputs --------------------------
User_Inputs.Doses_days_in_a_cycle = 0;
User_Inputs.Dose_Cycle_length = 504; % hours
User_Inputs.Dose_Target = 'Central_ADC.Cen_ADC_nM'; % Species where the drug is injected
User_Inputs.DAR_Dose_Target = 'Central_ADC.Cen_DummyDAR'; % For DAR dosing

%% MAIN SCRIPT EXECUTION


% Throw Error in case Exploratory_Params_Table contains non-mu statistic

Check_mu_Present(User_Inputs);


% Load and Prepare Data for Fitting

[gData, Dose_mpk_list, estimatedParam] = loadAndPrepareData(User_Inputs);


% Load and Configure Model for Fitting

[model, csObj_fit, dose2s] = loadAndConfigModelForFit(User_Inputs, gData);


% Get the combined parameter table before fitting

model_params_in = combine_params(User_Inputs);

% In case the run_type is auto_fitting, use the solver to fit, else use the
%  parameters in Exploratory_Params_Table for the simulation and subsequent
%  Plots

if strcmp(User_Inputs.run_type, 'auto_fitting')
    % Define Responses for Fitting
    
    % Get Solver Options
    opts = getSolverOptions(User_Inputs);
    
    % Run Parameter Estimation (Fitting)
    FitResults = runFitting(model, gData, estimatedParam, dose2s, opts, User_Inputs);

    % Process Fitting Results and save data to an excel file
    [model_params_out, allEstimates] = processFitResults(FitResults, User_Inputs, Dose_mpk_list);
    
    % Use the output parameters for simulation in case of auto_fitting
    model_params_sim = model_params_out;

else
    % Case: simulation or manual_fitting
    % Use the input parameters for simulation
    model_params_sim = model_params_in;
    allEstimates = []; % Not applicable for this path

end


% Run Virtual Patient Simulation and Plot Results. In case of auto_fitting uses
% the fitted params and in case of simulation, uses the input params from
% Exploratory_Params_Table and Non_Exploratory_Params_Table
runSimulationAndPlot(model_params_sim, allEstimates, Dose_mpk_list, User_Inputs);

%% FUNCTION DEFINITIONS

% Function to access and load the cleaned PK Data for fitting

function [gData, Dose_mpk_list, estimatedParam] = loadAndPrepareData(User_Inputs)
    % loadAndPrepareData : Loads the data for fitting
     %
     % Inputs:
     %   PK_Data_for_Calibration.xlsx           - File created from
     %   previous code that has PK data to be used for calibrations.
     %   User_Inputs.Exploratory_Params_Table   - Table that contains
     %   parameters to be estimated/explored
     %   User_Inputs.Non_Exploratory_Params_Table   - Table that contains
     %   parameters that are not estimated/explored
     %
     % Outputs:
     %  gData                   - Grouped Data that is fed into the
     %  optimiser for the parameter estimation
     %  Dose_mpk_list           - List of all doses (in milligram/kilogram)
     %  present the PK data
     %  estimatedParam          - Formats user input for the exploratory
     %  params into a MATLAB readable format that optimiser can read




PK_file = User_Inputs.PK_file;
Exploratory_Params_Table = User_Inputs.Exploratory_Params_Table;
Non_Exploratory_Params_Table = User_Inputs.Non_Exploratory_Params_Table;
    BodyWeight = Non_Exploratory_Params_Table.mu(strcmp(Non_Exploratory_Params_Table.Name, 'W'));
    MW_ADC = Non_Exploratory_Params_Table.mu(strcmp(Non_Exploratory_Params_Table.Name, 'MW_ADC'));
    % Create estimatedInfo array for fitting
    estimatedParam = estimatedInfo(Exploratory_Params_Table.Name, ...
        'InitialValue', Exploratory_Params_Table.InitialValue, ...
        'Bounds', [Exploratory_Params_Table.LowerBound Exploratory_Params_Table.UpperBound]);
    
    % Load the data
    format long g
    PK_Data_for_Calibration=readtable(PK_file);
    PK_Datai=PK_Data_for_Calibration(:,:);
    Cols_to_extract = [User_Inputs.Group_column_name, User_Inputs.Dose_column_name, ...
        User_Inputs.Time_column_name, User_Inputs.Observables_column_names];
    %PK_Datai = PK_Datai(:, [1 2 3 4 7 10]);
    PK_Datai = PK_Datai(:, Cols_to_extract);
    PK_Datai = PK_Datai(PK_Datai.(User_Inputs.Time_column_name) <= 504, :); % use this line of you want to consider only the first cycle
    PK_Datai.Dose_mpk(PK_Datai.(User_Inputs.Time_column_name) == 504) = NaN;
    PK_Datai.Dose_nmol = PK_Datai.(User_Inputs.Dose_column_name) * BodyWeight / 1000 / MW_ADC * 1e9;
    
    % Zero to NaN Handling
    colsToReplace = User_Inputs.Observables_column_names;
    for c = colsToReplace
        colData = PK_Datai{:, c};
        colData(colData == 0) = NaN;
        PK_Datai{:, c} = colData;
    end
    gData=groupedData(PK_Datai);
    gData.Properties.VariableNames = User_Inputs.groupdata_var_names;
    gData.Properties.VariableUnits = User_Inputs.groupdata_var_units;
    
    % Set gData properties
    gData.Properties.IndependentVariableName = User_Inputs.Time_column_name;
    gData.Properties.GroupVariableName = User_Inputs.Group_column_name;
    
    % Get Dose_mpk_list (needed for processing results)
    Dose_mpk_list = (PK_Datai.Dose_mpk(~isnan(PK_Datai.Dose_mpk)))';
end

% Function to load and configure the model for auto_fitting

function [model, csObj, dose2s] = loadAndConfigModelForFit(User_Inputs, gData)

 % loadAndConfigModelForFit : Loads and configures the model for fitting
     %
     % Inputs:
     %      User_Inputs.projfile           - The name of the Simbiology project
     %      User_Inputs.model_id           - Model index within Simbiology
     %      project.
     %      gData                          - Grouped data obtained from
     %      loadAndPrepareData function
     %
     % Outputs:
     %      model                           - Creates a model object that
     %      script can use for fitting
     %      csObj                            - Creates an object that stores the
     %      simulation settings
     %      dose2s                           - The dose object that is used
     %      for fitting

    projfile = User_Inputs.projfile;
    model_id = User_Inputs.model_id;
    % Load Model
    projStruct = sbioloadproject(projfile);
    model_var = "m" + string(model_id);   % creates "m1"
    model = projStruct.(model_var);
    csObj = getconfigset(model);
    csObj.TimeUnits = 'hour';
    csObj.SolverType = 'ode15s';
    csObj.SolverOptions.RelativeTolerance = 1e-4;
    csObj.SolverOptions.AbsoluteTolerance = 1e-6;
    csObj.SolverOptions.MaxStep = 0.5;
    csObj.SolverOptions.AbsoluteToleranceScaling = false;
    csObj.CompileOptions.UnitConversion = true;
    % DAR dose
    dose1 = sbiodose('dose1','TargetName',User_Inputs.DAR_Dose_Target);
    dose1.StartTime   = 0;
    dose1.Amount      = 1;
    dose1.Rate        = 0;
    dose1.AmountUnits = 'mole';
    dose1.TimeUnits   = 'hour';
    dose1.Interval = User_Inputs.Dose_Cycle_length;
    dose1.RepeatCount = round(User_Inputs.Sampling_time(end)/User_Inputs.Dose_Cycle_length) + 1;
    
    % Copy this object to all the dose days in a cycle
    DAR_doses = [];
    for dose_day = User_Inputs.Doses_days_in_a_cycle
        DAR_dose = copyobj(dose1);
        DAR_dose.StartTime = dose_day;
        DAR_doses = [DAR_doses,DAR_dose]; %#ok<AGROW>
    end
    dose1 = DAR_doses; clear DAR_doses DAR_dose %#ok<NASGU>

    % Dose from Data for fitting
    dose2 = sbiodose("Dose_nmole");
    dose2.TargetName = User_Inputs.Dose_Target;
    dose2s = createDoses(gData,"Dose_nmole","",dose2);
end

% Function having the optimiser settings

function opts = getSolverOptions(User_Inputs)
% getSolverOptions : Gives user the choice of different solvers
     %
     % Inputs:
     %     solvername                       - The name of the optimizer to be
     %     used for optimization
     %
     % Outputs:
     %      opts                            - The options used for the
     %      optimization

    solvername = User_Inputs.solvername;
    switch lower(solvername)
        % Local Solvers
        case 'lsqnonlin'
            disp('>>> Running LSQNONLIN optimization...');
            opts = optimoptions('lsqnonlin', ...
                    'Display', 'iter', ...           % Show iteration output
                    'MaxIterations', 1000, ...        % Maximum iterations
                    'MaxFunctionEvaluations', 5000,  ...
                    'StepTolerance', 1e-6, ...
                    'FunctionTolerance', 1e-6, ...
                    'OptimalityTolerance', 1e-6);
        case 'fminunc'
            disp('>>> Running FMINUNC optimization...');
            opts = optimoptions('fminunc', ...
                'Algorithm','quasi-newton','Display','iter', ...
                'MaxIterations',1000,'FunctionTolerance',1e-8);
        case 'lsqcurvefit'
            disp('>>> Running LSQCURVEFIT optimization...');
            opts = optimoptions('lsqcurvefit','Display','iter', ...
                'MaxIterations',1000,'MaxFunctionEvaluations',5000);
        case 'fminsearch'
            disp('>>> Running FMINSEARCH optimization...');
            opts = optimset('Display','iter','MaxIter',1000,'TolX',1e-8,'TolFun',1e-8);
        case 'fmincon'
            disp('>>> Running FMINCON optimization...');
            opts = optimoptions('fmincon', ...
                'Algorithm','interior-point','Display','iter', ...
                'MaxIterations',1000,'FunctionTolerance',1e-8);
        % Global Solvers
        case 'scattersearch'
            disp('>>> Running SCATTERSEARCH optimization (patternsearch variant)...');
            opts = optimoptions('patternsearch', ...
                'Display','iter','UseCompletePoll',true, ...
                'UseParallel',false,'MaxIterations',500, ...
                'FunctionTolerance',1e-6);
        case 'patternsearch'
            disp('>>> Running PATTERNSEARCH optimization...');
            opts = optimoptions('patternsearch', ...
                'Display','iter','PollMethod','GSSPositiveBasis2N', ...
                'UseParallel',false,'MaxIterations',500, ...
                'FunctionTolerance',1e-6);
        case 'ga'
            disp('>>> Running GENETIC ALGORITHM optimization...');
             opts = optimoptions('ga', ...
                 'Display','iter','MaxGenerations',200,...
                 'PopulationSize',48, ...
                 'PlotFcn',{@gaplotbestf, @gaplotscores, @gaplotdistance, @gaplotbestindiv},'UseParallel',true);
            
        case 'particleswarm'
            disp('>>> Running PARTICLESWARM optimization...');
            opts = optimoptions('particleswarm', ...
                'Display','iter','SwarmSize',50, ...
                'MaxIterations',1000,'FunctionTolerance',1e-6, ...
                'HybridFcn',@fmincon,'UseParallel',true);
        otherwise
            error('Unknown solver name: %s', solverName);
    end
end

% Function that performs the main fitting task

function FitResults = runFitting(model, gData, estimatedParam, dose2s, opts, User_Inputs)
    % runFitting : Runs the optimization and estimates the parameters
     %
     % Inputs:
     %     model                       - The model object
     %     gData                       - Grouped data that is fed for
     %     fitting and is created from loadAndPrepareData function
     %     dose2s                      - Dose objects to be used for the
     %     fitting. It is output of function loadAndConfigModelForFit
     %     opts                        - Optimizer settings from
     %     getSolverOptions function
     %     User_Inputs.fit_pooled      - True in case of pooled fitting
     %     false for unpooled fitting
     %     User_Inputs.errormodel                  - Error function to be used for
     %     parameter estimations
     %     User_Inputs.Non_Exploratory_Params_Table     - Table of
     %     parameters that are fixed and not to be estimated
     %     User_Inputs.responseMap                 - Mapping between the
     %     model species and column headers of gData
     %
     % Outputs:
     %      FitResults              - Raw unprocessed estimated params
    
    fit_pooled = User_Inputs.fit_pooled;
    errormodel = User_Inputs.errormodel;
    Non_Exploratory_Params_Table = User_Inputs. Non_Exploratory_Params_Table;
    responseMap = User_Inputs.responseMap;
    %Parameter set up
    % Assign fixed parameters
    for i = 1:height(Non_Exploratory_Params_Table)
        set(sbioselect(model, 'Type', 'parameter', 'Name', Non_Exploratory_Params_Table.Name{i}), 'Value', Non_Exploratory_Params_Table.mu(i));
    end

    % set the seed first before parameter estimation
    rng(User_Inputs.seed);
    % Use sbiofit for parameter estimation
    FitResults = sbiofit(model, gData, responseMap, estimatedParam, dose2s, ...
                         'Pooled', fit_pooled, ...
                         'ErrorModel', errormodel,...
                         'Options', opts);
   
end

% Function that processes the fitted params and outputs them in an excel
% file

function [model_params_out, allEstimates] = processFitResults(FitResults, User_Inputs, Dose_mpk_list)
% processFitResults : Processes the fit results and saves the fit estimates
% in excel files
     %
     % Inputs:
     %     FitResults                       - The raw fit estimates
     %     obtained from runFitting function
     %     User_Inputs.fit_pooled           - True in case the fit is
     %     pooled and False in case of Unpooled fitting
     %     User_Inputs. Non_Exploratory_Params_Table    - The parameters
     %     that are not estimated/explored
     %     User_Inputs.Exploratory_Params_Table         - The parameters
     %     that the user has entered and that are estimated
     %     Dose_mpk_list                                - List of doses in
     %     milligram/kilogram that are present in the PK data and is output
     %     of loadAndConfigModelForFit function
     %
     % Outputs:
     %      model_params_out              - Parameters that are estimated
     %      and to be further used in simulation and visualizations
     %      allEstimates                  - Contains all estimated
     %      parameters
     % 
     % Files generated
     %      Pooled_Params_fit.xlsx        - Pooled parameter estimates when
     %      User_Inputs.fit_pooled is true
     %      Unpooled_Params_fit.xlsx      - Unpooled dose separated
     %      parameter values are stored in this file. Only for the case
     %      when User_Inputs.fit_pooled is false
     %      Unpooled_Mean_Params_fit.xlsx   - Unpooled mean parameter
     %      estimated. Includes mean (geometric) and variability across
     %      different dose groups. Only for the case
     %      when User_Inputs.fit_pooled is false



    fit_pooled = User_Inputs.fit_pooled;
    Non_Exploratory_Params_Table = User_Inputs. Non_Exploratory_Params_Table;
    Exploratory_Params_Table = User_Inputs.Exploratory_Params_Table;
    % Get BodyWeight from fixedTable
    BodyWeight = Non_Exploratory_Params_Table.mu(strcmp(Non_Exploratory_Params_Table.Name, 'W'));  

    % Initialize outputs
    allEstimates = []; % Raw (unscaled) individual estimates (for unpooled)
    
    % Get combined parameters (fixed and exploratory/fitted)
    % NOTE: The implementation of combine_params is assumed to be available
    model_params_out = combine_params(User_Inputs);
    
    % Extract names of parameters that were estimated/fitted
    fittedParamNames = Exploratory_Params_Table.Name;
    
    if isscalar(FitResults)
        % --- Pooled fit ---
        Pooled_Table = FitResults.ParameterEstimates(:, {'Name', 'Estimate'});
        
        % Extract raw (unscaled) mean estimates (U)
        U = Pooled_Table.Estimate';   
        
        % Update the mu column of the combined table with the fitted values
        for i = 1:numel(fittedParamNames)
             idx = strcmp(model_params_out.Name, fittedParamNames{i});
             model_params_out.mu(idx) = U(i);
        end
        
        % Scaling for Excel output (multiply by BodyWeight unless 'K_P_dec')
        Result_Table_for_Excel = Pooled_Table;
        paramNames = Result_Table_for_Excel.Name;
        paramsToScale = ~strcmp(paramNames, 'K_P_dec');
        Result_Table_for_Excel.Estimate(paramsToScale) = ...
            Result_Table_for_Excel.Estimate(paramsToScale) * BodyWeight;
        
        Result_Table = Result_Table_for_Excel; % Assign to output
        
        % --- ADDED: Add BodyWeight (W) to Pooled Table ---
        BodyWeight_Row = table({'W'}, BodyWeight, 'VariableNames', {'Name', 'Estimate'});
        Result_Table = [Result_Table; BodyWeight_Row];
        % ---------------------------------------------
        
    else
        % --- Unpooled fit (Individual Fits) ---
        nGroups = numel(FitResults);
        paramNames = FitResults(1).ParameterEstimates.Name;
        nParams = numel(paramNames);
        
        % Initialize matrix to store all raw (unscaled) estimates
        allEstimates = zeros(nGroups, nParams);
        
        for i = 1:nGroups
            allEstimates(i, :) = FitResults(i).ParameterEstimates.Estimate';
        end
        
        % Compute mean and standard deviation of raw (unscaled) log-transformed data
        mu_log = mean(log(allEstimates), 1);
        sigma_log = std(log(allEstimates), 0, 1);
        U = exp(mu_log);
        
        % Update the mu column of the combined table with the geometric mean of fitted values
        for i = 1:numel(fittedParamNames)
             idx = strcmp(model_params_out.Name, fittedParamNames{i});
             model_params_out.mu(idx) = U(i);
        end
        
        % Create Excel Table (Unpooled Estimates - Individual Fits)
        paramNames_col = paramNames(:); 
        allEstimates_scaled = allEstimates; % Start with raw estimates
        paramsToScale = ~strcmp(paramNames_col, 'K_P_dec');
        
        % Multiply by Bodyweight for all those params that are per kilogram
        allEstimates_scaled(:, paramsToScale) = ...
            allEstimates_scaled(:, paramsToScale) * BodyWeight;
            
        doseNames = strcat(string(Dose_mpk_list), " mpk");
        Unpooled_Table_Excel = array2table(allEstimates_scaled', ...
            'VariableNames', doseNames);
        Unpooled_Table_Excel.Name = paramNames_col;
        Unpooled_Table_Excel = movevars(Unpooled_Table_Excel, 'Name', 'Before', doseNames(1));
        
        Result_Table = Unpooled_Table_Excel; % Assign main table (individual fits) to output
        
        % Create a row where the Name is 'W' and all dose estimates are BodyWeight
        BodyWeight_Row_data = [BodyWeight, repmat(BodyWeight, 1, nGroups)]; 
        BodyWeight_Row_Name = {'W'};
        
        BodyWeight_Row = array2table(BodyWeight_Row_data, 'VariableNames', [{'Name'}, doseNames]);
        BodyWeight_Row.Name = BodyWeight_Row_Name;
        
        Result_Table = [Result_Table; BodyWeight_Row];
        
        % Create and save Mean Unpooled Params Table 
        U_scaled = U'; % column vector
        U_scaled(paramsToScale') = U_scaled(paramsToScale') * BodyWeight;        
        Mean_Unpooled_Table_Excel = table(paramNames_col, U_scaled, sigma_log', ...
            'VariableNames', {'Name', 'Mean', 'Sigma'});
            
        % Add BodyWeight (W) to Mean Unpooled Table
        BodyWeight_Mean_Row = table({'W'}, BodyWeight, 0, ...
            'VariableNames', {'Name', 'Mean', 'Sigma'});
        Mean_Unpooled_Table_Excel = [Mean_Unpooled_Table_Excel; BodyWeight_Mean_Row];
        % -------------------------------------------------------------------------------------
        
        filename3 = 'Unpooled_Mean_Params_fit.xlsx';
        if isfile(filename3)
            delete(filename3);
        end
        writetable(Mean_Unpooled_Table_Excel, filename3);
        disp(['Mean PK Unpooled estimates saved to ', filename3]);
        
    end
    
    % Save data to Excel files 
    if fit_pooled == true
        filename = 'Pooled_Params_fit.xlsx';
        if isfile(filename)
            delete(filename);
        end
        writetable(Result_Table, filename);
        disp(['Pooled PK estimates saved to ', filename])
    else
        filename2 = 'Unpooled_Params_fit.xlsx';
        if isfile(filename2)
            delete(filename2);
        end
        writetable(Result_Table, filename2);
        disp(['Unpooled dose separated PK estimates saved to ', filename2]);
    end
    
    %disp('Data saved to Excel files.');
end

% Function to check and raise error in case statistic in the
% Exploratory_Params_Table is not mu
function Check_mu_Present(User_Inputs)
% Check_mu_Present : Raises error in case user wants to explore sigma
% parameter
     %
     % Inputs:
     %      User_Inputs.Exploratory_Params_Table        - The parameters to
     %      be estimated
     %
     % Outputs:
     %      Raises error in case user tries to estimate the sigma
     %      (variability) of parameters
    Exploratory_Params_Table = User_Inputs. Exploratory_Params_Table;
    for i = 1:height(Exploratory_Params_Table)
        if (~strcmp(Exploratory_Params_Table.Statistic{i}, 'mu'))
            error('For PK the statistic should be only mu.');
        end
    end
end


% Get a combined table of the Exploratory and Non-Exploratory Params
function model_params = combine_params(User_Inputs)
    % combine_params :   combines parameters that are for estimation with
    % those that are not estimated to be used for the simulation/fitting
        % Inputs
        %   User_Inputs.Exploratory_Params_Table      - Parameters that are
        %   to be explored and/or estimated
        %   User_Inputs.Non_Exploratory_Params_Table      - Parameters that
        %   are NOT to be explored and/or estimated
        
        % Outputs
        %   model_params        - Combined parameters that are used for
        %   simulations

    Exploratory_Params_Table = User_Inputs.Exploratory_Params_Table;
    Non_Exploratory_Params_Table = User_Inputs.Non_Exploratory_Params_Table;
    % Initialize the combined table with non-exploratory (fixed) parameters
    model_params = Non_Exploratory_Params_Table;
    
    % Get exploratory parameters names
    exploratory_names = Exploratory_Params_Table.Name;
    exploratory_mu_vals = Exploratory_Params_Table.InitialValue;
    
    % Iterate over exploratory parameters
    for i = 1:numel(exploratory_names)
        current_param_name = exploratory_names{i};
        
        % Check if the parameter is already in the combined table (i.e., non-exploratory)
        idx_in_non_exp = strcmp(model_params.Name, current_param_name);
        
        if any(idx_in_non_exp)
            % Parameter is also in Non_Exploratory_Params_Table
            % Keep the sigma from Non_Exploratory_Params_Table
            % Overwrite the mu with the InitialValue from Exploratory_Params_Table
            model_params.mu(idx_in_non_exp) = exploratory_mu_vals(i);
        else
            % Parameter is only in Exploratory_Params_Table
            % The sigma for a fitted parameter (if not fixed) should be 0.
            % Since it's not present, sigma is set to 0 .
            model_params = [model_params;{current_param_name, exploratory_mu_vals(i), 0.0}]; %#ok<AGROW>
        end
    end
    % Ensure no duplicate names by taking the last unique entry
    [~, ix] = unique(model_params.Name, 'stable', 'last');
    model_params = model_params(ix,:);
end


% Function to run simulation
function runSimulationAndPlot(model_params, allEstimates, Dose_mpk_list, User_Inputs)
% runSimulationAndPlot : Generates the VPoP (in case of run_type = 'simulation') or runs the
% simulation for single patient (in case of run_type = 'manual'
% or run_type = 'auto_fitting') and plots the PK profiles across different doeses in
% different figures. Each figure contains 3 subplots (ADC, TAB and Payload in that order)
        %
        % Inputs
        %
        %   model_params       - Model Parameters that are used for the simulation and further visualization. In case of run_type = 'auto_fitting'
        %   estimated parameters are used for the simulation. Otherwise
        %   model_params are combined form of Exploratory and
        %   Non-Exploratory Params
        %   allEstimates       - All parameter estimates which is output
        %   of processFitResults function
        %   Dose_mpk_list      - List of doses in
        %   milligram/kilogram that are present in the PK data and is output
        %   of loadAndConfigModelForFit function
        %   User_Inputs.projfile    - Simbiology file
        %   User_Inputs.model_id    - Model index in Simbiology file
        %   User_Inputs.fit_pooled  - True for pooled fitting and false for
        %   unpooled fitting
        %   User_Inputs.plot_separate_unpooled  -   True in case of dose
        %   separated parameter based plots false for mean parameter plot
        %   User_Inputs.N_VP    - Number of patients in VPoP
        %   User_Inputs.Sampling_time       - Times (in hours) where
        %   simulation outputs are calculated
        %   User_Inputs.Exploratory_Params_Table       - Parameters that
        %   are to be estimated/explored
        %   User_Inputs.run_type            -   'auto_fitting',
        %   'manual_fitting' or 'simulation'
        %  
        % Outputs
        %   Generated PK profiles. In case of run_type =
        %   'auto_fitting'/'manual_fitting', single patient trajectories.
        %   In case of run_type = 'simulation', simulates for N_VP patients
        %   and plots curve for single patient with estimated or user given
        %   params and plots variability across different patients.

    projfile = User_Inputs.projfile;
    model_id = User_Inputs.model_id;
    fit_pooled = User_Inputs.fit_pooled;
    plot_separate_unpooled = User_Inputs.plot_separate_unpooled;
    N_VP = User_Inputs.N_VP;
    Sampling_time = User_Inputs.Sampling_time;
    Exploratory_Params_Table = User_Inputs.Exploratory_Params_Table;
    run_type = User_Inputs.run_type;

    % Get names of estimated and fixed parameters
    estimatedNames = Exploratory_Params_Table.Name;
    
    % Extract fixed parameter values
    BodyWeight = model_params.mu(strcmp(model_params.Name, 'W'));
    MW_ADC = model_params.mu(strcmp(model_params.Name, 'MW_ADC'));
        
    % All parameter names (fixed and estimated) and their mu/sigma
    allParamNames = model_params.Name;
    
    % Get indices for estimated parameters within the allParamNames array
    [~, estimated_indices] = ismember(estimatedNames, allParamNames);
    
    % --- Adjust N_VP and Sigma for Fitting Runs ---
    
    % Determine the number of virtual patients (N_VP_eff) and the sigma values to use
    if strcmp(run_type, 'auto_fitting') || strcmp(run_type, 'manual_fitting')
        N_VP_eff = 1; % Run only one patient
        sigma_for_sampling = zeros(size(model_params.sigma)); % Set all sigma to zero
    else
        N_VP_eff = N_VP; % Use the input N_VP
        sigma_for_sampling = model_params.sigma; % Use the standard sigma values
    end
    
    % Load SimBiology Project
    projStruct = sbioloadproject(projfile);
    model_var = "m" + string(model_id);   % creates "m1"
    model = projStruct.(model_var);
    csObj = getconfigset(model);
    csObj.StopTime = max(Sampling_time);
    csObj.TimeUnits = 'hour';
    csObj.SolverOptions.OutputTimes = Sampling_time;
    csObj.SolverType = 'ode15s';
    csObj.SolverOptions.RelativeTolerance = 1e-4;
    csObj.SolverOptions.AbsoluteTolerance = 1e-4;
    csObj.SolverOptions.MaxStep = 0.5;
    csObj.RuntimeOptions.StatesToLog = User_Inputs.StatesToLog;
    % DAR dosing
    DAR_dose = sbiodose('DAR_Dose', 'repeat');
    DAR_dose.Amount = 1;
    DAR_dose.AmountUnits = 'mole';
    DAR_dose.TimeUnits = 'hour';
    DAR_dose.StartTime = 0;
    DAR_dose.Interval = User_Inputs.Dose_Cycle_length;
    DAR_dose.RepeatCount = round(User_Inputs.Sampling_time(end)/User_Inputs.Dose_Cycle_length) + 1;
    DAR_dose.TargetName = User_Inputs.DAR_Dose_Target;

     % Copy this object to all the dose days in a cycle
    DAR_doses_multiple = [];
    for dose_day = User_Inputs.Doses_days_in_a_cycle
        dose_DAR = copyobj(DAR_dose);
        dose_DAR.StartTime = dose_day;
        DAR_doses_multiple = [DAR_doses_multiple, dose_DAR]; %#ok<AGROW>
    end
    DAR_dose = DAR_doses_multiple; clear dose_DAR DAR_doses_multiple;
    
    %Run Simulation Loop
    
    num_doses = numel(Dose_mpk_list);
    doseobj = cell(1, num_doses);
    
    for dIdx = 1:numel(Dose_mpk_list)
        rng(User_Inputs.seed);
        
        % ADC dosing
        Dose_mpk = Dose_mpk_list(dIdx);
        Dose1 = Dose_mpk * BodyWeight / 1000 / MW_ADC * 1e9;  % Convert to nanomole
        
        % Generate dose name (d1, d2, etc.)
        doseName = sprintf('d%d', dIdx);
        
        % Create and configure dose
        doseobj{dIdx} = sbiodose(doseName, 'repeat');
        doseobj{dIdx}.Amount = Dose1;
        doseobj{dIdx}.AmountUnits = 'nanomole';
        doseobj{dIdx}.TimeUnits = 'hour';
        doseobj{dIdx}.StartTime = 0;
        doseobj{dIdx}.Interval = User_Inputs.Dose_Cycle_length;
        doseobj{dIdx}.RepeatCount = round(User_Inputs.Sampling_time(end)/User_Inputs.Dose_Cycle_length) + 1;
        doseobj{dIdx}.TargetName = User_Inputs.Dose_Target;

        % Copy this object to all the dose days in a cycle
        real_doses = [];
        for dose_day = User_Inputs.Doses_days_in_a_cycle
            real_dose = copyobj(doseobj{dIdx});
            real_dose.StartTime = dose_day;
            real_doses = [real_doses,real_dose]; %#ok<AGROW>
        end
        doseobj{dIdx} = real_doses; clear real_doses real_dose;
        
        % Preallocate - uses the effective N_VP
        All_ADC = NaN(N_VP_eff, numel(Sampling_time));
        All_tAb = NaN(N_VP_eff, numel(Sampling_time));
        All_PL  = NaN(N_VP_eff, numel(Sampling_time));
        
        % Set up the mu_values for this specific dose (Current Mean)
        current_mu_values = model_params.mu; 
        
        % Logic for Pooled/Unpooled Plotting
        if (fit_pooled == true) 
            % No change needed to current_mu_values as model_params.mu contains U.
        elseif (~fit_pooled && plot_separate_unpooled)
            % For unpooled fit, use the dose-specific fit estimates (allEstimates(dIdx, :))
            if ~isempty(allEstimates) && size(allEstimates, 1) >= dIdx
                U_dose = allEstimates(dIdx, :); % raw (unscaled) dose-specific estimates
                current_mu_values(estimated_indices) = U_dose';
            end
        end
        % The sigma values are now selected based on run_type
        current_sigma_values = sigma_for_sampling;
        
        % --- VPoP Loop (uses N_VP_eff) ---
        for i = 1:N_VP_eff
            
            % Sample all parameters
            for j = 1:numel(allParamNames)
                p_name = allParamNames{j};
                mu = current_mu_values(j);
                sigma = 0;
                if N_VP_eff ~= 1
                    sigma = current_sigma_values(j); 
                end               
                % Sample from log-normal distribution
                val = mu * exp(randn * sigma); 
                
                % Assign to model
                set(sbioselect(model, 'Type', 'parameter', 'Name', p_name), 'Value', val);
            end
            
            % Simulate and store results
            sd = sbiosimulate(model, csObj, [], [doseobj{dIdx} DAR_dose]);
            conc_ADC = selectbyname(sd, User_Inputs.ADC_Model_Species).Data;
            conc_tAb = selectbyname(sd, User_Inputs.TAB_Model_Species).Data;
            conc_PL = selectbyname(sd, User_Inputs.PL_Model_Species).Data;
            All_ADC(i,:) = conc_ADC;
            All_tAb(i,:) = conc_tAb;
            All_PL(i,:)  = conc_PL;
        end
        
         % Compute stats
        epsilon = 1e-6;
        mean_ADC = max(mean(All_ADC,1),epsilon);
        min_ADC  = max(min(All_ADC,[],1),epsilon);
        max_ADC  = max(max(All_ADC,[],1),epsilon);
        mean_tAb = max(mean(All_tAb,1),epsilon);
        min_tAb  = max(min(All_tAb,[],1),epsilon);
        max_tAb  = max(max(All_tAb,[],1),epsilon);
        mean_PL  = max(mean( All_PL,1),epsilon);
        min_PL   = max(min( All_PL,[],1),epsilon);
        max_PL   = max(max( All_PL,[],1),epsilon);
       
        % Plotting
        PK_Summary = readtable(User_Inputs.PK_Summary_file);
        % Extract unique doses and analyte names
        % Load simulation results
        % Assumes variables: Sampling_time, mean_ADC, min_ADC, max_ADC, mean_tAb, min_tAb, max_tAb, mean_PL, min_PL, max_PL
        % Define Colors 
        gem_hex_codes = rgb2hex(orderedcolors("gem"));
        c_ADC = gem_hex_codes(1);   % deep blue
        c_tAb = gem_hex_codes(7);   % deep red
        c_PL  = gem_hex_codes(3);     % black

        % Select dose to overlay
        SelectedDose = Dose_mpk;
        % Unique doses and analyte names
        AnalyteNamesAll = unique(PK_Summary.Analyte);
        % Auto-map analyte names 
        ADC_name = AnalyteNamesAll(contains(upper(AnalyteNamesAll), "CONJUGATED"));
        tAb_name = AnalyteNamesAll(contains(upper(AnalyteNamesAll), "TOTAL"));
        PL_name  = AnalyteNamesAll(contains(upper(AnalyteNamesAll), "GSK5757810"));
        % Function to compute mean, min, max for given analyte
        getStats = @(tbl) arrayfun(@(t) struct( ...
            'time', t, ...
            'mean', mean(tbl.Conc(tbl.Time==t), 'omitnan'), ...
            'min',  min(tbl.Conc(tbl.Time==t)), ...
            'max',  max(tbl.Conc(tbl.Time==t))), ...
            unique(tbl.Time));
        
        % Plot as Subplots 
        figure;
        % ADC Plot
        subplot(1, 3, 1); hold on;
        Sampling_time_i = Sampling_time / 24;
        
        % --- Data Spread (fill area) ---
        if N_VP_eff ~= 1
            h1 = fill([Sampling_time_i, fliplr(Sampling_time_i)], ...
                      [min_ADC, fliplr(max_ADC)], ...
                      c_ADC, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        % --- Simulated Mean ADC ---
        h2 = semilogy(Sampling_time_i, mean_ADC, 'LineWidth', 2, ...
                      'Color', c_ADC, 'DisplayName', 'Simulated Mean ADC');
        
        % --- Observed Data (if available) ---
        PK_ADC = PK_Summary(PK_Summary.Dose == SelectedDose & ...
            contains(PK_Summary.Analyte, ADC_name, 'IgnoreCase', true), :);
        
        if ~isempty(PK_ADC)
            stats = getStats(PK_ADC);
            T = ([stats.time]'); M = [stats.mean]'; Mi = [stats.min]'; Ma = [stats.max]';
            Ti = T / 24;
            h3 = errorbar(Ti, M, M - Mi, Ma - M, 'o', ...
                'MarkerFaceColor', c_ADC, 'MarkerEdgeColor', c_ADC, ...
                'Color', c_ADC, 'DisplayName', 'Observed Data');
        else
            h3 = NaN;
        end
        
        % --- LLOQ line ---
        yLLOQ_ADC = 2.5;   % µg/mL
        h4 = yline(yLLOQ_ADC, '--k', sprintf('LLOQ = %.1f µg/mL', yLLOQ_ADC), ...
                   'LabelHorizontalAlignment', 'left', ...
                   'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
        % --- Labels, limits, formatting ---
        xlabel('Time (Days)', 'FontSize', 12);
        ylabel('Serum Conc. (µg/mL)', 'FontSize', 12);
        title(['ADC: ' num2str(Dose_mpk) ' mpk'], 'FontSize', 12);
        ylim([0.1 1000]);
        pbaspect([4 3 1]);
        set(gca, 'YScale', 'log');
        grid off; box on;
        
        % --- Legend (custom order) ---
        % if N_VP_eff ~= 1
        %     legend([h1 h2 h3 h4], ...
        %         {'Data Spread','Sim: Mean ADC', 'Observed Data', 'LLOQ'}, ...
        %         'Location', 'best', 'NumColumns', 2, 'FontSize', 8);
        % else
        %     legend([h2 h3 h4], ...
        %     {'Sim: Mean ADC', 'Observed Data', 'LLOQ'}, ...
        %     'Location', 'best', 'NumColumns', 2, 'FontSize', 8, 'Box','off');
        % end


        % --- tAb Plot ---
        subplot(1, 3, 2); hold on;
        Sampling_time_i = Sampling_time / 24;
        
        % Data Spread
        if N_VP_eff ~= 1
            h1_tAb = fill([Sampling_time_i, fliplr(Sampling_time_i)], ...
                          [min_tAb, fliplr(max_tAb)], ...
                          c_tAb, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % Simulated Mean tAb
        h2_tAb = semilogy(Sampling_time_i, mean_tAb, 'LineWidth', 2, ...
                          'Color', c_tAb, 'DisplayName', 'Simulated Mean tAb');
        
        % Observed Data
        PK_tAb = PK_Summary(PK_Summary.Dose == SelectedDose & ...
            contains(PK_Summary.Analyte, tAb_name, 'IgnoreCase', true), :);
        
        if ~isempty(PK_tAb)
            stats = getStats(PK_tAb);
            T = [stats.time]'; M = [stats.mean]'; Mi = [stats.min]'; Ma = [stats.max]';
            Ti = T / 24;
            h3_tAb = errorbar(Ti, M, M - Mi, Ma - M, 'o', ...
                'MarkerFaceColor', c_tAb, 'MarkerEdgeColor', c_tAb, ...
                'Color', c_tAb, 'DisplayName', 'Observed Data');
        else
            h3_tAb = NaN;
        end
        
        % LLOQ line
        yLLOQ_tAb = 2.5;   % µg/mL
        h4_tAb = yline(yLLOQ_tAb, '--k', sprintf('LLOQ = %.1f µg/mL', yLLOQ_tAb), ...
                       'LabelHorizontalAlignment', 'left', ...
                       'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
        
        % Formatting
        xlabel('Time (Days)', 'FontSize', 12);
        ylabel('Serum Conc. (µg/mL)', 'FontSize', 12);
        title(['tAb: ' num2str(Dose_mpk) ' mpk'], 'FontSize', 12);
        ylim([0.1 1000]);
        pbaspect([4 3 1]);
        set(gca, 'YScale', 'log');
        grid off; box on;
        
        % Custom Legend Order
        % if N_VP_eff ~= 1
        % legend([h1_tAb h2_tAb h3_tAb h4_tAb], ...
        %     {'Data Spread', 'Sim: Mean tAb', 'Observed Data', 'LLOQ'}, ...
        %     'Location', 'best', 'NumColumns', 2, 'FontSize', 8);
        % else
        %     legend([h2_tAb h3_tAb h4_tAb], ...
        %     {'Sim: Mean tAb', 'Observed Data', 'LLOQ'}, ...
        %     'Location', 'best', 'NumColumns', 2, 'FontSize', 8, 'Box', 'off');
        % end
        
        
        % --- PL Plot ---
        subplot(1, 3, 3); hold on;
        
        % Data Spread
        if N_VP_eff ~= 1
            h1_PL = fill([Sampling_time_i, fliplr(Sampling_time_i)], ...
                         [min_PL, fliplr(max_PL)], ...
                         c_PL, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        % Simulated Mean PL
        h2_PL = semilogy(Sampling_time_i, mean_PL, 'LineWidth', 2, ...
                         'Color', c_PL, 'DisplayName', 'Simulated Mean PL');
        
        % Observed Data
        PK_PL = PK_Summary(PK_Summary.Dose == SelectedDose & ...
            contains(PK_Summary.Analyte, PL_name, 'IgnoreCase', true), :);
        
        if ~isempty(PK_PL)
            stats = getStats(PK_PL);
            T = [stats.time]'; M = [stats.mean]'; Mi = [stats.min]'; Ma = [stats.max]';
            Ti = T / 24;
            h3_PL = errorbar(Ti, M, M - Mi, Ma - M, 'o', ...
                'MarkerFaceColor', c_PL, 'MarkerEdgeColor', c_PL, ...
                'Color', c_PL, 'DisplayName', 'Observed Data');
        else
            h3_PL = NaN;
        end
        
        % LLOQ line
        yLLOQ_PL = 5e-5;   % µg/mL
        h4_PL = yline(yLLOQ_PL, '--k', sprintf('LLOQ = %.5f µg/mL', yLLOQ_PL), ...
                      'LabelHorizontalAlignment', 'left', ...
                      'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
        
        % Formatting
        xlabel('Time (Days)', 'FontSize', 12);
        ylabel('Serum Conc. (µg/mL)', 'FontSize', 12);
        title(['PL: ' num2str(Dose_mpk) ' mpk'], 'FontSize', 12);
        ylim([1e-5 1e-1]);
        pbaspect([4 3 1]);
        set(gca, 'YScale', 'log');
        grid off; box on;
        
        % % Custom Legend Order
        % if N_VP_eff ~= 1
        %     legend([h1_PL h2_PL h3_PL h4_PL], ...
        %         {'Data Spread', 'Sim: Mean PL', 'Observed Data', 'LLOQ'}, ...
        %         'Location', 'best', 'NumColumns', 2, 'FontSize', 8);
        % else
        %     legend([h2_PL h3_PL h4_PL], ...
        %         {'Sim: Mean PL', 'Observed Data', 'LLOQ'}, ...
        %         'Location', 'best', 'NumColumns', 2, 'FontSize', 8, 'Box','off');
        % end
        fontsize(15,"points");
        sgtitle('Using Mean 10 mpk PK Parameters','FontSize',25,'FontWeight','bold')
        fontname('Consolas');

    end
end