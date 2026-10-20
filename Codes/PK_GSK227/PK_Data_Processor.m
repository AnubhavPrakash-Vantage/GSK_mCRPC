% This code is to get the raw .csv file in SDTM format and convert it to
% simpler .mat/.xlsx table to be used for PK calibrations. The raw file is stored
% in the google drive.

clear; 

% Enter the raw data file name in SDTM format to be cleaned up
raw_pk_data = "G:\Shared drives\GSK_Oncology\B7_H3_Pilot_Project\06 Slides\Shared_by_client_(No Downloading)\GSK FTIH 14thDec2025\PK_3analytes_GSK227_223054_DCOOct25.csv";

% Call the funtion to clean the data
PK_Extractor(raw_pk_data);





















function PK_Extractor(raw_pk_data)
     % PK_Extractor: Extracts data from the SDTM .csv file and generates two
     % files to be used for calibrating data and for simulations
     %
     % Inputs:
     %   raw_pk_data           - location of the SDTM csv raw file
     %
     % Generates and Saves:
     %   PK_Data_for_Calibration.xlsx          - Saves the mean profiles of
     %   all three analytes: ADC, TAB and PL along with their min and max vals
     %   PK_Summary.xlsx        - Creates a summary file that is used for
     %   plotting in next step


    % Import and prepare PK data table
    warning off;
    % Read the PK dataset
    PK_Datai = readtable(raw_pk_data);
    % Subset relevant columns
    PK_Data = PK_Datai(:, {'USUBJIDN','AFRLT','EVID','DVID','DVIDN','DV','DOSEP','CANCERC'});
    PK_Data = PK_Data(PK_Data.USUBJIDN == 223054000121,:);
    PK_Data = PK_Data(strcmp(PK_Data.CANCERC,'PROSTATE_CANCER'),:);
    % Find unique dosing times
    DosingTime_i = PK_Data.EVID ~= 0;
    DosingTime = unique(PK_Data.AFRLT(DosingTime_i));
    % Define your unique doses and analyte identifiers
    Dose_u        = unique(PK_Data.DOSEP(~isnan(PK_Data.DOSEP)));
    Analyte_u     = unique(PK_Data.DVID);
    Analyte_u(1, :)=[];
    AnalyteNames  = cellstr(Analyte_u);  % readable names
    % Create and save summary table
    % Initialize an empty table to store processed PK data
    PK_Summary = table();
    % Loop over doses
    for d = 1:numel(Dose_u)
        DoseVal = Dose_u(d);
        PK_Dose = PK_Data(PK_Data.DOSEP == DoseVal, :);
        if isempty(PK_Dose)
            continue;
        end
        
        % Loop over analytes
        for a = 1:numel(Analyte_u)
            analyteID   = Analyte_u{a};
            analyteName = AnalyteNames{a};
            
            PK_Analyte = PK_Dose(strcmp(PK_Dose.DVID, analyteID), :);
            if isempty(PK_Analyte)
                continue;
            end
            
            % Create temporary table for this subset
            tmpTbl = table();
            tmpTbl.Time       = PK_Analyte.AFRLT;
            tmpTbl.Conc       = PK_Analyte.DV;
            tmpTbl.Dose       = PK_Analyte.DOSEP;
            tmpTbl.Dosing_Time= PK_Analyte.EVID;
            tmpTbl.SUBID      = PK_Analyte.USUBJIDN;
            tmpTbl.CancerType = PK_Analyte.CANCERC;
            tmpTbl.Analyte    = repmat({analyteName}, height(PK_Analyte), 1);
            
            % Append to summary table
            PK_Summary = [PK_Summary; tmpTbl]; %#ok<AGROW>
        end
    end
    % Find analytes containing 'GSK5757810' bcz this is in ng/ml...this to be
    % be converted to ug/ml
    idx_convert = contains(PK_Summary.Analyte, 'GSK5757810');
    % Convert their concentrations from ng/mL to ug/mL (divide by 1000)
    PK_Summary.Conc(idx_convert) = PK_Summary.Conc(idx_convert) / 1000;
    % Save processed data
    save('PK_Summary.mat', 'PK_Summary');
    writetable(PK_Summary, 'PK_Summary.xlsx');
    disp('PK_Summary.mat has been saved successfully.');
    %Load saved table
    % Load saved table
    load('PK_Summary.mat'); %#ok<LOAD>
    % Define unique doses and analytes again (from saved data)
    Dose_u       = unique(PK_Summary.Dose);
    AnalyteNames = unique(PK_Summary.Analyte);
    
    % Create a summary table for all doses and analytes
    % Preallocate an empty table
    PK_Table = table();
    for d = 1:numel(Dose_u)
        DoseVal = Dose_u(d);
        PK_Dose = PK_Summary(PK_Summary.Dose == DoseVal, :);
        if isempty(PK_Dose)
            continue;
        end
        % Get all unique time points for this dose
        Time_u = unique(PK_Dose.Time);
        for t = 1:numel(Time_u)
            TimeVal = Time_u(t);
            row.DoseID = d;
            if ismember(TimeVal, DosingTime)
                row.Dose_mpk = DoseVal;
            else
                row.Dose_mpk = NaN;
            end
            row.Time = TimeVal;
            % Loop through analytes
            for a = 1:numel(AnalyteNames)
                analyteName = AnalyteNames{a};
                PK_Analyte = PK_Dose(strcmp(PK_Dose.Analyte, analyteName) & PK_Dose.Time == TimeVal, :);
                if isempty(PK_Analyte)
                    row.(['Mean_' analyteName]) = NaN;
                    row.(['Min_'  analyteName]) = NaN;
                    row.(['Max_'  analyteName]) = NaN;
                else
                    vals = PK_Analyte.Conc;
                    row.(['Mean_' analyteName]) = mean(vals);
                    row.(['Min_'  analyteName]) = min(vals);
                    row.(['Max_'  analyteName]) = max(vals);
                end
            end
            % Convert struct to table and append
            PK_Table = [PK_Table; struct2table(row)]; %#ok<AGROW>
        end
    end
    save('PK_Data_for_Calibration.mat', 'PK_Table');
    writetable(PK_Table, 'PK_Data_for_Calibration.xlsx')
    disp('PK summary table created.');
    writetable(PK_Table, 'my_PK_table.xlsx');
end
