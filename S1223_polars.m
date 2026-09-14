    %% --- REBUILD POLARS ROBUSTLY WITH AUTO NUMERIC CONVERSION ---
clear; clc;

% Reynolds numbers corresponding to your CSV files
Re_vals = [30000 35000 40000 45000 50000 55000 60000 65000 70000 75000 ...
           80000 85000 90000 95000 100000 105000 110000 115000 120000 125000 ...
           130000 135000 140000 145000 150000];

% List of all CSVs (ensure these files exist in the current folder)
file_names = {
   "S1223_T3_Re0.030_M0.06_N9.0.csv"
   "S1223_T3_Re0.035_M0.06_N9.0.csv"
   "S1223_T3_Re0.040_M0.06_N9.0.csv"
   "S1223_T3_Re0.045_M0.06_N9.0.csv"
   "S1223_T3_Re0.050_M0.06_N9.0.csv"
   "S1223_T3_Re0.055_M0.06_N9.0.csv"
   "S1223_T3_Re0.060_M0.06_N9.0.csv"
   "S1223_T3_Re0.065_M0.06_N9.0.csv"
   "S1223_T3_Re0.070_M0.06_N9.0.csv"
   "S1223_T3_Re0.075_M0.06_N9.0.csv"
   "S1223_T3_Re0.080_M0.06_N9.0.csv"
   "S1223_T3_Re0.085_M0.06_N9.0.csv"
   "S1223_T3_Re0.090_M0.06_N9.0.csv"
   "S1223_T3_Re0.095_M0.06_N9.0.csv"
   "S1223_T3_Re0.100_M0.06_N9.0.csv"
   "S1223_T3_Re0.105_M0.06_N9.0.csv"
   "S1223_T3_Re0.110_M0.06_N9.0.csv"
   "S1223_T3_Re0.115_M0.06_N9.0.csv"
   "S1223_T3_Re0.120_M0.06_N9.0.csv"
   "S1223_T3_Re0.125_M0.06_N9.0.csv"
   "S1223_T3_Re0.130_M0.06_N9.0.csv"
   "S1223_T3_Re0.135_M0.06_N9.0.csv"
   "S1223_T3_Re0.140_M0.06_N9.0.csv"
   "S1223_T3_Re0.145_M0.06_N9.0.csv"
   "S1223_T3_Re0.150_M0.06_N9.0.csv"
   };

polar_data = cell(1, numel(file_names));

%% --- Loop through each file ---
for k = 1:numel(file_names)
    fname = file_names{k};

    % read entire CSV with auto detection
    Traw = readtable(fname, 'PreserveVariableNames', true);

    
    % safer to force case-insensitive match:
    vn = lower(Traw.Properties.VariableNames);

    idxA = find(strcmp(vn,'alpha'),1);
    idxCl = find(strcmp(vn,'cl'),1);
    idxCd = find(strcmp(vn,'cd'),1);
    
    if isempty(idxA) || isempty(idxCl) || isempty(idxCd)
        error('Missing required col in %s',fname);
    end
    
    Alpha = Traw.(Traw.Properties.VariableNames{idxA});
    Cl    = Traw.(Traw.Properties.VariableNames{idxCl});
    Cd    = Traw.(Traw.Properties.VariableNames{idxCd});



    T = table(Alpha,Cl,Cd);

    % clean + sort
    T = sortrows(T(isfinite(T.Alpha) & isfinite(T.Cl) & isfinite(T.Cd),:),"Alpha");

    polar_data{k} = T;
    fprintf("Loaded %s (%d rows)\n", fname, height(T));
end

save('S1223_polars.mat','polar_data','Re_vals');

disp('Rebuilt S1223_polars.mat successfully (all numeric data)');
