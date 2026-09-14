function [CL, CD] = interpolatePolar(alpha_deg, Re, Re_vals, polar_data)
    % CLamp Re
    Re = max(min(Re, max(Re_vals)), min(Re_vals));
    iLow = find(Re_vals <= Re, 1, 'last'); if isempty(iLow), iLow = 1; end
    iHigh = find(Re_vals >= Re, 1, 'first'); if isempty(iHigh), iHigh = length(Re_vals); end

    % Get tables
    tblL = polar_data{iLow};
    tblH = polar_data{iHigh};

    % CLean and prepare each table (helper)
    tblL = CLeanPolarTable(tblL, iLow);
    tblH = CLeanPolarTable(tblH, iHigh);

    % Interpolate alpha within each Re table
    % interp1 requires finite numeric vectors and no NaNs in sample points
    CLL = interp1(tblL.alpha, tblL.CL, alpha_deg, 'linear', 'extrap');
    CLH = interp1(tblH.alpha, tblH.CL, alpha_deg, 'linear', 'extrap');
    CDL = interp1(tblL.alpha, tblL.CD, alpha_deg, 'linear', 'extrap');
    CDH = interp1(tblH.alpha, tblH.CD, alpha_deg, 'linear', 'extrap');

    % linear interp in Re
    if iLow == iHigh
        w = 0;
    else
        w = (Re - Re_vals(iLow)) / (Re_vals(iHigh) - Re_vals(iLow));
    end
    CL = (1 - w) * CLL + w * CLH;
    CD = (1 - w) * CDL + w * CDH;
end

function T = CLeanPolarTable(T, idx)
    % Ensure expected column names (case-insensitive)
    wantNames = {'alpha','CL','CD'};
    vars = T.Properties.VariableNames;

    for iWant = 1:numel(wantNames)
        want = wantNames{iWant};
        % find a column whose name matches want (case-insensitive)
        j = find(strcmpi(vars, want), 1);
        if ~isempty(j)
            T.Properties.VariableNames{j} = want;  % normalize to canonical name
        else
            % try to auto-detect common variants
            % e.g., 'alpha', 'alpha', 'a' etc. (you can expand heuristics if needed)
            % if missing, we'll leave it and catch it later
        end
    end

    % After renaming attempt, check we have required columns
    if ~all(ismember(wantNames, T.Properties.VariableNames))
        % show helpful diagnostic before failing
        error(['polar table %d missing required columns. Found: %s. ', ...
               'Expected alpha, CL, CD (case-insensitive).'], ...
               idx, strjoin(T.Properties.VariableNames, ','));
    end

    % Convert to numeric if necessary (handles tables with text cells)
    T.alpha = double(T.alpha);
    T.CL    = double(T.CL);
    T.CD    = double(T.CD);

    % Remove rows with non-finite entries
    finiteMask = isfinite(T.alpha) & isfinite(T.CL) & isfinite(T.CD);
    if ~all(finiteMask)
        warning('polar table %d: removing %d non-finite rows from polar.', idx, sum(~finiteMask));
        T = T(finiteMask,:);
    end

    % Sort by alpha (interp1 expects X to be monotonic)
    [~, ia] = sort(T.alpha);
    T = T(ia,:);

    % Remove duplicate alpha values (average CL/CD for duplicates)
    [uniqueA, ~, ic] = unique(T.alpha, 'stable');
    if numel(uniqueA) < height(T)
        CL_agg = accumarray(ic, T.CL)./accumarray(ic, 1);
        CD_agg = accumarray(ic, T.CD)./accumarray(ic, 1);
        T = table(uniqueA, CL_agg, CD_agg, 'VariableNames', {'alpha','CL','CD'});
    end

    % final sanity check
    if isempty(T) || ~all(isfinite(T.alpha))
        error('polar table %d is empty or invalid after CLeaning — check CSV import.', idx);
    end
end
