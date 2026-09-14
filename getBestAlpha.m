function Alpha_opt = getBestAlpha(Re, Re_vals, polar_data)
    % Clamp Re to valid range
    Re = max(min(Re, max(Re_vals)), min(Re_vals));

    % Safe index resolution
    iLow = find(Re_vals <= Re, 1, 'last');
    if isempty(iLow), iLow = 1; end
    iHigh = find(Re_vals >= Re, 1, 'first');
    if isempty(iHigh), iHigh = length(Re_vals); end

    % Make sure indices are within polar_data range
    iLow = max(1, min(iLow, length(polar_data)));
    iHigh = max(1, min(iHigh, length(polar_data)));

    % Proceed as before
    if iLow == iHigh
        data = polar_data{iLow};
        [~, idx] = max(data.Cl);
        Alpha_opt = data.Alpha(idx);
    else
        dataL = polar_data{iLow};
        dataH = polar_data{iHigh};
        [~, idxL] = max(dataL.Cl);
        [~, idxH] = max(dataH.Cl);
        AlphaL = dataL.Alpha(idxL);
        AlphaH = dataH.Alpha(idxH);
        w = (Re - Re_vals(iLow)) / (Re_vals(iHigh) - Re_vals(iLow));
        Alpha_opt = (1 - w) * AlphaL + w * AlphaH;
    end
end
