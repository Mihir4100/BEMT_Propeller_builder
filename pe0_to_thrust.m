function [TotalThrust, TotalTorque, T] = pe0_to_thrust(pe0_filename, Vinf, RPM, polar_matfile, USE_TIP_LOSS)
% pe0_to_thrust  Read a .pe0 file and compute thrust & torque using the
% exact BEM formulation from your design code.
%
% Inputs:
%   pe0_filename   - string, path to .pe0 (or .dat) file
%   Vinf           - freestream speed (m/s)
%   RPM            - rotor speed (rpm)
%   polar_matfile  - .mat file containing polar_data and Re_vals
%   USE_TIP_LOSS   - boolean, true to include Prandtl tip/root losses
%
% Outputs:
%   TotalThrust    - total thrust (N)
%   TotalTorque    - total torque (N*m)
%   T              - table with node results (mirrors your script)
%
% Example:
%   [Tthrust, Ttorq, T] = pe0_to_thrust("4x4E-3.dat", 20, 14000, "S1223_polars.mat", false);

% --- Constants & defaults
if nargin < 5
    USE_TIP_LOSS = false;
end

% fluid/material constants (same as your script)
rho = 1.225;          % kg/m^3
mu  = 1.81e-5;        % Pa*s

% convert RPM to rad/s (same formula)
Omega = (RPM*(pi/30));

% Unit conversions
inch2m = 0.0254;

%% --- PARSE .pe0 FILE ---
fid = fopen(pe0_filename,'r');
if fid == -1
    error("Cannot open file %s", pe0_filename);
end

station_vals = [];
R_radius_in = NaN;
hubtra_in = NaN;
B = NaN;
density_sg = NaN;   % specific gravity if present

tline = fgetl(fid);
while ischar(tline)
    % try extract RADIUS, HUBTRA, BLADES, DENSITY if present
    if isempty(R_radius_in)
        tokens = regexp(tline, 'RADIUS:\s*([0-9.+-Ee]+)', 'tokens');
        if ~isempty(tokens); R_radius_in = str2double(tokens{1}{1}); end
    end
    if isempty(hubtra_in)
        tokens = regexp(tline, 'HUBTRA:\s*([0-9.+-Ee]+)', 'tokens');
        if ~isempty(tokens); hubtra_in = str2double(tokens{1}{1}); end
    end
    if isnan(B)
        tokens = regexp(tline, 'BLADES:\s*([0-9]+)', 'tokens');
        if ~isempty(tokens); B = str2double(tokens{1}{1}); end
    end
    if isnan(density_sg)
        tokens = regexp(tline, 'DENSITY \(SPECIFIC GRAVITY, INPUT FILE\)\s*=\s*([0-9.+-Ee]+)', 'tokens');
        if ~isempty(tokens); density_sg = str2double(tokens{1}{1}); end
    end

    % Try read a station row: 13 floats
    nums = sscanf(tline, '%f %f %f %f %f %f %f %f %f %f %f %f %f')';
    if numel(nums) == 13
        station_vals(end+1, :) = nums; %#ok<AGROW>
    end

    tline = fgetl(fid);
end
fclose(fid);

if isempty(station_vals)
    error("No station rows found in %s", pe0_filename);
end
% station_vals columns (as in your file):
% 1: station_in, 2: chord_in, 3: pitch_quoted_deg, 4: pitch_le-te_deg,
% 5: pitch_prather_deg, 6: sweep_in, 7: thickness_ratio, 8: twist_deg,
% 9: max_thick_in, 10: cross_section_area_in2, 11: zhigh_in, 12: cgy_in, 13: cgz_in

% extract arrays
station_in = station_vals(:,1);
chord_in   = station_vals(:,2);
% twist_deg etc if needed:
twist_deg  = station_vals(:,8);
F_arr_file = ones(size(station_in));  % we'll compute or overwrite

% convert to SI
r_m = station_in * inch2m;
chord_m = chord_in * inch2m;
% If hub radius from file exists, convert:
if ~isnan(hubtra_in)
    hubRadius = hubtra_in * inch2m;
else
    % fallback: use minimum station radius as hub (small)
    hubRadius = min(r_m);
end
% If R stored, convert
if ~isnan(R_radius_in)
    R = R_radius_in * inch2m;
else
    R = max(r_m);
end

% if B was missing, default to 3 (but .pe0 generally includes)
if isnan(B); B = 3; end

%% --- set r vector to station nodes from file (use those nodes directly) ---
r = r_m;   % node locations in meters (same length as station rows)

% initialize variables consistent with your script
a = 0.25 * ones(size(r));
ap = 0.01 * ones(size(r));
phi = atan2(Vinf, (Omega .* r));    % [rad]
theta = zeros(size(r));

Re_arr = zeros(size(r));
F_arr = ones(size(r));

maxIter = 1000;
tol = 1e-6;
relax = 0.12;

% load polars (same as your script)
if exist(polar_matfile,'file')
    load(polar_matfile, 'polar_data', 'Re_vals'); %#ok<NODEF>
else
    error("Polar matfile %s not found", polar_matfile);
end

% precompute min/max Alpha bounds as in your script
Alpha_min = cellfun(@(T) min(T.Alpha), polar_data);
Alpha_max = cellfun(@(T) max(T.Alpha), polar_data);

%% --- MAIN LOOP (exact math) ---
for i = 1:length(r)
    if i > 1
        a(i) = a(i-1);
        ap(i) = ap(i-1);
        phi(i) = phi(i-1);
    end

    c = chord_m(i);   % local chord from file

    for kIter = 1:maxIter
        phi_old = phi(i);

        Va = Vinf * (1 + a(i));
        Vt = Omega * r(i) * (1 - ap(i));
        W = sqrt(Va^2 + Vt^2);
        Re = rho * W * c / mu;

        % --- Get best design Alpha for this element (deg) ---
        Alpha_best = getBestAlpha(Re, Re_vals, polar_data);   % user-supplied helper

        % set twist (theta) so target AoA is Alpha_best
        theta(i) = deg2rad(Alpha_best) + phi(i);   % theta in radians

        % compute actual local AoA (deg) for interpolation
        Alpha_local = rad2deg(theta(i) - phi(i));

        % clamp Alpha_local to polar bounds (prevents crazy extrapolation)
        iLow = find(Re_vals <= Re, 1, 'last'); if isempty(iLow), iLow = 1; end
        iHigh = find(Re_vals >= Re, 1, 'first'); if isempty(iHigh), iHigh = length(Re_vals); end
        global_Alpha_min = min(Alpha_min(iLow), Alpha_min(iHigh));
        global_Alpha_max = max(Alpha_max(iLow), Alpha_max(iHigh));
        Alpha_local = max(global_Alpha_min, min(global_Alpha_max, Alpha_local));

        % --- Interpolate CL and CD at (Alpha_local, Re) ---
        [CL_i, CD_i] = interpolatePolar(Alpha_local, Re, Re_vals, polar_data);

        % --- Blade-element forces ---
        Cn = CL_i * cos(phi(i)) - CD_i * sin(phi(i));
        Ct = CL_i * sin(phi(i)) + CD_i * cos(phi(i));
        sigma = B * c / (2 * pi * r(i));

        % --- Prandtl loss ---
        if USE_TIP_LOSS
            sphi = sin(phi(i));
            f_tip = (B*(R - r(i))) / (2 * R * sphi);
            f_root = (B*(r(i) - hubRadius)) / (2 * r(i) * sphi);
            val_tip = max(0, min(1, exp(-f_tip)));
            val_root = max(0, min(1, exp(-f_root)));
            F_tip = (2/pi) * acos(val_tip);
            F_root = (2/pi) * acos(val_root);
            F = max(F_tip * F_root, 0);
        else
            F = 1;
        end
        F_arr(i) = F;

        CT = sigma * Cn / (sin(phi(i))^2);

        % Standard Glauert correction threshold
        a_Glauert_thresh = 0.32;

        % Compute intermediate values
        K = 4 * F * sin(phi(i))^2 / (sigma * Cn);
        a_ideal = 1 / (K - 1);

        if a_ideal <= a_Glauert_thresh
            anew = a_ideal;
        else
            % Glauert's empirical high-thrust correction (see Glauert 1935)
            CT1 = 1.816;  % empirical max value before non-physical region
            CT2 = 2 * F * (1 - a_Glauert_thresh)^2;

            if CT < CT1
                anew = (CT2 / CT) * (1 - sqrt(1 - CT / CT2));
            else
                anew = 0.95;  % Limit to physical max
            end
        end

        % Clamp
        anew = max(0.0, min(0.95, real(anew)));

        % --- Tangential induction factor ap ---
        denom_ap = 4 * F * sin(phi(i)) * cos(phi(i));
        if abs(denom_ap) < 1e-8
            apnew = ap(i);  % Avoid division by zero
        else
            apnew = sigma * Ct / (denom_ap + sigma * Ct);
        end

        % Clamp ap to stay within physical range
        apnew = max(0.0, min(0.5, real(apnew)));

        % Update phi, a, ap with relaxation
        Va = Vinf * (1 + anew);
        Vt = Omega * r(i) * (1 - apnew);
        phi_new = atan2(Va, Vt);

        phi(i) = (1 - relax) * phi_old + relax * phi_new;
        a(i)   = (1 - relax) * a(i)   + relax * anew;
        ap(i)  = (1 - relax) * ap(i)  + relax * apnew;

        if abs(phi(i) - phi_old) < tol
            break;
        end
    end % kIter

    Re_arr(i) = Re;
end % i

theta = smoothdata(theta, 'movmean', 5);

% --- Final output table and thrust calc (unchanged formulas, SI units) ---
T = table(r, r - hubRadius, rad2deg(theta), a', ap', rad2deg(phi'), chord_m, Re_arr', F_arr', ...
    'VariableNames', {'r_m','r_for_Qblade','theta_deg','a','aPrime','phi_deg','chord_m','Re','F'});

% --- FIX 3: CORRECT THRUST CALCULATION --- (same exact form as your script)
dT_dr = 4 * pi * rho * Vinf^2 .* T.a .* (1 + T.a) .* T.r_m .* T.F;

% store node values
T.dT_dr = dT_dr;

% --- thrust distribution per section (trapezoid) ---
r_nodes = T.r_m;
dr = diff(r_nodes);
dT_section = 0.5*(dT_dr(1:end-1) + dT_dr(2:end)) .* dr;

% create a separate table column (pad size)
T.dT_section = [dT_section; NaN];

% === total thrust (momentum) ===
TotalThrust = trapz(r_nodes, dT_dr);    % integrates same formulation
fprintf('\nTotal Calculated Thrust: %.2f N\n', TotalThrust);

% === also test sum of sections ===
TotalThrust_sections = sum(dT_section);
fprintf('Total via summing sections: %.2f N\n', TotalThrust_sections);

% --- TORQUE (momentum) (same as your script) ---
dQ_dr = 4 * pi * rho * Vinf * Omega .* (T.r_m.^3) .* T.aPrime .* (1 + T.a) .* T.F;

T.dQ_dr = dQ_dr;

dQ_section = 0.5 .* (dQ_dr(1:end-1) + dQ_dr(2:end)) .* dr;
T.dQ_section = [dQ_section; NaN];

TotalTorque = trapz(r_nodes, dQ_dr);
fprintf('Total Calculated Torque: %.3f N·m\n', TotalTorque);

TotalTorque_sections = sum(dQ_section);
fprintf('Total via summing sections: %.3f N·m\n', TotalTorque_sections);

% final outputs
if nargout < 1
    clear TotalThrust;
end
if nargout < 2
    clear TotalTorque;
end

end
