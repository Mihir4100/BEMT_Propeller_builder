clear; clc;

% --- Setup Parameters ---
% axial velocity at blade inlet
Vinf = 20; 
% rotational velocity of propeller (rad/s)
Omega = 1200; 
% Number of blades
B = 6; 
% Radius of propeller
R = 0.075;
% density of air
rho = 1.225;
% kinematic coefficient of viscosity
mu = 1.81e-5; 
USE_TIP_LOSS = false;
hubRadius = 0.015;
% radius of section being analysed
r = linspace(hubRadius, R, 20); 
% initialize
a = 0.25*ones(size(r));
ap = 0.01*ones(size(r));
c = zeros(size(r));
phi = atan2(Vinf, (Omega.*r));    % [rad]
theta = zeros(size(r));

Re_arr = zeros(size(r));
F_arr = ones(size(r));

maxIter = 1000;
tol = 1e-6;
relax = 0.12;

% --- LOAD POLARS ONCE ---
load("S1223_polars.mat");   % provides: polar_data (cell of tables), Re_vals

% (Optional) precompute min/max Alpha in each polar to clamp later:
Alpha_min = cellfun(@(T) min(T.Alpha), polar_data);
Alpha_max = cellfun(@(T) max(T.Alpha), polar_data);

%% --- MAIN LOOP ---
for i = 1:length(r)
    if i > 1
        a(i) = a(i-1);
        ap(i) = ap(i-1);
        phi(i) = phi(i-1);
    end

    c(i) = 0.0275 - 0.01* r(i) / R;

    for kIter = 1:maxIter
        phi_old = phi(i);

        Va = Vinf * (1 + a(i));
        Vt = Omega * r(i) * (1 - ap(i));
        W = sqrt(Va^2 + Vt^2);
        Re = rho * W * c(i) / mu;

        % --- Get best design Alpha for this element (deg) ---
        Alpha_best = getBestAlpha(Re, Re_vals, polar_data);   % separate file getBestAlpha.m

        % set twist (theta) so target AoA is Alpha_best
        theta(i) = deg2rad(Alpha_best) + phi(i);   % theta in radians

        % compute actual local AoA (deg) for interpolation
        Alpha_local = rad2deg(theta(i) - phi(i));

        % clamp Alpha_local to polar bounds (prevents crazy extrapolation)
        % find nearest Re index for bounds (we only need min/max Alpha range)
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
        sigma = B * c(i) / (2 * pi * r(i));

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
        a_Glauert_thresh = 0.33;

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
% From angular momentum balance

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
    end %

    Re_arr(i) = Re;
end % i

theta = smoothdata(theta, 'movmean', 5);

% --- Final output table and thrust calc (unchanged) ---
T = table(r', rad2deg(theta)', a', ap', rad2deg(phi'), c', Re_arr', F_arr', ...
    'VariableNames', {'r_m','twist angle','a','aPrime','inflow angle','chord','Re','F'});

% --- FIX 3: CORRECT THRUST CALCULATION ---
% Must include the Prandtl loss factor F in the momentum equation
dT_dr = 4 * pi * rho * Vinf^2 .* T.a .* (1 + T.a) .* T.r_m .* T.F;
% --- END FIX 3 ---

% --- thrust distribution per section (trapezoid) ---
r = T.r_m;
dr = diff(r);
dT_section = 0.5*(dT_dr(1:end-1) + dT_dr(2:end)) .* dr;

% create a separate table for section results OR append as last N-1 rows
T.dT_section = [dT_section; NaN];   % pad last row with NaN so table sizes match



% === total thrust (momentum) ===
TotalThrust = trapz(r, dT_dr);    % integrates exactly same formulation
fprintf('\nTotal Calculated Thrust: %.2f N\n', TotalThrust);

% === also test sum of sections ===
TotalThrust_sections = sum(dT_section);
fprintf('Total via summing sections: %.2f N\n', TotalThrust_sections);

% --- TORQUE (momentum) ---
% Using: dQ/dr = 4*pi*rho*Vinf*Omega * r^3 * a' * (1 - a) * F
lambda_r = Omega .* T.r_m ./ Vinf;    % (only if you want to keep the explicit form)
dQ_dr = 4 * pi * rho * Vinf * Omega .* (T.r_m.^3) .* T.aPrime .* (1 + T.a) .* T.F;

% --- torque per section (trapezoid over r) ---
r  = T.r_m;
dr = diff(r);
dQ_section = 0.5 * (dQ_dr(1:end-1) + dQ_dr(2:end)) .* dr;

% append (pad last row to keep table sizes aligned)
T.dQ_section = [dQ_section; NaN];

% === total torque (momentum) ===
TotalTorque = trapz(r, dQ_dr);                % N·m
fprintf('Total Calculated Torque: %.3f N·m\n',TotalTorque);

% === also test sum of sections ===
TotalTorque_sections = sum(dQ_section);       % N·m
fprintf('Total via summing sections: %.3f N·m\n',TotalTorque_sections);

disp(T);
