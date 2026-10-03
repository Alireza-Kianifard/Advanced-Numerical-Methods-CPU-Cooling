clear; clc; close all;

%% 1. Parameters Definition
% Geometry
Lx          = 0.02;             % Width (m)
Ly          = 0.01;             % Height (m)

% Material Properties (Aluminum)
k           = 237;              % Thermal conductivity (W/m.K)
rho         = 2700;             % Density (kg/m^3)
cp          = 900;              % Specific heat (J/kg.K)
alpha       = k / (rho * cp);   % Thermal diffusivity (m^2/s)

% Boundary & Initial Conditions
q_flux      = 200000;           % Bottom constant heat flux (W/m^2)
h_side      = 64;               % Side natural convection (W/m^2.K)
T_inf       = 25;               % Ambient temperature (C)
T_init      = 25;               % Initial temperature (C)

% Cooling Parameters (Smart Fan)
h_base      = 1000;             % Base convection coefficient (fan idle)
T_target    = 55;               % Temperature threshold to trigger cooling (C)
Kp          = 100;              % Proportional gain for controller
h_max_limit = 4000;             % Maximum capacity of the cooling fan

% Mesh Grid (for Explicit and ADI)
Nx          = 81;
Ny          = 41;
dx          = Lx / (Nx - 1);
dy          = Ly / (Ny - 1);
x           = linspace(0, Lx, Nx);
y           = linspace(0, Ly, Ny);
[X, Y]      = meshgrid(x, y);

% Time Settings
t_final     = 50;               % Total simulation time (s)
dt_common   = 0.1;              % Common time step for ADI, PDE, and evaluation
t_common    = 0:dt_common:t_final;

%% 2. Run Explicit FTCS Solver
disp('1/3: Running Explicit FTCS Solver');
% Explicit Stability Setup
dt_max_cond = 1 / (2 * alpha * (1/dx^2 + 1/dy^2));
dt_max_conv = (rho * cp * dy) / (2 * h_max_limit + 2 * k / dy); 
dt_exp      = 0.9 * min(dt_max_cond, dt_max_conv);

% Pack parameters
params.Nx = Nx; params.Ny = Ny; params.dx = dx; params.dy = dy;
params.alpha = alpha; params.dt = dt_exp; params.Nt = ceil(t_final / dt_exp);
params.k = k; params.q_flux = q_flux; params.h_side = h_side; 
params.T_inf = T_inf; params.T_init = T_init;
params.h_base = h_base; params.T_target = T_target;
params.Kp = Kp; params.h_max_limit = h_max_limit;

tic;
[T_explicit, ~, ~, t_exp_hist, Tmax_exp_hist, htop_exp_hist] = solve_Explicit(params);
time_explicit = toc;
disp(['Explicit completed in ', num2str(time_explicit, '%.2f'), ' seconds.']);

% Interpolate Explicit results to common time grid
Tmax_exp_interp = interp1(t_exp_hist, Tmax_exp_hist, t_common, 'linear', 'extrap');
htop_exp_interp = interp1(t_exp_hist, htop_exp_hist, t_common, 'linear', 'extrap');

%% 3. Run ADI Solver
disp('2/3: Running ADI Solver');
T_ADI         = T_init * ones(Nx, Ny);
Tmax_adi_hist = zeros(1, length(t_common));
htop_adi_hist = zeros(1, length(t_common));

Tmax_adi_hist(1) = T_init;
htop_adi_hist(1) = h_base;

tic;
for n = 1:(length(t_common)-1)
    T_ADI = solve_ADI(T_ADI, Nx, Ny, dx, dy, dt_common, alpha, k, q_flux, h_side, h_base, T_inf, Kp, T_target, h_max_limit);
    T_top_avg = mean(T_ADI(:, Ny));
    h_top_current = min(h_base + Kp * max(0, T_top_avg - T_target), h_max_limit);
    
    Tmax_adi_hist(n+1) = max(T_ADI(:));
    htop_adi_hist(n+1) = h_top_current;
end
time_ADI = toc;
disp(['ADI completed in ', num2str(time_ADI, '%.2f'), ' seconds.']);

%% 4. Run Reference Solver (PDE Toolbox)
disp('3/3: Running MATLAB PDE Toolbox Solver');
thermalmodel = createpde('thermal', 'transient');
geom = decsg([3; 4; 0; Lx; Lx; 0; 0; 0; Ly; Ly]);
geometryFromEdges(thermalmodel, geom);

thermalProperties(thermalmodel, 'ThermalConductivity', k, 'MassDensity', rho, 'SpecificHeat', cp);
thermalBC(thermalmodel, 'Edge', 1, 'HeatFlux', q_flux);
thermalBC(thermalmodel, 'Edge', [2, 4], 'ConvectionCoefficient', h_side, 'AmbientTemperature', T_inf);
thermalBC(thermalmodel, 'Edge', 3, 'ConvectionCoefficient', @(region, state) smartCoolingFunc(region, state, h_base, T_target, Kp, h_max_limit), 'AmbientTemperature', T_inf);
thermalIC(thermalmodel, T_init);
generateMesh(thermalmodel, 'Hmax', min(Lx/40, Ly/20));

tic;
result = solve(thermalmodel, t_common);
time_PDE = toc;
disp(['PDE Toolbox completed in ', num2str(time_PDE, '%.2f'), ' seconds.']);

T_all_PDE     = result.Temperature; 
Tmax_pde_hist = max(T_all_PDE, [], 1);
tol           = 1e-9;
topNodes      = find(abs(thermalmodel.Mesh.Nodes(2,:) - Ly) < tol);
T_top_avg_pde = mean(T_all_PDE(topNodes, :), 1); 
htop_pde_hist = min(h_base + Kp .* max(0, T_top_avg_pde - T_target), h_max_limit);

%% 5. Error Calculation
% Relative Error (%) = |Numerical - Reference| / Reference * 100
Error_Explicit = abs(Tmax_exp_interp - Tmax_pde_hist) ./ Tmax_pde_hist * 100;
Error_ADI      = abs(Tmax_adi_hist - Tmax_pde_hist) ./ Tmax_pde_hist * 100;

%% 6. Comparative Plotting
disp('Generating Comparative Plots');

% Figure 1: Comparison of Max Temperature
figure('Name', 'Comparison: Max Temperature', 'Color', 'w', 'Position', [100, 100, 800, 500]);
plot(t_common, Tmax_pde_hist, 'k-', 'LineWidth', 2); hold on;
plot(t_common, Tmax_exp_interp, 'r--', 'LineWidth', 2);
plot(t_common, Tmax_adi_hist, 'b-.', 'LineWidth', 2);
yline(90, 'k--', 'Safe Limit (90^\circC)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
title('Comparison of Maximum Temperature over Time');
xlabel('Time (s)');
ylabel('Max Temperature (^\circC)');
legend('Reference (PDE)', 'Explicit FTCS', 'ADI Method', 'Location', 'southeast');
grid on;

% Figure 2: Comparison of Smart Fan Behavior (h_top)
figure('Name', 'Comparison: Smart Fan Behavior', 'Color', 'w', 'Position', [150, 150, 800, 500]);
plot(t_common, htop_pde_hist, 'k-', 'LineWidth', 2); hold on;
plot(t_common, htop_exp_interp, 'r--', 'LineWidth', 2);
plot(t_common, htop_adi_hist, 'b-.', 'LineWidth', 2);
title('Comparison of Smart Fan Operation (h_{top})');
xlabel('Time (s)');
ylabel('h_{top} (W/m^2.K)');
legend('Reference (PDE)', 'Explicit FTCS', 'ADI Method', 'Location', 'southeast');
grid on;

% Figure 3: Relative Error of Numerical Solvers
figure('Name', 'Relative Error Analysis', 'Color', 'w', 'Position', [200, 200, 800, 500]);
plot(t_common, Error_Explicit, 'r-', 'LineWidth', 2); hold on;
plot(t_common, Error_ADI, 'b-', 'LineWidth', 2);
title('Relative Error of Numerical Methods vs Reference PDE');
xlabel('Time (s)');
ylabel('Relative Error (%)');
legend('Error Explicit', 'Error ADI', 'Location', 'northeast');
grid on;

% Figure 4: Final Temperature Contours (Side-by-Side)
figure('Name', 'Final Temperature Distributions', 'Color', 'w', 'Position', [100, 300, 1200, 350]);

subplot(1, 3, 1);
contourf(X', Y', T_explicit, 30, 'LineColor', 'none'); colormap('jet'); 
c = colorbar; c.Label.String = 'Temperature (^\circC)';
title('Explicit FTCS'); xlabel('x (m)'); ylabel('y (m)'); axis equal; axis tight;

subplot(1, 3, 2);
contourf(X', Y', T_ADI, 30, 'LineColor', 'none'); colormap('jet'); 
c = colorbar; c.Label.String = 'Temperature (^\circC)';
title('ADI Method'); xlabel('x (m)'); ylabel('y (m)'); axis equal; axis tight;

subplot(1, 3, 3);
pdeplot(thermalmodel, 'XYData', T_all_PDE(:, end), 'Contour', 'on', 'ColorMap', 'jet');
title('Reference PDE Toolbox'); xlabel('x (m)'); ylabel('y (m)'); axis equal; axis tight;

%% Local Function for PDE Toolbox Smart Boundary Condition
function h = smartCoolingFunc(region, state, h_base, T_target, Kp, h_max_limit)
    if any(isnan(state.u))
        h = NaN(1, numel(region.x));
        return;
    end
    h_current = h_base + Kp * max(0, state.u - T_target);
    h = min(h_current, h_max_limit);
end
